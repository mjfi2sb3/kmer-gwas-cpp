#ifndef BIN_STORE_HPP
#define BIN_STORE_HPP

// ===========================================================================
//  BinStore — bounded-memory k-mer accumulation for Stage 1.
// ===========================================================================
//
//  Replaces the old scheme, which for every accession created
//      1 dir + num_bins dirs + 2 x num_bins files   (4,501 inodes at 1500 bins)
//  wrote every hash-map entry into `bin_i/keys.dat` + `bin_i/values.dat`, and
//  then re-read the lot in a second pass to deduplicate. Because the counting
//  chunks were small, a k-mer seen 30x was written to disk 30 times: measured
//  0.13 TB written + 0.13 TB read back per rice accession, 5.8 TB + 5.8 TB for
//  wheat.
//
//  HOW IT WORKS
//
//  Workers append packed Records into per-bin in-memory buffers. Each bin has a
//  share of a global memory budget. When a bin exceeds its share:
//
//     1. COMPACT  - sort the buffer and run-length-merge equal keys, summing
//                   counts. At 30x coverage this reclaims ~97% of the space, so
//                   in the common case nothing ever reaches disk.
//     2. SPILL    - only if the bin is STILL over its share after compaction,
//                   append the compacted block to ONE shared spill file and
//                   record (offset, count). This is a safety valve for genomes
//                   whose distinct k-mer set genuinely exceeds the budget.
//
//  At the end, each bin merges its in-memory remainder with any spilled blocks,
//  applies the count filter, and writes one sorted slice to the pack.
//
//  Peak memory is the budget, by construction, whatever the genome — that is
//  the property this class exists to provide.
//
//  Inodes per accession: 1 pack (+1 spill file, only if spilling happened),
//  down from 4,501.
//
//  WHY SORTING RATHER THAN A HASH MAP: measured at 16.3-16.9 M k-mer/s for
//  per-bin partition+sort+compact versus 9.3-10.1 M/s for unordered_map
//  insertion — sorting is ~1.6x FASTER, because random inserts into a large map
//  are cache- and TLB-miss bound while per-bin sorts stay cache-resident. The
//  compaction itself is ~1% of the time. Sorting also produces exactly the
//  ordering the pack format and Stage 2's k-way merge require.
// ===========================================================================

#include "kmer_key.hpp"
#include "pack_io.hpp"

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <fstream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace kmer {

// Sort a buffer and merge equal keys, summing counts (saturating). Returns the
// compacted length; the buffer is truncated to it.
inline void sort_and_compact(std::vector<Record>& v)
{
    if (v.size() < 2) return;
    std::sort(v.begin(), v.end(),
              [](const Record& a, const Record& b) { return a.key < b.key; });
    size_t w = 0;
    for (size_t r = 0; r < v.size(); ) {
        size_t s = r;
        uint32_t acc = 0;
        while (r < v.size() && v[r].key == v[s].key) { acc += v[r].count; r++; }
        v[w].key   = v[s].key;
        v[w].count = (acc > COUNT_MAX) ? COUNT_MAX : (Count)acc;
        w++;
    }
    v.resize(w);
}

class BinStore {
public:
    // budget_bytes is the TOTAL across all bins; each bin gets an equal share.
    BinStore(size_t num_bins, size_t budget_bytes, std::string spill_path)
        : bins_(num_bins), spill_path_(std::move(spill_path))
    {
        // Keep a sane floor so tiny budgets or huge bin counts still make progress.
        per_bin_limit_ = budget_bytes / (num_bins ? num_bins : 1) / RECORD_BYTES;
        if (per_bin_limit_ < 4096) per_bin_limit_ = 4096;
    }

    ~BinStore() { close_spill(); }

    // Bulk-append a run of records that all belong to `bin`. Appending in runs
    // means one lock acquisition per (worker, bin) rather than one per record.
    void append(size_t bin, const Record* recs, size_t n)
    {
        Bin& b = bins_[bin];
        std::lock_guard<std::mutex> lk(b.m);
        b.buf.insert(b.buf.end(), recs, recs + n);
        if (b.buf.size() >= per_bin_limit_) reduce(bin, b);
    }

    // Finalize `bin` IN PLACE: after this returns, bins_[bin].buf holds the one
    // sorted, deduplicated, count-filtered slice ready to write to the pack.
    //
    // Safe to call concurrently for DIFFERENT bins: each call touches only its
    // own Bin, and the spill merge reads through a per-call ifstream with its own
    // file offset, so no mutable state is shared between bins. This is what lets
    // the finalize phase run across the thread pool instead of one bin at a time;
    // the per-bin sorts dominate this stage, so parallelising them is the win.
    // (The pack WRITE stays serial in the caller — bins must be written in
    // ascending order — but it is cheap I/O next to the sorting.)
    size_t finalize_bin_inplace(size_t bin, Count min_count)
    {
        Bin& b = bins_[bin];
        sort_and_compact(b.buf);

        if (!b.spills.empty()) {
            // Concatenate the spilled blocks with the in-memory remainder and
            // compact once. Per accession a single bin holds only
            // (total / num_bins) records, so this stays small even when the
            // accession as a whole was far too large to hold at once.
            std::ifstream in(spill_path_, std::ios::binary);
            if (!in) throw std::runtime_error("BinStore: cannot reopen spill file " + spill_path_);
            size_t total = b.buf.size();
            for (auto& s : b.spills) total += s.second;
            std::vector<Record> merged;
            merged.reserve(total);
            merged.insert(merged.end(), b.buf.begin(), b.buf.end());
            std::vector<Record> tmp;
            for (auto& s : b.spills) {
                tmp.resize(s.second);
                in.seekg((std::streamoff)s.first, std::ios::beg);
                in.read(reinterpret_cast<char*>(tmp.data()),
                        (std::streamsize)(s.second * RECORD_BYTES));
                if (!in) throw std::runtime_error("BinStore: short read from spill file");
                merged.insert(merged.end(), tmp.begin(), tmp.end());
            }
            sort_and_compact(merged);
            b.buf.swap(merged);
        }

        // Drop k-mers below the error threshold. Filtered in place.
        //
        // Note there is deliberately NO is_zero() check here. The all-zero key
        // is a legitimate k-mer (poly-A: A encodes as 00), and filtering it out
        // was silently discarding it. Invalid windows are now skipped at
        // counting time rather than being funnelled into the zero key, so
        // nothing spurious reaches this point. See task #6.
        size_t w = 0;
        for (size_t i = 0; i < b.buf.size(); i++)
            if (b.buf[i].count >= min_count)
                b.buf[w++] = b.buf[i];
        b.buf.resize(w);
        return b.buf.size();
    }

    // Read-only view of a finalized bin's slice, for writing to the pack.
    const std::vector<Record>& bin_slice(size_t bin) const { return bins_[bin].buf; }

    // Release a bin's memory once it has been written to the pack.
    void release_bin(size_t bin) { std::vector<Record>().swap(bins_[bin].buf); }

    // Serial convenience wrapper (finalize in place, then hand the slice out).
    void finalize_bin(size_t bin, std::vector<Record>& out, Count min_count)
    {
        finalize_bin_inplace(bin, min_count);
        out.swap(bins_[bin].buf);
    }

    bool   spilled()      const { return spill_rounds_ > 0; }
    size_t spill_rounds() const { return spill_rounds_; }

    void close_spill()
    {
        if (spill_out_.is_open()) spill_out_.close();
        if (spill_started_) { std::remove(spill_path_.c_str()); spill_started_ = false; }
    }

private:
    struct Bin {
        std::vector<Record>                          buf;
        std::vector<std::pair<uint64_t, size_t>>     spills;  // (byte offset, record count)
        std::mutex                                   m;
    };

    // Called with b.m held.
    void reduce(size_t bin, Bin& b)
    {
        sort_and_compact(b.buf);
        if (b.buf.size() < per_bin_limit_) return;   // compaction was enough

        // Still over budget: spill this bin's compacted block.
        std::lock_guard<std::mutex> lk(spill_m_);
        if (!spill_started_) {
            spill_out_.open(spill_path_, std::ios::binary | std::ios::trunc);
            if (!spill_out_) throw std::runtime_error("BinStore: cannot create spill file " + spill_path_);
            spill_started_ = true;
        }
        uint64_t off = spill_pos_;
        spill_out_.write(reinterpret_cast<const char*>(b.buf.data()),
                         (std::streamsize)(b.buf.size() * RECORD_BYTES));
        if (!spill_out_) throw std::runtime_error("BinStore: spill write failed");
        spill_out_.flush();
        spill_pos_ += b.buf.size() * RECORD_BYTES;
        b.spills.emplace_back(off, b.buf.size());
        std::vector<Record>().swap(b.buf);
        spill_rounds_++;
        (void)bin;
    }

    std::vector<Bin>   bins_;
    size_t             per_bin_limit_ = 0;   // in RECORDS
    std::string        spill_path_;
    std::ofstream      spill_out_;
    std::mutex         spill_m_;
    uint64_t           spill_pos_     = 0;
    bool               spill_started_ = false;
    std::atomic<size_t> spill_rounds_{0};
};

} // namespace kmer

#endif
