#ifndef PACK_IO_HPP
#define PACK_IO_HPP

// ===========================================================================
//  Pack file — the interface between Stage 1 (kmer_count) and Stage 2
//  (matrix_merge).
// ===========================================================================
//
//  ONE file per accession, holding that accession's k-mers for every bin:
//
//     +------------------+ 0
//     | bin 0 records    |
//     | bin 1 records    |
//     | ...              |
//     | bin B-1 records  |
//     +------------------+ index_offset
//     | (B+1) x uint64   |   byte offset of each bin; [B] is the end sentinel,
//     | offset table     |   so bin i occupies [off[i], off[i+1])
//     +------------------+
//     | PackFooter       |   fixed 24 bytes, magic last => at end of file
//     +------------------+ EOF
//
//  Reading bin i is: read the 24-byte tail, seek to off[i], read
//  off[i+1]-off[i] bytes. No scan, no extraction, no temporary files.
//
//  One pack file per accession holds every bin, so Stage 2 reads one file per
//  accession and seeks within it, rather than materialising a separate file per
//  (accession, bin) pair. The cost that matters here is inodes and metadata
//  operations, not read volume.
//
//  HARD CONTRACT: records within each bin slice are SORTED by key and
//  deduplicated. Stage 2's k-way merge depends on it; validate_pack() checks it.
//
//  SELF-DESCRIBING: the footer records k and the record width. A pack written
//  at k=31 (10-byte records) opened by a binary built for k=51 (18-byte records)
//  is rejected with a clear error instead of decoding to silent nonsense.
// ===========================================================================

#include "kmer_key.hpp"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace kmer {

static constexpr uint16_t PACK_VERSION = 1;

#pragma pack(push, 1)
struct PackFooter {
    uint64_t index_offset;    // byte offset of the offset table
    uint32_t num_bins;        // B
    uint16_t version;         // PACK_VERSION
    uint16_t kmer_k;          // k the records were written at
    uint16_t record_bytes;    // sizeof(Key) + sizeof(Count)
    uint16_t reserved;        // 0; keeps the struct 8-byte friendly
    char     magic[4];        // "KGWS" — last, so it sits at EOF
};
#pragma pack(pop)

static_assert(sizeof(PackFooter) == 24, "PackFooter must be exactly 24 bytes");

inline const char* pack_magic() { return "KGWS"; }

// ---------------------------------------------------------------------------
// PackWriter — bins must be written in ascending order, each exactly once.
// ---------------------------------------------------------------------------
class PackWriter {
public:
    PackWriter(const std::string& path, uint32_t num_bins)
        : path_(path), num_bins_(num_bins)
    {
        out_.open(path, std::ios::binary | std::ios::trunc);
        if (!out_) throw std::runtime_error("PackWriter: cannot create " + path);
        offsets_.reserve(num_bins + 1);
        offsets_.push_back(0);
    }

    // Append one bin's slice. Records must already be sorted and deduplicated —
    // that is the format's contract, not something this class enforces on the
    // hot path (validate_pack() verifies it).
    void write_bin(uint32_t bin, const Record* recs, size_t n)
    {
        if (bin != next_bin_)
            throw std::runtime_error("PackWriter: bins must be written in order, expected "
                                     + std::to_string(next_bin_) + " got " + std::to_string(bin));
        if (n)
            out_.write(reinterpret_cast<const char*>(recs), (std::streamsize)(n * RECORD_BYTES));
        if (!out_) throw std::runtime_error("PackWriter: write failed on " + path_);
        pos_ += n * RECORD_BYTES;
        offsets_.push_back(pos_);
        next_bin_++;
    }

    void write_bin(uint32_t bin, const std::vector<Record>& recs)
    { write_bin(bin, recs.data(), recs.size()); }

    // Write the offset table and footer. Must be called exactly once.
    void finish()
    {
        if (next_bin_ != num_bins_)
            throw std::runtime_error("PackWriter: only " + std::to_string(next_bin_) + " of "
                                     + std::to_string(num_bins_) + " bins written");
        PackFooter f{};
        f.index_offset = pos_;
        f.num_bins     = num_bins_;
        f.version      = PACK_VERSION;
        f.kmer_k       = (uint16_t)K;
        f.record_bytes = (uint16_t)RECORD_BYTES;
        f.reserved     = 0;
        std::memcpy(f.magic, pack_magic(), 4);

        out_.write(reinterpret_cast<const char*>(offsets_.data()),
                   (std::streamsize)(offsets_.size() * sizeof(uint64_t)));
        out_.write(reinterpret_cast<const char*>(&f), sizeof(f));
        out_.flush();
        if (!out_) throw std::runtime_error("PackWriter: failed to finalise " + path_);
        out_.close();
        finished_ = true;
    }

    ~PackWriter() { if (!finished_ && out_.is_open()) out_.close(); }

private:
    std::string           path_;
    std::ofstream         out_;
    uint32_t              num_bins_;
    uint32_t              next_bin_ = 0;
    uint64_t              pos_      = 0;
    std::vector<uint64_t> offsets_;
    bool                  finished_ = false;
};

// ---------------------------------------------------------------------------
// PackReader — random access to any bin by seek.
// ---------------------------------------------------------------------------
class PackReader {
public:
    explicit PackReader(const std::string& path) : path_(path)
    {
        in_.open(path, std::ios::binary);
        if (!in_) throw std::runtime_error("PackReader: cannot open " + path);

        in_.seekg(0, std::ios::end);
        std::streamoff size = in_.tellg();
        if (size < (std::streamoff)sizeof(PackFooter))
            throw std::runtime_error("PackReader: " + path + " is too small to be a pack file");

        in_.seekg(size - (std::streamoff)sizeof(PackFooter), std::ios::beg);
        in_.read(reinterpret_cast<char*>(&f_), sizeof(f_));
        if (!in_) throw std::runtime_error("PackReader: cannot read footer of " + path);

        if (std::memcmp(f_.magic, pack_magic(), 4) != 0)
            throw std::runtime_error("PackReader: " + path + " is not a pack file (bad magic)");
        if (f_.version != PACK_VERSION)
            throw std::runtime_error("PackReader: " + path + " has pack version "
                                     + std::to_string(f_.version) + ", this build expects "
                                     + std::to_string(PACK_VERSION));
        // The check that makes a k mismatch loud instead of silent.
        if (f_.kmer_k != (uint16_t)K)
            throw std::runtime_error("PackReader: " + path + " was written with k="
                                     + std::to_string(f_.kmer_k) + " but this binary was built for k="
                                     + std::to_string(K) + " — rebuild with -DKMER_K="
                                     + std::to_string(f_.kmer_k));
        if (f_.record_bytes != (uint16_t)RECORD_BYTES)
            throw std::runtime_error("PackReader: " + path + " has "
                                     + std::to_string(f_.record_bytes)
                                     + "-byte records, this build expects "
                                     + std::to_string(RECORD_BYTES));

        offsets_.resize(f_.num_bins + 1);
        in_.seekg((std::streamoff)f_.index_offset, std::ios::beg);
        in_.read(reinterpret_cast<char*>(offsets_.data()),
                 (std::streamsize)(offsets_.size() * sizeof(uint64_t)));
        if (!in_) throw std::runtime_error("PackReader: cannot read offset table of " + path);

        for (size_t i = 1; i < offsets_.size(); i++)
            if (offsets_[i] < offsets_[i-1])
                throw std::runtime_error("PackReader: corrupt (non-monotonic) offset table in " + path);
        if (offsets_.back() != f_.index_offset)
            throw std::runtime_error("PackReader: offset table end does not match data size in " + path);
    }

    uint32_t num_bins() const { return f_.num_bins; }
    int      kmer_k()   const { return f_.kmer_k; }

    size_t bin_bytes(uint32_t bin) const
    {
        check_bin(bin);
        return (size_t)(offsets_[bin + 1] - offsets_[bin]);
    }
    size_t bin_records(uint32_t bin) const { return bin_bytes(bin) / RECORD_BYTES; }

    uint64_t bin_offset(uint32_t bin) const { check_bin(bin); return offsets_[bin]; }

    // Read one bin's slice. This is a seek plus a single sequential read.
    void read_bin(uint32_t bin, std::vector<Record>& out)
    {
        check_bin(bin);
        size_t n = bin_records(bin);
        out.resize(n);
        if (!n) return;
        in_.seekg((std::streamoff)offsets_[bin], std::ios::beg);
        in_.read(reinterpret_cast<char*>(out.data()), (std::streamsize)(n * RECORD_BYTES));
        if (!in_) throw std::runtime_error("PackReader: short read on bin "
                                           + std::to_string(bin) + " of " + path_);
    }

private:
    void check_bin(uint32_t bin) const
    {
        if (bin >= f_.num_bins)
            throw std::runtime_error("PackReader: bin " + std::to_string(bin)
                                     + " out of range (pack has " + std::to_string(f_.num_bins)
                                     + ") in " + path_);
    }

    std::string           path_;
    mutable std::ifstream in_;
    PackFooter            f_{};
    std::vector<uint64_t> offsets_;
};

// ---------------------------------------------------------------------------
// BinCursor — streams ONE bin's slice out of a pack, a buffer at a time.
//
// This is what makes Stage 2 memory O(number of accessions) rather than
// O(union k-mers x accessions): the merge holds one small buffer per accession
// instead of the whole bin's union in a hash map.
//
// The cursor keeps its file handle open for its lifetime, so a merge over N
// accessions needs N descriptors. That is deliberate — the alternative
// (reopening on every refill) turns a bin task into millions of open() calls.
// Callers should raise RLIMIT_NOFILE accordingly.
// ---------------------------------------------------------------------------
class BinCursor {
public:
    BinCursor(const std::string& path, uint32_t bin, size_t buf_records = 4096)
        : path_(path), buf_cap_(buf_records ? buf_records : 1)
    {
        PackReader r(path);                       // validates magic/version/k
        if (bin >= r.num_bins())
            throw std::runtime_error("BinCursor: bin " + std::to_string(bin)
                                     + " out of range in " + path);
        remaining_ = r.bin_records(bin);
        data_off_  = r.bin_offset(bin);

        in_.open(path, std::ios::binary);
        if (!in_) throw std::runtime_error("BinCursor: cannot open " + path);
        in_.seekg((std::streamoff)data_off_, std::ios::beg);
        buf_.resize(buf_cap_);
        refill();
    }

    bool valid() const { return pos_ < filled_; }
    const Record& current() const { return buf_[pos_]; }

    void advance()
    {
        if (++pos_ >= filled_) refill();
    }

private:
    void refill()
    {
        pos_ = 0;
        filled_ = 0;
        if (!remaining_) return;
        size_t n = remaining_ < buf_cap_ ? remaining_ : buf_cap_;
        in_.read(reinterpret_cast<char*>(buf_.data()), (std::streamsize)(n * RECORD_BYTES));
        if (!in_) throw std::runtime_error("BinCursor: short read from " + path_);
        filled_     = n;
        remaining_ -= n;
    }

    std::string         path_;
    std::ifstream       in_;
    std::vector<Record> buf_;
    size_t              buf_cap_   = 0;
    size_t              pos_       = 0;
    size_t              filled_    = 0;
    size_t              remaining_ = 0;
    uint64_t            data_off_  = 0;
};

// ---------------------------------------------------------------------------
// validate_pack — structural check plus the sorted/deduplicated contract.
// Returns the total record count; throws with a specific message on failure.
// ---------------------------------------------------------------------------
inline size_t validate_pack(const std::string& path)
{
    PackReader r(path);
    std::vector<Record> buf;
    size_t total = 0;
    for (uint32_t b = 0; b < r.num_bins(); b++) {
        if (r.bin_bytes(b) % RECORD_BYTES != 0)
            throw std::runtime_error("validate_pack: bin " + std::to_string(b) + " of " + path
                                     + " is not a whole number of records");
        r.read_bin(b, buf);
        for (size_t i = 1; i < buf.size(); i++) {
            if (buf[i].key < buf[i-1].key)
                throw std::runtime_error("validate_pack: bin " + std::to_string(b) + " of " + path
                                         + " is not sorted at record " + std::to_string(i));
            if (buf[i].key == buf[i-1].key)
                throw std::runtime_error("validate_pack: bin " + std::to_string(b) + " of " + path
                                         + " has a duplicate key at record " + std::to_string(i));
        }
        total += buf.size();
    }
    return total;
}

} // namespace kmer

#endif
