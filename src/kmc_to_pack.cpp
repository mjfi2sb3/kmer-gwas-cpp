// ===========================================================================
//  kmc_to_pack — convert a KMC database into a Stage-1 .kbin pack.
// ===========================================================================
//  Licence: GNU GPL 3 (this file links the KMC API in src/kmc_api/; see src/kmc_api/README).
//
//  Reads the KMC database directly in binary through the KMC C++ API
//  (CKMCFile::OpenForListing + ReadNextKmer), so there is no text conversion:
//  each k-mer arrives in its native 2-bit form and is reshaped straight into a Key.
//
//  The pack data is byte-compatible with what kmer_count emits from reads: same
//  canonical key, same `legacy_bin_hash(key) % num_bins` binning, the same
//  16-bit count saturation at 65535, and the same sorted/deduplicated per-bin
//  slices. A matrix built from these packs equals one built from the fastq path
//  when the database and the reads were filtered at the same minimum count.
//
//  COUNT FILTER. By default the database is taken as it is: every k-mer in it is
//  kept, and the database's own -ci is the filter. A min_count argument larger
//  than the database's -ci enforces that filter in addition. A min_count equal to
//  or below -ci needs no filter: KMC already removed those k-mers. (Applying it
//  anyway would empty a database built with -cs1 or `kmc_tools compact`, where
//  every count is stored as 1.) If an enforced filter would remove every k-mer,
//  the job stops: the database's counts never reach min_count.
//
//  Every k-mer is RE-CANONICALISED here (canonical = min(forward, revcomp)),
//  built directly from the API's 2-bit symbols. This is idempotent for a
//  canonical database (KMC default) and correct for a non-canonical one
//  (`kmc -b`): a k-mer and its reverse complement fold to the same key and their
//  counts are summed. k is compiled in via -DKMER_K and MUST match the
//  database's k (checked against CKMCFileInfo; a mismatch aborts).
//
//  With --headers, the program reads only the database headers of a list of
//  accessions, writes them as a table and checks that the databases were built
//  with the same settings (see headers_main).
//
//  LICENCE: this program links the KMC C++ API (src/kmc_api/, GNU GPL 3), so the
//  resulting kmc_to_pack binary is GPL. It reuses the pipeline's own pack/key
//  code; nothing else in the pipeline links KMC.
//
//  Build (k baked in, KMC API compiled in, no external KMC install needed):
//      g++ -std=c++17 -O3 -march=native -pthread -DKMER_K=<k> -o kmc_to_pack
//          kmc_to_pack.cpp kmc_api/kmc_file.cpp kmc_api/kmer_api.cpp kmc_api/mmer.cpp
// ===========================================================================

#include "kmc_api/kmc_file.h"
#include "kmc_api/kmer_api.h"

#include "kmer_key.hpp"
#include "pack_io.hpp"
#include "bin_store.hpp"
#include "mem_limit.hpp"
#include "progress.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

using namespace kmer;

// Lower bound on the accumulation budget. Overridable at compile time only so the
// spill path can be exercised in tests; production always uses the 1 GB default.
#ifndef KMC2PACK_BUDGET_FLOOR_BYTES
#define KMC2PACK_BUDGET_FLOOR_BYTES (1ULL << 30)
#endif

namespace {

[[noreturn]] void die(const std::string& msg) {
    std::cerr << "kmc_to_pack: ERROR: " << msg << "\n";
    std::exit(1);
}
void warn(const std::string& msg) { std::cerr << "kmc_to_pack: WARNING: " << msg << "\n"; }
void note(const std::string& msg) { std::cerr << "kmc_to_pack: " << msg << "\n"; }

// Checks that `path` exists and starts and ends with the 4-byte `marker`.
// Returns an error message, or an empty string if the file is valid.
std::string kmc_file_error(const std::string& path, const char* marker, uint64_t* size_out = nullptr) {
    std::error_code ec;
    if (!std::filesystem::is_regular_file(path, ec)) return "missing file " + path;
    const uint64_t size = std::filesystem::file_size(path, ec);
    char head[4] = {}, tail[4] = {};
    std::ifstream f(path, std::ios::binary);
    if (!ec && size >= 8 && f.read(head, 4) && f.seekg(-4, std::ios::end) && f.read(tail, 4)
        && std::memcmp(head, marker, 4) == 0 && std::memcmp(tail, marker, 4) == 0) {
        if (size_out) *size_out = size;
        return std::string();
    }
    return path + " is not a KMC database file (no " + std::string(marker, 4) + " markers)";
}

// Reads the database format version from a .kmc_pre file (the uint32 that
// precedes the header offset and the trailing marker): 0 = KMC1 format,
// 0x200 = KMC2 format. Returns an error message for any other format.
std::string kmc_format_error(const std::string& pre_path, uint32_t& format) {
    uint64_t size = 0;
    const std::string err = kmc_file_error(pre_path, "KMCP", &size);
    if (!err.empty()) return err;
    format = 0;
    std::ifstream f(pre_path, std::ios::binary);
    if (size < 16 || !f.seekg(-12, std::ios::end)
        || !f.read(reinterpret_cast<char*>(&format), sizeof(format)))
        return pre_path + " is not a KMC database file (too short)";
    if (format != 0 && format != 0x200) {
        std::ostringstream hex;
        hex << "0x" << std::hex << format;
        return "unsupported KMC database format " + hex.str() + " in " + pre_path
               + "; supported: 0 (KMC1), 0x200 (KMC2)";
    }
    return std::string();
}

// Checks both files of a database and reads its header through the API.
// Returns an error message, or an empty string on success.
std::string open_kmc_header(const std::string& db_prefix, uint32_t& format, CKMCFileInfo& info) {
    std::string err = kmc_format_error(db_prefix + ".kmc_pre", format);
    if (err.empty()) err = kmc_file_error(db_prefix + ".kmc_suf", "KMCS");
    if (!err.empty()) return err;
    CKMCFile db;
    if (!db.OpenForListing(db_prefix))
        return "cannot open KMC database " + db_prefix;
    const bool ok = db.Info(info);
    db.Close();
    return ok ? std::string() : "cannot read the KMC database header of " + db_prefix;
}

// The KMC API skips every k-mer whose stored count is below the header's -ci,
// and SetMinCount cannot go below it. In a database built with -cs1 or made by
// `kmc_tools compact`, every stored count is 1 although each k-mer passed -ci,
// so the API would return nothing. This reader lowers the API's own filter to 1;
// for any other database it changes nothing, because every stored count is
// already at least -ci.
class KmcReader : public CKMCFile {
public:
    void read_all_counts() { min_count = 1; }
};

const char* format_name(uint32_t format) { return format == 0x200 ? "KMC2" : "KMC1"; }

void usage(const char* prog) {
    std::cerr
        << "usage: " << prog << " <kmc_db_prefix> <num_bins> <output_dir> [out_prefix] [min_count] [budget_gb]\n"
        << "       " << prog << " --headers <accessions_file> <kmc_dir> <out_tsv> [min_count] [keep_counts y|n] [allow_mixed y|n]\n"
        << "\n"
        << "  <kmc_db_prefix> : path to a KMC database, WITHOUT the .kmc_pre/.kmc_suf suffix\n"
        << "  writes <output_dir>/<out_prefix>.kbin and <out_prefix>.pack_stats.tsv\n"
        << "  (out_prefix defaults to the db basename)\n"
        << "  k is compiled in (-DKMER_K=" << K << ") and MUST match the database's k.\n"
        << "  min_count : 0 (default) keeps every k-mer in the database; its own -ci applies.\n"
        << "              A value above the database's -ci is enforced in addition.\n"
        << "  budget_gb : accumulation budget; 0 = auto (default). Affects spilling only.\n"
        << "\n"
        << "  --headers reads the header of <kmc_dir>/<accession> for every accession,\n"
        << "  writes them to <out_tsv>, and stops if the databases cannot be combined\n"
        << "  (different k, -ci, -cx, or counter size when counts are kept).\n"
        << "  allow_mixed y reports different -ci, -cx or counter sizes as warnings instead.\n"
        << "\n"
        << "  A progress line is written every 600 s. The environment variable\n"
        << "  KMER_GWAS_PROGRESS_SEC changes the interval; 0 turns it off.\n";
}

// ---------------------------------------------------------------------------
// --headers: read every database header, write a table, and check that the
// databases can be combined in one matrix.
// ---------------------------------------------------------------------------
int headers_main(int argc, char** argv) {
    if (argc < 5) { usage(argv[0]); return 1; }
    const std::string acc_file    = argv[2];
    const std::string kmc_dir     = argv[3];
    const std::string out_tsv     = argv[4];
    const long        min_count   = (argc > 5) ? std::strtol(argv[5], nullptr, 10) : 0;
    const bool        keep_counts = (argc > 6) && std::string(argv[6]) == "y";
    const bool        allow_mixed = (argc > 7) && std::string(argv[7]) == "y";

    std::ifstream in(acc_file);
    if (!in) die("cannot read accessions file " + acc_file);
    std::vector<std::string> accessions;
    for (std::string line; std::getline(in, line); ) {
        const size_t a = line.find_first_not_of(" \t\r");
        const size_t b = line.find_last_not_of(" \t\r");
        if (a != std::string::npos) accessions.push_back(line.substr(a, b - a + 1));
    }
    if (accessions.empty()) die("no accessions in " + acc_file);

    std::ofstream out(out_tsv);
    out << "accession\tformat\tk\tci\tcx\tcounter_bytes\tcanonical\tdistinct_kmers\n";

    // value -> number of databases, for each setting that must agree
    std::map<std::string, size_t> by_k, by_ci, by_cx, by_counter;
    std::vector<std::string> errors;
    uint32_t max_ci = 0;
    for (const auto& acc : accessions) {
        const std::string prefix = kmc_dir + "/" + acc;
        uint32_t format = 0;
        CKMCFileInfo info{};
        const std::string err = open_kmc_header(prefix, format, info);
        if (!err.empty()) { errors.push_back(acc + ": " + err); continue; }
        out << acc << '\t' << format_name(format) << '\t' << info.kmer_length << '\t'
            << info.min_count << '\t' << info.max_count << '\t' << info.counter_size << '\t'
            << (info.both_strands ? "yes" : "no (-b)") << '\t' << info.total_kmers << '\n';
        by_k[std::to_string(info.kmer_length)]++;
        by_ci[std::to_string(info.min_count)]++;
        by_cx[std::to_string(info.max_count)]++;
        by_counter[std::to_string(info.counter_size)]++;
        max_ci = std::max(max_ci, (uint32_t)info.min_count);
    }
    out.close();
    if (!out) die("cannot write " + out_tsv);

    // "-ci1: 2 databases; -ci2: 1 database"
    auto groups = [](const std::string& label, const std::map<std::string, size_t>& m) {
        std::string s;
        for (auto& kv : m)
            s += (s.empty() ? "" : "; ") + label + kv.first + ": " + std::to_string(kv.second)
                 + (kv.second == 1 ? " database" : " databases");
        return s;
    };

    bool stop = false;
    for (auto& e : errors) { std::cerr << "kmc_to_pack: ERROR: " << e << "\n"; stop = true; }
    if (by_k.size() > 1) {
        std::cerr << "kmc_to_pack: ERROR: the databases have different k-mer lengths ("
                  << groups("k=", by_k) << ")\n";
        stop = true;
    } else if (!by_k.empty() && by_k.begin()->first != std::to_string(K)) {
        std::cerr << "kmc_to_pack: ERROR: the databases have k=" << by_k.begin()->first
                  << " but the pipeline runs at k=" << K << "\n";
        stop = true;
    }

    // A setting that differs stops the run unless allow_mixed is given.
    auto mixed = [&](const std::string& what, const std::string& label,
                     const std::map<std::string, size_t>& m, const std::string& extra) {
        if (m.size() <= 1) return;
        const std::string msg = "the databases were built with different " + what + " ("
                                + groups(label, m) + "; see " + out_tsv + ")" + extra;
        if (allow_mixed) std::cout << "WARNING: " << msg << "\n";
        else { std::cerr << "kmc_to_pack: ERROR: " << msg
                         << ". Rebuild them with the same settings, or allow it with --kmc_allow_mixed true.\n";
               stop = true; }
    };
    if (by_ci.size() > 1 && min_count > 0 && (uint32_t)min_count >= max_ci)
        std::cout << "NOTE: the databases were built with different -ci (" << groups("-ci", by_ci)
                  << "); --min_kmer_count " << min_count << " is applied to all of them, so they are consistent\n";
    else
        mixed("-ci", "-ci", by_ci, min_count > 0 ? std::string() : std::string(
              "; a common --min_kmer_count at least as large as the largest -ci makes them consistent"));
    mixed("-cx", "-cx", by_cx, std::string());
    // Counter sizes matter only when the counts go into the matrix.
    if (keep_counts) mixed("counter sizes", "counter bytes ", by_counter, std::string());

    std::cerr << "kmc_to_pack: read " << accessions.size() - errors.size() << " of "
              << accessions.size() << " database headers into " << out_tsv << "\n";
    return stop ? 1 : 0;
}

// ---------------------------------------------------------------------------
// Conversion of one database.
// ---------------------------------------------------------------------------
int convert_main(int argc, char** argv) {
    if (argc < 4) { usage(argv[0]); return 1; }

    const std::string db_prefix  = argv[1];
    const long long   num_bins_ll = std::strtoll(argv[2], nullptr, 10);
    std::string       output_dir = argv[3];
    std::string       out_prefix = (argc > 4) ? argv[4] : std::filesystem::path(db_prefix).filename().string();
    const long        min_count_in = (argc > 5) ? std::strtol(argv[5], nullptr, 10) : 0;
    const double      budget_gb    = (argc > 6) ? std::strtod(argv[6], nullptr) : 0.0;

    if (num_bins_ll <= 0) die("num_bins must be a positive integer (got '" + std::string(argv[2]) + "')");
    if (min_count_in < 0 || min_count_in > COUNT_MAX)
        die("min_count must be between 0 and " + std::to_string(COUNT_MAX));
    const size_t num_bins = (size_t)num_bins_ll;
    if (out_prefix.empty()) out_prefix = "out";

    if (output_dir.empty())        output_dir = ".";
    if (output_dir.back() != '/')  output_dir += '/';
    std::error_code ec;
    std::filesystem::create_directories(output_dir, ec);

    // --- Open the KMC database and read its header ---------------------------
    // The API reports any open failure the same way, so check the files first
    // to give a precise error.
    uint32_t format = 0;
    {
        std::string err = kmc_format_error(db_prefix + ".kmc_pre", format);
        if (err.empty()) err = kmc_file_error(db_prefix + ".kmc_suf", "KMCS");
        if (!err.empty()) die(err);
    }

    KmcReader db;
    if (!db.OpenForListing(db_prefix))
        die("cannot open KMC database '" + db_prefix + "' (expected " + db_prefix
            + ".kmc_pre and " + db_prefix + ".kmc_suf)");

    CKMCFileInfo info;
    if (!db.Info(info)) die("cannot read KMC database header for '" + db_prefix + "'");

    std::cerr << "kmc_to_pack: db=" << db_prefix
              << "  format=" << format_name(format)
              << "  k=" << info.kmer_length
              << "  min_count(-ci)=" << info.min_count
              << "  max_count=" << info.max_count
              << "  counter_size=" << info.counter_size
              << "  both_strands=" << (info.both_strands ? "true" : "false")
              << "  distinct_kmers=" << info.total_kmers << "\n";

    // k mismatch is fatal: the k-mers live in a different space entirely.
    if ((int)info.kmer_length != K)
        die("database k=" + std::to_string(info.kmer_length)
            + " but this binary was built for k=" + std::to_string(K)
            + " — rebuild with -DKMER_K=" + std::to_string(info.kmer_length));

    db.read_all_counts();

    // --- Count filter ----------------------------------------------------------
    const Count db_ci   = (Count)std::max<uint32_t>(1, std::min<uint32_t>(info.min_count, COUNT_MAX));
    const bool  enforce = min_count_in > (long)db_ci;
    const Count filter  = enforce ? (Count)min_count_in : (Count)1;
    if (enforce)
        note("min_count filter: applied (" + std::to_string(min_count_in) + ", above the database -ci "
             + std::to_string(db_ci) + ")");
    else if (min_count_in > 0)
        note("min_count filter: skipped (database -ci " + std::to_string(db_ci) + " >= "
             + std::to_string(min_count_in) + ")");
    else
        note("min_count filter: none (the database -ci " + std::to_string(db_ci) + " applies)");
    if (!enforce && db_ci == 1)
        note("the database keeps k-mers seen once (-ci1); they are kept in the pack");
    if (info.counter_size == 0)
        note("the database stores no counts; every k-mer has count 1");
    else if (info.counter_size == 1)
        note("the database stores 1-byte counts; counts are capped at 255 or less (-cs)");
    // The largest count the database can hold: 1 with no counts, 255 or less
    // with a 1-byte counter (KMC does not store the exact -cs).
    const Count count_limit = info.counter_size == 0 ? (Count)1
                            : info.counter_size == 1 ? (Count)255 : COUNT_MAX;

    // both_strands == true means KMC ran WITHOUT -b (canonical, the default);
    // false means -b (non-canonical). Re-canonicalisation handles both, but with
    // -b the DB's own -ci filtered per orientation, so split low-count k-mers can
    // be missing.
    if (!info.both_strands) {
        warn("database was built non-canonical (kmc -b); k-mers are re-canonicalised here");
        if (info.min_count > 1)
            warn("non-canonical database with -ci" + std::to_string(info.min_count)
                 + ": a canonical k-mer whose occurrences split across strands may have been "
                   "dropped below the threshold and cannot be recovered");
    }

    // --- Budget (spill threshold only; never affects output) -----------------
    const MemLimit ml = detect_memory_limit();
    size_t budget_bytes;
    std::string budget_note;
    if (budget_gb > 0.0) {
        budget_bytes = (size_t)(budget_gb * (double)(1ULL << 30));
        budget_note  = "explicit";
    } else if (ml.bytes) {
        budget_bytes = (size_t)((double)ml.bytes * 0.7);
        budget_note  = "auto (70% of enforced limit)";
    } else {
        budget_bytes = (size_t)8 * (1ULL << 30);
        budget_note  = "no limit found, default";
    }
    if (budget_bytes < KMC2PACK_BUDGET_FLOOR_BYTES) {
        budget_bytes = KMC2PACK_BUDGET_FLOOR_BYTES;
        budget_note += ", raised to the minimum";
    }
    const double GiB = (double)(1ULL << 30);
    std::cerr << "kmc_to_pack: memory: enforced limit = "
              << (ml.bytes ? std::to_string((long)(ml.bytes / GiB)) + " GB" : std::string("unknown"))
              << " [" << ml.source << "]\n"
              << "kmc_to_pack: memory: budget = " << (long)(budget_bytes / GiB)
              << " GB (" << budget_note << ")\n";

    const std::string pack_path  = output_dir + out_prefix + ".kbin";
    const std::string stats_path = output_dir + out_prefix + ".pack_stats.tsv";
    const std::string spill_path = output_dir + out_prefix + ".kmc2pack.spill";

    BinStore store(num_bins, budget_bytes, spill_path);

    // Per-bin staging, so each append covers a contiguous run (one lock per run).
    const size_t FLUSH = 1024;
    std::vector<std::vector<Record>> pending(num_bins);
    for (auto& v : pending) v.reserve(FLUSH);

    // --- Stream the database, building the canonical key from 2-bit symbols ---
    const uint32 kl = info.kmer_length;
    CKmerAPI kmer(kl);
    uint64 count = 0;
    uint64 max_in = 0;
    size_t n_in = 0;
    // Counters for the progress line only. The read count is published every
    // 2^20 k-mers, so the loop does not pay for an atomic per k-mer.
    std::atomic<size_t> seen{0}, bins_done{0};
    const uint64 total = info.total_kmers;
    ProgressReporter progress(std::cerr, "kmc_to_pack: ", [&](double) {
        const double n = (double)seen.load(std::memory_order_relaxed);
        char pct[16];
        std::snprintf(pct, sizeof pct, "%.0f%%", total ? 100.0 * n / (double)total : 0.0);
        return "reading: " + ProgressReporter::count_text(n) + " of "
               + ProgressReporter::count_text((double)total) + " k-mers (" + pct + "), "
               + std::to_string(store.spill_rounds()) + " spills";
    });
    while (db.ReadNextKmer(kmer, count)) {
        Key f, r;
        for (uint32 i = 0; i < kl; i++) {
            const uint8_t v = kmer.get_num_symbol(i);   // 0..3, A=0 C=1 G=2 T=3
            f.push_fwd(v);
            r.push_rev(v);
        }
        const Key key = (r < f) ? r : f;                 // canonical = min(fwd, revcomp)
        const Count cnt = (count > (uint64)COUNT_MAX) ? COUNT_MAX : (Count)count;
        if (count > max_in) max_in = count;

        const uint32_t bin = (uint32_t)(legacy_bin_hash(key) % num_bins);
        pending[bin].push_back(Record{key, cnt});
        if (pending[bin].size() >= FLUSH) {
            store.append(bin, pending[bin].data(), pending[bin].size());
            pending[bin].clear();
        }
        n_in++;
        if ((n_in & 0xFFFFF) == 0) seen.store(n_in, std::memory_order_relaxed);
    }
    db.Close();

    if (n_in != info.total_kmers)
        die("read " + std::to_string(n_in) + " k-mers but the database header lists "
            + std::to_string(info.total_kmers) + "; the database may be damaged");

    // An enforced filter that no count can pass would leave an empty pack. For
    // a canonical database this is known before sorting; for a -b database the
    // counts of both strands are summed first, so it is checked after writing.
    const std::string never_reached =
        "the counts in this database never reach --min_kmer_count " + std::to_string(min_count_in)
        + " (largest count " + std::to_string(max_in) + "); it was probably built with -cs below "
        + std::to_string(min_count_in) + " or with kmc_tools compact. Rebuild it with a larger -cs, "
        "or omit --min_kmer_count to keep the database as it is";
    if (enforce && n_in > 0 && info.both_strands && max_in < (uint64)min_count_in) die(never_reached);

    // Flush staging remainders.
    for (size_t b = 0; b < num_bins; b++)
        if (!pending[b].empty()) {
            store.append(b, pending[b].data(), pending[b].size());
            std::vector<Record>().swap(pending[b]);
        }

    // Finalize + write in ascending bin order (the pack format requires it).
    PackWriter pack(pack_path, (uint32_t)num_bins);
    PackStats stats(count_limit);
    uint64_t distinct_before = 0;
    progress.set_line([&](double) {
        return "finalizing: " + std::to_string(bins_done.load()) + " of " + std::to_string(num_bins) + " bins";
    });
    for (size_t b = 0; b < num_bins; b++) {
        size_t before = 0;
        store.finalize_bin_inplace(b, filter, &before);
        distinct_before += before;
        const std::vector<Record>& slice = store.bin_slice(b);
        pack.write_bin((uint32_t)b, slice);
        stats.add_bin(slice.data(), slice.size());
        store.release_bin(b);
        bins_done++;
    }
    progress.stop();

    PackMeta meta = pack_meta_base(out_prefix, "kmc", "kmc_to_pack");
    meta.set("kmc_api_version", KMC_VER);
    meta.set("num_bins", (uint64_t)num_bins);
    meta.set("min_kmer_count", (uint64_t)(enforce ? filter : db_ci));
    meta.set("min_kmer_count_source", enforce ? "enforced" : "database");
    stats.to_meta(meta);
    meta.set("distinct_kmers_before_filter", distinct_before);
    pack.finish(meta.text());

    const bool   spilled = store.spilled();
    const size_t rounds  = store.spill_rounds();
    store.close_spill();

    if (enforce && n_in > 0 && stats.distinct == 0) {
        std::filesystem::remove(pack_path, ec);
        die(never_reached);
    }
    write_meta_tsv(stats_path, meta);

    std::cerr << "kmc_to_pack: wrote " << pack_path
              << "  (in=" << n_in << " kmers, out=" << stats.distinct
              << " records, min_kmer_count=" << meta.get("min_kmer_count")
              << ", bins=" << num_bins
              << (spilled ? (", spilled x" + std::to_string(rounds)) : std::string(", no spill"))
              << ")\n";
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    if (argc >= 2 && std::string(argv[1]) == "--headers") return headers_main(argc, argv);
    return convert_main(argc, argv);
}
