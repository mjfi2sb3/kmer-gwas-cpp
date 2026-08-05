// ===========================================================================
//  kbin_dump — inspect or export a .kbin pack file.
// ===========================================================================
//
//  A .kbin pack is one accession's k-mers for every bin, produced by
//  kmer_count. This tool reads a pack directly from its self-describing footer,
//  so ONE binary reads a pack of ANY k (it does not need to be built with a
//  matching -DKMER_K). It either prints a summary or writes the k-mers as text,
//  one "<kmer><TAB><count>" line per record.
//
//  Usage:
//    kbin_dump <file.kbin>                 all k-mers, every bin, to stdout
//    kbin_dump <file.kbin> --info          header summary only (k, bins, totals)
//    kbin_dump <file.kbin> --bin N         only bin N, to stdout
//    kbin_dump <file.kbin> --all_bins DIR  one gzipped file per bin (parallel)
//    kbin_dump <file.kbin> --all_bins DIR --no-gzip   plain-text bin files
//
//  --all_bins writes DIR/<accession>/<N>.tsv.gz (or <N>.tsv with --no-gzip) for
//  every bin, where <accession> is the pack's .kbin prefix, so several accessions
//  can share one DIR without colliding. The bins are independent, so the decode
//  runs across --threads workers (all cores by default), and each worker also
//  gzips its own bin, so both decode and compression scale with cores in one
//  process, with no external pigz pipe.
//
//  The pack layout (see src/pack_io.hpp) is: bin records back to back, then a
//  (num_bins + 1) x uint64 offset table, then a 24-byte footer ending in the
//  magic "KGWS". Each record is a 2-bit-per-base packed key (8 bytes for k <= 32,
//  16 bytes otherwise) followed by a uint16 count. The first base of the k-mer
//  sits in the highest bits, base values A=0 C=1 G=2 T=3, matching decode() in
//  src/kmer_key.hpp. Words are stored little-endian.
// ===========================================================================

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include <zlib.h>

#pragma pack(push, 1)
struct Footer {
    uint64_t index_offset;
    uint32_t num_bins;
    uint16_t version;
    uint16_t kmer_k;
    uint16_t record_bytes;
    uint16_t reserved;
    char     magic[4];
};
#pragma pack(pop)

static void die(const std::string& m) {
    std::cerr << "kbin_dump: " << m << "\n";
    std::exit(1);
}

static void usage(std::ostream& o, const char* prog) {
    o << "usage: " << prog << " <file.kbin> [--info] [--bin N]\n"
      << "                 [--all_bins DIR [--no-gzip] [--threads T]]\n"
      << "\n"
      << "  Reads a .kbin pack (kmer_count's Stage 1 output) using its self-describing\n"
      << "  footer, so one binary reads a pack of any k. With no option it writes every\n"
      << "  k-mer as '<kmer><TAB><count>' to stdout. The pack must be UNCOMPRESSED (it is\n"
      << "  seeked by the footer index); decompress a .kbin.gz first with 'pigz -dc'.\n"
      << "\n"
      << "    --info          print a header summary (k, bins, record width, totals) and exit\n"
      << "    --bin N         dump only bin N to stdout\n"
      << "    --all_bins DIR  write one file per bin into DIR/<accession>/, decoded in\n"
      << "                    parallel; each bin is gzipped by default (<N>.tsv.gz)\n"
      << "    --no-gzip       with --all_bins, write plain text (<N>.tsv)\n"
      << "    --threads T     workers for --all_bins (default: all cores)\n"
      << "\n"
      << "  For a single stdout dump, pipe to gzip or pigz to compress. For a whole pack,\n"
      << "  --all_bins compresses in parallel in-process (no pigz pipe needed).\n";
}

// Append a uint16 count as decimal, no allocation.
static inline void append_count(std::string& s, uint16_t v) {
    char tmp[5];
    int i = 5;
    if (v == 0) { s += '0'; return; }
    while (v) { tmp[--i] = char('0' + v % 10); v /= 10; }
    s.append(tmp + i, 5 - i);
}

int main(int argc, char** argv) {
    std::string path, all_bins_dir;
    bool info = false, all_bins = false, gzip = true;   // --all_bins gzips by default
    long only_bin = -1;
    unsigned threads = 0;   // 0 => default to all cores

    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--info" || a == "-i") info = true;
        else if (a == "--gzip" || a == "-z") gzip = true;
        else if (a == "--no-gzip" || a == "--plain") gzip = false;
        else if (a == "--bin" || a == "-b") {
            if (i + 1 >= argc) die("--bin needs a bin number");
            only_bin = std::atol(argv[++i]);
        }
        else if (a == "--all_bins") {
            if (i + 1 >= argc) die("--all_bins needs an output directory");
            all_bins = true;
            all_bins_dir = argv[++i];
        }
        else if (a == "--threads" || a == "-t") {
            if (i + 1 >= argc) die("--threads needs a number");
            threads = (unsigned)std::strtoul(argv[++i], nullptr, 10);
        }
        else if (a == "-h" || a == "--help") { usage(std::cout, argv[0]); return 0; }
        else if (!a.empty() && a[0] == '-') die("unknown option: " + a);
        else if (path.empty()) path = a;
        else die("more than one input file given");
    }
    if (path.empty()) { usage(std::cerr, argv[0]); return 1; }
    if (all_bins && only_bin >= 0) die("--all_bins and --bin are mutually exclusive");

    std::ifstream in(path, std::ios::binary);
    if (!in) die("cannot open " + path);

    // A pack is read by seeking to its footer index, so it must be the raw,
    // uncompressed .kbin. If the pipeline compressed the pack (--compress_kbin_packs
    // leaves a .kbin.gz), tell the user to decompress it first rather than failing
    // later with a cryptic "bad magic". Detect gzip by its two-byte magic (1f 8b).
    unsigned char m[2] = {0, 0};
    in.read(reinterpret_cast<char*>(m), 2);
    if (in.gcount() == 2 && m[0] == 0x1f && m[1] == 0x8b) {
        std::string out = path;
        if (out.size() > 3 && out.compare(out.size() - 3, 3, ".gz") == 0)
            out.resize(out.size() - 3);
        else
            out += ".kbin";
        die(path + " is gzip-compressed. kbin_dump needs the uncompressed pack "
            "(it seeks by the footer index). Decompress it first, e.g.:\n"
            "    pigz -dc " + path + " > " + out);
    }
    in.clear();   // reset EOF/fail flags from the 2-byte probe before seeking

    in.seekg(0, std::ios::end);
    long long size = in.tellg();
    if (size < (long long)sizeof(Footer)) die(path + " is too small to be a pack file");

    Footer f{};
    in.seekg(size - (long long)sizeof(Footer), std::ios::beg);
    in.read(reinterpret_cast<char*>(&f), sizeof(f));
    if (!in) die("cannot read footer of " + path);
    if (std::memcmp(f.magic, "KGWS", 4) != 0)
        die(path + " is not a pack file (bad magic); expected a .kbin from kmer_count");

    const int    k         = f.kmer_k;
    const size_t rb        = f.record_bytes;
    const size_t key_bytes = rb >= 2 ? rb - 2 : 0;
    const int    BITS      = 2 * k;
    if (key_bytes != 8 && key_bytes != 16)
        die("unexpected record width " + std::to_string(rb) + " (expected 10 or 18 bytes)");

    std::vector<uint64_t> off(f.num_bins + 1);
    in.seekg((long long)f.index_offset, std::ios::beg);
    in.read(reinterpret_cast<char*>(off.data()),
            (std::streamsize)(off.size() * sizeof(uint64_t)));
    if (!in) die("cannot read the offset table of " + path);
    const size_t total = off.back() / rb;

    if (info) {
        std::cout << "file         : " << path << "\n"
                  << "format       : KGWS pack, version " << f.version << "\n"
                  << "k            : " << k << "\n"
                  << "record bytes : " << rb << " (" << key_bytes
                  << "-byte key + 2-byte count)\n"
                  << "bins         : " << f.num_bins << "\n"
                  << "total k-mers : " << total << "\n";
        if (f.num_bins) {
            size_t mn = SIZE_MAX, mx = 0;
            for (uint32_t b = 0; b < f.num_bins; b++) {
                size_t n = (off[b + 1] - off[b]) / rb;
                mn = std::min(mn, n);
                mx = std::max(mx, n);
            }
            std::cout << "k-mers/bin   : min " << mn << ", mean " << (total / f.num_bins)
                      << ", max " << mx << "\n";
        }
        return 0;
    }

    static const char BASE[4] = {'A', 'C', 'G', 'T'};
    const size_t FLUSH_AT   = 1u << 20;                        // flush text every ~1 MB
    const size_t CHUNK_RECS = std::max<size_t>(1, (1u << 20) / rb); // read ~1 MB of records at a time

    // Stream one bin from `fin` (seeking to its slice) and decode it to text,
    // flushing through `sink(data, len)`. Both the input read and the text output
    // are chunked at ~1 MB, so peak memory is a few MB regardless of bin size.
    // `kmer`/`txt`/`buf` are reusable scratch owned by the caller.
    auto stream_bin = [&](std::ifstream& fin, uint32_t b, std::string& kmer,
                          std::string& txt, std::vector<unsigned char>& buf, auto&& sink) {
        size_t    remaining = (off[b + 1] - off[b]) / rb;
        long long rpos      = (long long)off[b];
        while (remaining) {
            size_t take = std::min(remaining, CHUNK_RECS);
            buf.resize(take * rb);
            fin.seekg(rpos, std::ios::beg);
            fin.read(reinterpret_cast<char*>(buf.data()), (std::streamsize)(take * rb));
            if (!fin) throw std::runtime_error("short read on bin " + std::to_string(b));
            for (size_t r = 0; r < take; r++) {
                const unsigned char* rec = buf.data() + r * rb;
                uint64_t w0 = 0, w1 = 0;
                std::memcpy(&w0, rec, 8);
                if (key_bytes == 16) std::memcpy(&w1, rec + 8, 8);
                uint16_t count;
                std::memcpy(&count, rec + key_bytes, 2);
                for (int i = 0; i < k; i++) {
                    int p = BITS - 2 - 2 * i;             // first base in the highest bits
                    uint64_t bits = (p >= 64) ? ((w1 >> (p - 64)) & 3ULL) : ((w0 >> p) & 3ULL);
                    kmer[(size_t)i] = BASE[bits];
                }
                txt.append(kmer);
                txt += '\t';
                append_count(txt, count);
                txt += '\n';
                if (txt.size() >= FLUSH_AT) { sink(txt.data(), txt.size()); txt.clear(); }
            }
            rpos      += (long long)(take * rb);
            remaining -= take;
        }
        if (!txt.empty()) { sink(txt.data(), txt.size()); txt.clear(); }
    };

    // -------- --all_bins: one file per bin, decoded (and optionally gzipped) in parallel -----
    if (all_bins) {
        // Each accession's bins go in their own subfolder named after the pack
        // (the .kbin prefix), so extracting several accessions into one --all_bins
        // directory never collides.
        std::string acc = std::filesystem::path(path).stem().string();
        std::filesystem::path out_dir = std::filesystem::path(all_bins_dir) / acc;
        std::error_code ec;
        std::filesystem::create_directories(out_dir, ec);
        if (ec) die("cannot create directory " + out_dir.string() + ": " + ec.message());

        unsigned T = threads ? threads : std::max(1u, std::thread::hardware_concurrency());
        if (T > f.num_bins) T = std::max(1u, f.num_bins);

        std::atomic<uint32_t> next{0};
        std::exception_ptr err;
        std::mutex err_mx;

        auto worker = [&]() {
            std::ifstream fin(path, std::ios::binary);   // one handle per worker
            std::vector<unsigned char> buf;
            std::string kmer((size_t)k, 'A'), txt;
            txt.reserve(FLUSH_AT + 128);
            try {
                if (!fin) throw std::runtime_error("cannot open " + path);
                uint32_t b;
                while ((b = next.fetch_add(1)) < f.num_bins) {
                    std::string base = (out_dir / (std::to_string(b) + ".tsv")).string();
                    if (gzip) {
                        std::string outp = base + ".gz";
                        gzFile gf = gzopen(outp.c_str(), "wb");
                        if (!gf) throw std::runtime_error("cannot create " + outp);
                        gzbuffer(gf, 1u << 20);
                        auto sink = [&](const char* d, size_t len) {
                            if (len && gzwrite(gf, d, (unsigned)len) == 0)
                                throw std::runtime_error("gzwrite failed on " + outp);
                        };
                        stream_bin(fin, b, kmer, txt, buf, sink);
                        gzclose(gf);
                    } else {
                        std::ofstream ofs(base, std::ios::binary);
                        if (!ofs) throw std::runtime_error("cannot create " + base);
                        auto sink = [&](const char* d, size_t len) { ofs.write(d, (std::streamsize)len); };
                        stream_bin(fin, b, kmer, txt, buf, sink);
                        ofs.close();
                    }
                }
            } catch (...) {
                std::lock_guard<std::mutex> lk(err_mx);
                if (!err) err = std::current_exception();
            }
        };

        if (T <= 1) {
            worker();
        } else {
            std::vector<std::thread> pool;
            pool.reserve(T);
            for (unsigned t = 0; t < T; t++) pool.emplace_back(worker);
            for (auto& th : pool) th.join();
        }
        if (err) { try { std::rethrow_exception(err); } catch (const std::exception& e) { die(e.what()); } }
        return 0;
    }

    // -------- default / --bin: single stream to stdout -----------------------
    uint32_t b0 = 0, b1 = f.num_bins;
    if (only_bin >= 0) {
        if (only_bin >= (long)f.num_bins)
            die("bin " + std::to_string(only_bin) + " out of range (pack has "
                + std::to_string(f.num_bins) + ")");
        b0 = (uint32_t)only_bin;
        b1 = (uint32_t)only_bin + 1;
    }

    std::vector<unsigned char> buf;
    std::string kmer((size_t)k, 'A'), txt;
    txt.reserve(FLUSH_AT + 128);
    auto sink = [&](const char* d, size_t len) { std::fwrite(d, 1, len, stdout); };
    try {
        for (uint32_t b = b0; b < b1; b++) stream_bin(in, b, kmer, txt, buf, sink);
    } catch (const std::exception& e) { die(e.what()); }
    return 0;
}
