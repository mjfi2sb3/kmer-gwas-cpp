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
//    kbin_dump <file.kbin>              all k-mers, every bin, to stdout
//    kbin_dump <file.kbin> --info       header summary only (k, bins, totals)
//    kbin_dump <file.kbin> --bin N       only bin N
//    kbin_dump <file.kbin> --with-bin    prepend a bin-index column
//    kbin_dump <file.kbin> > out.tsv     redirect; pipe to gzip/pigz to compress
//
//  The pack layout (see src/pack_io.hpp) is: bin records back to back, then a
//  (num_bins + 1) x uint64 offset table, then a 24-byte footer ending in the
//  magic "KGWS". Each record is a 2-bit-per-base packed key (8 bytes for k <= 32,
//  16 bytes otherwise) followed by a uint16 count. The first base of the k-mer
//  sits in the highest bits, base values A=0 C=1 G=2 T=3, matching decode() in
//  src/kmer_key.hpp. Words are stored little-endian.
// ===========================================================================

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

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
    o << "usage: " << prog << " <file.kbin> [--info] [--bin N] [--with-bin]\n"
      << "\n"
      << "  Reads a .kbin pack (kmer_count's Stage 1 output) using its self-describing\n"
      << "  footer, so one binary reads a pack of any k. With no option it writes every\n"
      << "  k-mer as '<kmer><TAB><count>' to stdout.\n"
      << "\n"
      << "    --info       print a header summary (k, bins, record width, totals) and exit\n"
      << "    --bin N      dump only bin N\n"
      << "    --with-bin   prepend a '<bin><TAB>' column to each row\n"
      << "\n"
      << "  Output is plain text; pipe to gzip or pigz to compress a large dump.\n";
}

int main(int argc, char** argv) {
    std::string path;
    bool info = false, with_bin = false;
    long only_bin = -1;

    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--info" || a == "-i") info = true;
        else if (a == "--with-bin") with_bin = true;
        else if (a == "--bin" || a == "-b") {
            if (i + 1 >= argc) die("--bin needs a bin number");
            only_bin = std::atol(argv[++i]);
        }
        else if (a == "-h" || a == "--help") { usage(std::cout, argv[0]); return 0; }
        else if (!a.empty() && a[0] == '-') die("unknown option: " + a);
        else if (path.empty()) path = a;
        else die("more than one input file given");
    }
    if (path.empty()) { usage(std::cerr, argv[0]); return 1; }

    std::ifstream in(path, std::ios::binary);
    if (!in) die("cannot open " + path);

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

    uint32_t b0 = 0, b1 = f.num_bins;
    if (only_bin >= 0) {
        if (only_bin >= (long)f.num_bins)
            die("bin " + std::to_string(only_bin) + " out of range (pack has "
                + std::to_string(f.num_bins) + ")");
        b0 = (uint32_t)only_bin;
        b1 = (uint32_t)only_bin + 1;
    }

    static const char BASE[4] = {'A', 'C', 'G', 'T'};
    std::string kmer((size_t)k, 'A');
    std::vector<unsigned char> buf;
    std::ios::sync_with_stdio(false);

    for (uint32_t b = b0; b < b1; b++) {
        size_t bytes = off[b + 1] - off[b];
        size_t n     = bytes / rb;
        if (!n) continue;
        buf.resize(bytes);
        in.seekg((long long)off[b], std::ios::beg);
        in.read(reinterpret_cast<char*>(buf.data()), (std::streamsize)bytes);
        if (!in) die("short read on bin " + std::to_string(b) + " of " + path);

        for (size_t r = 0; r < n; r++) {
            const unsigned char* rec = buf.data() + r * rb;
            uint64_t w0 = 0, w1 = 0;
            std::memcpy(&w0, rec, 8);
            if (key_bytes == 16) std::memcpy(&w1, rec + 8, 8);
            uint16_t count;
            std::memcpy(&count, rec + key_bytes, 2);

            for (int i = 0; i < k; i++) {
                int p = BITS - 2 - 2 * i;               // first base in the highest bits
                uint64_t bits = (p >= 64) ? ((w1 >> (p - 64)) & 3ULL) : ((w0 >> p) & 3ULL);
                kmer[(size_t)i] = BASE[bits];
            }
            if (with_bin) std::cout << b << '\t';
            std::cout << kmer << '\t' << count << '\n';
        }
    }
    return 0;
}
