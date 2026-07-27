// bits_to_text — convert a k-mer matrix between the bit-packed and text forms.
//
// The two forms of a presence/absence matrix row are
//     text:  <k-mer><TAB><v0><delim><v1>...<v(N-1)>     (each v is 0 or 1)
//     bits:  <k-mer><TAB><hex>                          (N bits, hex-encoded)
// where in the bits form <hex> is ceil(N/8) bytes, one bit per accession,
// LSB-first within each byte (accession i is bit (i & 7) of byte (i >> 3)).
//
// Runs both ways:
//     --decode  (default)  bits -> text   (expand)
//     --encode             text -> bits   (pack; a matrix "compressor")
//
// A C++ companion to tools/bits_to_text.py, for large matrices where the
// per-line overhead of Python matters; same interface and byte-identical
// output. Input/output may be plain or gzip-compressed transparently (zlib);
// both default to stdin/stdout, and an output name ending in .gz is written
// gzip-compressed.
//
//   # decode (bits -> text)
//   bits_to_text -a accessions.txt bin0_bits.tsv.gz bin0_text.tsv.gz
//   # encode (text -> bits)
//   bits_to_text --encode -a accessions.txt bin0_text.tsv.gz bin0_bits.tsv.gz
//
// The accession count is required (via the accessions file, or -n): when
// decoding it trims padding bits; when encoding it validates the row width.
// --delimiter is the separator of the TEXT side in both directions (how the
// values are written when decoding, and how they are read when encoding).

#include <zlib.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <unistd.h>

using namespace std;

static void die(const string& msg) { fprintf(stderr, "bits_to_text: %s\n", msg.c_str()); exit(1); }

static void usage(const char* prog) {
    fprintf(stderr,
        "usage: %s [--decode|--encode] (-a <accessions> | -n <N>)\n"
        "          [--delimiter tab|space|none] [--header] [infile] [outfile]\n"
        "\n"
        "  Convert a presence/absence matrix between bit-packed and text forms.\n"
        "  --decode        bits -> text (default): <kmer><TAB><hex> -> <kmer><TAB><vals>\n"
        "  --encode        text -> bits: <kmer><TAB><vals> -> <kmer><TAB><hex>\n"
        "  -a <file>       accessions file; non-blank lines give N and column order\n"
        "  -n <N>          accession count, if the accessions file is unavailable\n"
        "  --delimiter     separator of the TEXT side: tab (default) | space | none\n"
        "  --header        (decode only) prefix an accession-name header row; needs -a\n"
        "  infile/outfile  default stdin/stdout; .gz is handled transparently\n", prog);
}

// Read one line (without the trailing newline) from a gzFile, any length.
static bool gz_getline(gzFile f, string& line) {
    line.clear();
    char buf[65536];
    while (gzgets(f, buf, sizeof(buf))) {
        line += buf;
        if (!line.empty() && line.back() == '\n') { line.pop_back(); return true; }
    }
    return !line.empty();
}

int main(int argc, char** argv) {
    string accessions_path, infile, outfile, delim = "\t";
    long n_acc = -1;
    bool header = false;
    bool encode = false;             // default: decode (bits -> text)

    vector<string> pos;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        auto need = [&](const char* what) -> string {
            if (i + 1 >= argc) die(string("missing value for ") + what);
            return argv[++i];
        };
        if      (a == "--decode") encode = false;
        else if (a == "--encode") encode = true;
        else if (a == "-a" || a == "--accessions") accessions_path = need("-a");
        else if (a == "-n" || a == "--num-accessions") n_acc = atol(need("-n").c_str());
        else if (a == "--delimiter") {
            string d = need("--delimiter");
            if      (d == "tab")   delim = "\t";
            else if (d == "space") delim = " ";
            else if (d == "none")  delim = "";
            else die("--delimiter must be tab, space or none");
        }
        else if (a == "--header") header = true;
        else if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
        else if (!a.empty() && a[0] == '-' && a != "-") die("unknown option: " + a);
        else pos.push_back(a);
    }
    if (pos.size() > 0) infile  = pos[0];
    if (pos.size() > 1) outfile = pos[1];
    if (pos.size() > 2) die("too many positional arguments");
    if (header && encode) die("--header applies to --decode only");
    // "-" is the usual convention for stdin/stdout; treat it as such.
    if (infile == "-")  infile.clear();
    if (outfile == "-") outfile.clear();

    // Resolve the accession count (and names, if a file was given).
    vector<string> names;
    if (!accessions_path.empty()) {
        ifstream f(accessions_path);
        if (!f) die("cannot open accessions file: " + accessions_path);
        string line;
        while (getline(f, line)) {
            size_t b = line.find_first_not_of(" \t\r\n");
            if (b == string::npos) continue;
            size_t e = line.find_last_not_of(" \t\r\n");
            names.push_back(line.substr(b, e - b + 1));
        }
        n_acc = (long)names.size();
    }
    if (n_acc <= 0) { usage(argv[0]); die("accession count is required and must be positive (use -a or -n)"); }
    if (header && names.empty()) die("--header requires -a (to know the accession names)");
    const size_t N = (size_t)n_acc;
    const size_t n_bytes = (N + 7) / 8;

    gzFile in = infile.empty() ? gzdopen(dup(STDIN_FILENO), "rb")
                               : gzopen(infile.c_str(), "rb");
    if (!in) die("cannot open input" + (infile.empty() ? "" : ": " + infile));

    bool out_gz = outfile.size() > 3 && outfile.compare(outfile.size() - 3, 3, ".gz") == 0;
    gzFile gout = nullptr; FILE* fout = nullptr;
    if (out_gz) { gout = gzopen(outfile.c_str(), "wb"); if (!gout) die("cannot open output: " + outfile); }
    else        { fout = outfile.empty() ? stdout : fopen(outfile.c_str(), "wb"); if (!fout) die("cannot open output: " + outfile); }

    auto put = [&](const string& s) {
        if (out_gz) { if (!s.empty() && gzwrite(gout, s.data(), (unsigned)s.size()) == 0) die("write failed"); }
        else        { if (fwrite(s.data(), 1, s.size(), fout) != s.size()) die("write failed"); }
    };

    string obuf; obuf.reserve(64 * 1024);

    if (header && !encode) {
        obuf = "kmer\t";
        for (size_t i = 0; i < names.size(); i++) { if (i) obuf += delim; obuf += names[i]; }
        obuf += '\n';
        put(obuf);
    }

    // hex-nibble lookup for decoding
    unsigned char hexval[256];
    memset(hexval, 0xFF, sizeof(hexval));
    for (int c = '0'; c <= '9'; c++) hexval[c] = c - '0';
    for (int c = 'a'; c <= 'f'; c++) hexval[c] = 10 + c - 'a';
    for (int c = 'A'; c <= 'F'; c++) hexval[c] = 10 + c - 'A';
    static const char HEX[] = "0123456789abcdef";

    vector<unsigned char> bytes(n_bytes);
    string line;
    size_t lineno = 0;

    while (gz_getline(in, line)) {
        lineno++;
        if (line.empty()) continue;
        size_t tab = line.find('\t');
        if (tab == string::npos)
            die("line " + to_string(lineno) + ": expected '<kmer><TAB>...'");
        const char* vals = line.data() + tab + 1;
        size_t vlen = line.size() - tab - 1;

        if (!encode) {
            // ---- decode: hex -> 0/1 values -------------------------------
            if (vlen != n_bytes * 2)
                die("line " + to_string(lineno) + ": row has " + to_string(vlen / 2) +
                    " bytes but " + to_string(N) + " accessions need " + to_string(n_bytes) +
                    ". Wrong -a/-n, or the input is not bits-encoded.");
            for (size_t b = 0; b < n_bytes; b++) {
                unsigned char hi = hexval[(unsigned char)vals[2 * b]];
                unsigned char lo = hexval[(unsigned char)vals[2 * b + 1]];
                if (hi == 0xFF || lo == 0xFF)
                    die("line " + to_string(lineno) + ": invalid hex digit");
                bytes[b] = (unsigned char)((hi << 4) | lo);
            }
            obuf.assign(line.data(), tab);
            obuf += '\t';
            if (delim.empty())
                for (size_t i = 0; i < N; i++)
                    obuf += ((bytes[i >> 3] >> (i & 7)) & 1) ? '1' : '0';
            else {
                char d = delim[0];
                for (size_t i = 0; i < N; i++) {
                    if (i) obuf += d;
                    obuf += ((bytes[i >> 3] >> (i & 7)) & 1) ? '1' : '0';
                }
            }
        } else {
            // ---- encode: 0/1 values -> hex -------------------------------
            // A value is "present" unless it is exactly "0"; this collapses a
            // count matrix to presence/absence too. Values are read according
            // to --delimiter (single chars for 'none', delimited fields else).
            for (size_t b = 0; b < n_bytes; b++) bytes[b] = 0;
            size_t count = 0;
            auto set_bit = [&](bool present) {
                if (count >= N)
                    die("line " + to_string(lineno) + ": more than " + to_string(N) + " values");
                if (present) bytes[count >> 3] |= (unsigned char)(1u << (count & 7));
                count++;
            };
            if (delim.empty()) {
                for (size_t i = 0; i < vlen; i++) set_bit(vals[i] != '0');
            } else {
                char d = delim[0];
                size_t start = 0;
                for (size_t i = 0; i <= vlen; i++) {
                    if (i == vlen || vals[i] == d) {
                        // one field is [start, i)
                        bool present = !(i - start == 1 && vals[start] == '0');
                        set_bit(present);
                        start = i + 1;
                    }
                }
            }
            if (count != N)
                die("line " + to_string(lineno) + ": row has " + to_string(count) +
                    " values but " + to_string(N) + " accessions were expected. Wrong "
                    "-a/-n or --delimiter.");
            obuf.assign(line.data(), tab);
            obuf += '\t';
            for (size_t b = 0; b < n_bytes; b++) { obuf += HEX[bytes[b] >> 4]; obuf += HEX[bytes[b] & 15]; }
        }
        obuf += '\n';
        put(obuf);
    }

    gzclose(in);
    if (out_gz) gzclose(gout);
    else if (fout != stdout) fclose(fout);
    return 0;
}
