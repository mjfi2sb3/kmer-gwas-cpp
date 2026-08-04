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
// per-line overhead of Python matters; same interface and byte-identical output.
//
// Input is read transparently whether plain or gzip-compressed. Output is plain
// by default; --gzip (or a ".gz" output name) gzips it, but that gzip is
// single-threaded, so for large matrices "plain | pigz" is faster. Both default
// to stdin/stdout.
//
//   # decode (bits -> text): N is required to trim the byte padding
//   bits_to_text --decode -a accessions.txt -i bin0_bits.tsv.gz -o bin0_text.tsv
//   # encode (text -> bits): N and the delimiter are inferred from the row
//   bits_to_text --encode -i bin0_text.tsv -o bin0_bits.tsv

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
        "usage: %s [--decode|--encode] [-i in] [-o out] [options]\n"
        "\n"
        "  Convert a presence/absence matrix between bit-packed and text forms.\n"
        "  Input may be plain or gzip (auto-detected). Output is plain by default.\n"
        "\n"
        "  Common:\n"
        "    -i <file>       input  (default: stdin; plain or gzip, auto-detected)\n"
        "    -o <file>       output (default: stdout)\n"
        "    --gzip, -z      gzip the output (also implied by a .gz output name).\n"
        "                    Single-threaded; for big files, prefer plain | pigz.\n"
        "    --delimiter     text-side value separator: tab | space | none\n"
        "                    (default: none; on --encode it is auto-detected per row)\n"
        "\n"
        "  --decode  bits -> text (default): <kmer><TAB><hex> -> <kmer><TAB><vals>\n"
        "    -a <file>       accessions file; its non-blank line count gives N\n"
        "    -n <N>          accession count, if the accessions file is unavailable\n"
        "                    (N is REQUIRED for --decode: it trims the byte padding)\n"
        "\n"
        "  --encode  text -> bits: <kmer><TAB><vals> -> <kmer><TAB><hex>\n"
        "                    N is inferred from the row width; -a/-n are optional and,\n"
        "                    if given, validate every row against that count.\n",
        prog);
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
    string accessions_path, infile, outfile, delim = "";   // default delimiter: none
    long n_acc = -1;
    bool encode = false;        // default: decode (bits -> text)
    bool gzip   = false;        // default: plain (gzip here is single-threaded; see --gzip)
    bool delim_set = false;     // did the user pass --delimiter explicitly?

    vector<string> pos;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        auto need = [&](const char* what) -> string {
            if (i + 1 >= argc) die(string("missing value for ") + what);
            return argv[++i];
        };
        if      (a == "--decode") encode = false;
        else if (a == "--encode") encode = true;
        else if (a == "-i" || a == "--input")  infile  = need("-i");
        else if (a == "-o" || a == "--output") outfile = need("-o");
        else if (a == "--gzip" || a == "-z") gzip = true;
        else if (a == "-a" || a == "--accessions") accessions_path = need("-a");
        else if (a == "-n" || a == "--num-accessions") n_acc = atol(need("-n").c_str());
        else if (a == "--delimiter") {
            string d = need("--delimiter");
            if      (d == "tab")   delim = "\t";
            else if (d == "space") delim = " ";
            else if (d == "none")  delim = "";
            else die("--delimiter must be tab, space or none");
            delim_set = true;
        }
        else if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
        else if (!a.empty() && a[0] == '-' && a != "-") die("unknown option: " + a);
        else pos.push_back(a);
    }
    // Positional args fill in any -i/-o not given (back-compat: in [out]).
    size_t p = 0;
    if (infile.empty()  && p < pos.size()) infile  = pos[p++];
    if (outfile.empty() && p < pos.size()) outfile = pos[p++];
    if (p < pos.size()) die("too many positional arguments");
    if (infile == "-")  infile.clear();
    if (outfile == "-") outfile.clear();

    // Resolve the accession count. Required for --decode; optional for --encode.
    if (!accessions_path.empty()) {
        ifstream f(accessions_path);
        if (!f) die("cannot open accessions file: " + accessions_path);
        string line; long cnt = 0;
        while (getline(f, line)) {
            if (line.find_first_not_of(" \t\r\n") != string::npos) cnt++;
        }
        n_acc = cnt;
    }
    if (!encode && n_acc <= 0) {
        usage(argv[0]);
        die("--decode requires the accession count (use -a or -n) to trim byte padding");
    }

    gzFile in = infile.empty() ? gzdopen(dup(STDIN_FILENO), "rb")
                               : gzopen(infile.c_str(), "rb");
    if (!in) die("cannot open input" + (infile.empty() ? "" : ": " + infile));

    // Output: plain by default; --gzip forces it, and a ".gz" output name implies it.
    if (!outfile.empty() && outfile.size() >= 3 && outfile.compare(outfile.size() - 3, 3, ".gz") == 0)
        gzip = true;
    gzFile gout = nullptr; FILE* fout = nullptr;
    if (gzip) {
        if (outfile.empty()) gout = gzdopen(dup(STDOUT_FILENO), "wb");
        else {
            if (outfile.size() < 3 || outfile.compare(outfile.size() - 3, 3, ".gz") != 0)
                outfile += ".gz";
            gout = gzopen(outfile.c_str(), "wb");
        }
        if (!gout) die("cannot open output" + (outfile.empty() ? "" : ": " + outfile));
    } else {
        fout = outfile.empty() ? stdout : fopen(outfile.c_str(), "wb");
        if (!fout) die("cannot open output: " + outfile);
    }
    auto put = [&](const string& s) {
        if (gzip) { if (!s.empty() && gzwrite(gout, s.data(), (unsigned)s.size()) == 0) die("write failed"); }
        else      { if (fwrite(s.data(), 1, s.size(), fout) != s.size()) die("write failed"); }
    };

    // hex-nibble lookup for decoding
    unsigned char hexval[256];
    memset(hexval, 0xFF, sizeof(hexval));
    for (int c = '0'; c <= '9'; c++) hexval[c] = c - '0';
    for (int c = 'a'; c <= 'f'; c++) hexval[c] = 10 + c - 'a';
    for (int c = 'A'; c <= 'F'; c++) hexval[c] = 10 + c - 'A';
    static const char HEX[] = "0123456789abcdef";

    size_t N       = n_acc > 0 ? (size_t)n_acc : 0;   // 0 = infer (encode only)
    size_t n_bytes = N ? (N + 7) / 8 : 0;
    vector<unsigned char> bytes(n_bytes);

    string obuf; obuf.reserve(64 * 1024);
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
            // On the first row, auto-detect the delimiter (unless given) from the
            // values, and infer N from the row width (unless -a/-n fixed it).
            if (!delim_set && lineno == 1) {
                if      (memchr(vals, '\t', vlen)) delim = "\t";
                else if (memchr(vals, ' ',  vlen)) delim = " ";
                else                               delim = "";
            }
            // A value is "present" unless it is exactly "0"; this collapses a
            // count matrix to presence/absence too.
            size_t upper = vlen + 1;                 // >= number of values in this row
            if (bytes.size() < (upper + 7) / 8) bytes.assign((upper + 7) / 8, 0);
            else memset(bytes.data(), 0, bytes.size());
            size_t count = 0;
            auto set_bit = [&](bool present) {
                if (N && count >= N)
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
                        bool present = !(i - start == 1 && vals[start] == '0');
                        set_bit(present);
                        start = i + 1;
                    }
                }
            }
            if (N == 0) { N = count; n_bytes = (N + 7) / 8; }   // infer from first row
            else if (count != N)
                die("line " + to_string(lineno) + ": row has " + to_string(count) +
                    " values but " + to_string(N) + " accessions were expected. "
                    "Wrong -a/-n or --delimiter.");
            obuf.assign(line.data(), tab);
            obuf += '\t';
            for (size_t b = 0; b < n_bytes; b++) { obuf += HEX[bytes[b] >> 4]; obuf += HEX[bytes[b] & 15]; }
        }
        obuf += '\n';
        put(obuf);
    }

    gzclose(in);
    if (gzip) gzclose(gout);
    else if (fout != stdout) fclose(fout);
    return 0;
}
