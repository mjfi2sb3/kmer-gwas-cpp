// bits_to_text — convert a bits-encoded k-mer matrix (--encoding bits) to text.
//
// The bits encoding stores each row as
//     <k-mer><TAB><hex>
// where <hex> is ceil(N/8) bytes, one bit per accession, LSB-first within each
// byte (accession i is bit (i & 7) of byte (i >> 3)). This expands it back to
//     <k-mer><TAB><v0><delim><v1>...<v(N-1)>
//
// A C++ companion to tools/bits_to_text.py, for large matrices where the
// per-line overhead of the Python version matters. Same interface and
// byte-identical output. Input/output may be plain or gzip-compressed
// (transparently, via zlib); both default to stdin/stdout.
//
//   bits_to_text -a accessions.txt bin0_matrix.tsv.gz bin0_text.tsv.gz
//   pigz -dc bin0_matrix.tsv.gz | bits_to_text -n 12600 --delimiter tab | pigz > out.tsv.gz
//
// The number of accessions is required (via the accessions file, or -n),
// because the last hex byte may contain padding bits that must not be emitted.

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
        "usage: %s (-a <accessions> | -n <N>) [--delimiter tab|space|none]\n"
        "          [--header] [infile] [outfile]\n"
        "\n"
        "  Expand a bits-encoded matrix (<kmer><TAB><hex>) to text.\n"
        "  -a <file>       accessions file; non-blank lines give N and column order\n"
        "  -n <N>          accession count, if the accessions file is unavailable\n"
        "  --delimiter     value separator: tab (default) | space | none\n"
        "  --header        prefix an accession-name header row (requires -a)\n"
        "  infile/outfile  default stdin/stdout; .gz is handled transparently\n", prog);
}

// Read one line (without the trailing newline) from a gzFile of any length.
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

    vector<string> pos;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        auto need = [&](const char* what) -> string {
            if (i + 1 >= argc) die(string("missing value for ") + what);
            return argv[++i];
        };
        if      (a == "-a" || a == "--accessions") accessions_path = need("-a");
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

    // Resolve the accession count (and names, if a file was given).
    vector<string> names;
    if (!accessions_path.empty()) {
        ifstream f(accessions_path);
        if (!f) die("cannot open accessions file: " + accessions_path);
        string line;
        while (getline(f, line)) {
            // trim
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

    // Input: gzopen reads plain and gzip transparently. stdin via gzdopen.
    gzFile in = infile.empty() ? gzdopen(dup(STDIN_FILENO), "rb")
                               : gzopen(infile.c_str(), "rb");
    if (!in) die("cannot open input" + (infile.empty() ? "" : ": " + infile));

    // Output: gzip if the name ends in .gz, else plain.
    bool out_gz = outfile.size() > 3 && outfile.compare(outfile.size() - 3, 3, ".gz") == 0;
    gzFile   gout = nullptr;
    FILE*    fout = nullptr;
    if (out_gz) { gout = gzopen(outfile.c_str(), "wb"); if (!gout) die("cannot open output: " + outfile); }
    else        { fout = outfile.empty() ? stdout : fopen(outfile.c_str(), "wb"); if (!fout) die("cannot open output: " + outfile); }

    string obuf;                       // reused per line
    obuf.reserve(64 * 1024);
    auto flush_line = [&](const string& s) {
        if (out_gz) { if (gzwrite(gout, s.data(), (unsigned)s.size()) == 0) die("write failed"); }
        else        { if (fwrite(s.data(), 1, s.size(), fout) != s.size()) die("write failed"); }
    };

    if (header) {
        obuf = "kmer\t";
        for (size_t i = 0; i < names.size(); i++) { if (i) obuf += delim; obuf += names[i]; }
        obuf += '\n';
        flush_line(obuf);
    }

    // hex-nibble lookup: '0'-'9','a'-'f','A'-'F' -> 0..15, else 0xFF
    unsigned char hexval[256];
    memset(hexval, 0xFF, sizeof(hexval));
    for (int c = '0'; c <= '9'; c++) hexval[c] = c - '0';
    for (int c = 'a'; c <= 'f'; c++) hexval[c] = 10 + c - 'a';
    for (int c = 'A'; c <= 'F'; c++) hexval[c] = 10 + c - 'A';

    vector<unsigned char> bytes(n_bytes);
    string line;
    size_t lineno = 0;
    while (gz_getline(in, line)) {
        lineno++;
        if (line.empty()) continue;
        size_t tab = line.find('\t');
        if (tab == string::npos)
            die("line " + to_string(lineno) + ": expected '<kmer><TAB><hex>'");
        const char* hex = line.data() + tab + 1;
        size_t hexlen = line.size() - tab - 1;
        if (hexlen != n_bytes * 2)
            die("line " + to_string(lineno) + ": row has " + to_string(hexlen / 2) +
                " bytes but " + to_string(N) + " accessions need " + to_string(n_bytes) +
                ". Wrong -a/-n, or the file is not bits-encoded.");
        for (size_t b = 0; b < n_bytes; b++) {
            unsigned char hi = hexval[(unsigned char)hex[2 * b]];
            unsigned char lo = hexval[(unsigned char)hex[2 * b + 1]];
            if (hi == 0xFF || lo == 0xFF)
                die("line " + to_string(lineno) + ": invalid hex digit");
            bytes[b] = (unsigned char)((hi << 4) | lo);
        }
        obuf.assign(line.data(), tab);        // the k-mer
        obuf += '\t';
        if (delim.empty()) {
            for (size_t i = 0; i < N; i++)
                obuf += ((bytes[i >> 3] >> (i & 7)) & 1) ? '1' : '0';
        } else {
            char d = delim[0];
            for (size_t i = 0; i < N; i++) {
                if (i) obuf += d;
                obuf += ((bytes[i >> 3] >> (i & 7)) & 1) ? '1' : '0';
            }
        }
        obuf += '\n';
        flush_line(obuf);
    }

    gzclose(in);
    if (out_gz) gzclose(gout);
    else if (fout != stdout) fclose(fout);
    return 0;
}
