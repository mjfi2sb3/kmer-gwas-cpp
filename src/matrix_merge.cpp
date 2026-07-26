//#include "thread_pool.hpp"
#include <iostream>
#include <istream>
#include <streambuf>
#include <fstream>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <vector>
#include <bitset>
#include <mutex>
#include <thread>
#include <string>
#include <cstring>
#include <chrono>
#include <stdexcept>
#include "kmer_key.hpp"
#include "pack_io.hpp"
#include <queue>
#include <memory>
#include <sys/resource.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <filesystem>
#include <getopt.h>
#include <unistd.h>


using namespace std;

const int k = kmer::K;

using Key = kmer::Key;

inline string bit_decode(const Key &code) { return kmer::decode(code); }

// const size_t NUM_ACC = 100;
vector<string> split(const string &str, const char &sep)
{
    string segment;
    vector<string> ret;
    istringstream ss(str);
    while (getline(ss, segment, sep))
        ret.push_back(segment);
    return ret;
}

/*vector<string> get_accessions(string accessions_path)
{
        vector<string> accessions;
	cout << "accessions path: " << accessions_path << endl;
	try {
        	ifstream stream(accessions_path);
        	string line;
        	while (getline(stream, line))
        	    accessions.push_back(line);
        	stream.close();
	}
    catch (exception const &e)
        	{cout << e.what() << endl;}
        	
	if (accessions.size() == 0) 
		cout << "could not read accessions, make sure file is not empty." << endl;
	return accessions;
}*/

vector<string> get_accessions(const string& accessions_path) {
    vector<string> accessions;
    cout << "Accessions path: " << accessions_path << endl;

    try {
        ifstream stream(accessions_path);
        if (!stream) {
            throw runtime_error("Failed to open accessions file: " + accessions_path);
        }

        string line;
        while (getline(stream, line)) {
            accessions.push_back(line);
        }

        stream.close();
    } catch (const exception& e) {
        cerr << e.what() << endl;
        throw;
    }

    if (accessions.empty()) {
        throw runtime_error("Could not read accessions. Make sure the file is not empty.");
    }

    return accessions;
}

// ===========================================================================
// merge_chunk — build one bin's matrix by k-way merging the accessions'
// sorted pack slices.
//
// MEMORY: O(number of accessions), NOT O(union k-mers x accessions).
//
// The previous implementation loaded the whole bin into
// `unordered_map<Key, vector<Count>>`, so it needed union_rows x (N+1) x 2
// bytes — 25.2 KB per distinct k-mer at 12,600 accessions, i.e. 8 GB per bin
// task at a 0.5 B union and 82 GB at 5 B. That is why num_bins had to grow with
// genome size, and num_bins is what drove the inode count.
//
// Because Stage 1 now emits each bin slice SORTED, the same result can be
// produced by advancing one cursor per accession through a heap: pop the
// smallest key, gather the counts from every accession holding it, emit the
// row. Resident memory is one small buffer per accession (N x 4096 records)
// plus a single reusable row.
//
// Note the row itself is inherently O(N) to WRITE — a dense matrix row is
// 25.2 KB at 12,600 accessions regardless of how it was computed. That is the
// output format's cost, addressed separately in task #15.
// ===========================================================================
void merge_chunk(const uint file_index, const uint min_occur, string input_path,
                 string accessions_path, string delimiter, bool show_count,
                 string ouput_dir, bool write_core, bool bit_packed, uint num_threads = 1)
{
    (void)num_threads;   // the merge is sequential by nature; see note below

    auto accessions = get_accessions(accessions_path);
    size_t NUM_ACC = accessions.size();

    // The accession-occurrence counter is a 16-bit kmer::Count.
    if (NUM_ACC > kmer::COUNT_MAX)
        throw runtime_error("Too many accessions (" + to_string(NUM_ACC) + "): the "
                            "occurrence counter is 16-bit, maximum is "
                            + to_string(kmer::COUNT_MAX) + ". Widen kmer::Count to raise this.");

    // One descriptor per accession stays open for the whole merge. Raise the
    // soft limit to the hard limit so large panels work without the caller
    // having to remember `ulimit -n`.
    {
        struct rlimit rl;
        if (getrlimit(RLIMIT_NOFILE, &rl) == 0 && rl.rlim_cur < rl.rlim_max) {
            rl.rlim_cur = rl.rlim_max;
            setrlimit(RLIMIT_NOFILE, &rl);
        }
        if (getrlimit(RLIMIT_NOFILE, &rl) == 0 && rl.rlim_cur < NUM_ACC + 16)
            throw runtime_error("Need at least " + to_string(NUM_ACC + 16)
                                + " file descriptors for " + to_string(NUM_ACC)
                                + " accessions, but the limit is " + to_string(rl.rlim_cur)
                                + ". Raise it with `ulimit -n`.");
    }

    // Check every pack exists before opening any, so a missing accession fails
    // immediately rather than part-way through a long merge.
    vector<string> paths(NUM_ACC);
    for (size_t i = 0; i < NUM_ACC; i++) {
        paths[i] = input_path + accessions[i] + ".kbin";
        if (!filesystem::exists(paths[i]))
            throw runtime_error("Pack file does not exist: " + paths[i]);
    }

    ofstream m_stream(ouput_dir + to_string(file_index) + "_matrix.tsv");
    ofstream ck_stream;
    if (write_core) ck_stream.open(ouput_dir + to_string(file_index) + "_core.txt");

    // Reusable row. Only the entries actually touched by a k-mer are written
    // and then cleared, so per-row work is O(accessions carrying the k-mer),
    // not O(N) — zeroing 25.2 KB per row would dominate everything else.
    vector<kmer::Count> row(NUM_ACC, 0);
    vector<unsigned char> packed;
    vector<uint32_t>    touched;
    touched.reserve(256);

    auto write_row = [&](const Key& key, size_t occ) {
        // Core k-mers (present in every accession) are invariant across the
        // panel and carry no association signal: recorded separately and
        // dropped from the matrix. The file was originally named _discard.txt.
        if (write_core && occ == NUM_ACC)
            ck_stream << bit_decode(key) << "\n";
        // Two-sided MAF filter: discard both rare and near-ubiquitous k-mers.
        // min_occur == 0 disables it. Intentional — see commit cafcc67.
        if (occ < min_occur || occ > NUM_ACC - min_occur) return;

        // Row layout: <kmer> TAB <v0> D <v1> D ... D <vN-1>
        // The k-mer is always tab-separated from the values so the first column
        // can be split off even when --delimiter none packs the values
        // together; `delimiter` separates values from EACH OTHER only.
        m_stream << bit_decode(key) << "\t";
        if (bit_packed) {
            // One bit per accession, LSB-first within each byte, written as
            // lowercase hex. 12,600 accessions become 1,575 bytes -> 3,150 hex
            // characters, against 25,200 characters for the tab form: 8x smaller.
            // Presence/absence only, so there is nothing to lose by packing.
            packed.assign((NUM_ACC + 7) / 8, 0);
            for (size_t i = 0; i < NUM_ACC; i++)
                if (row[i]) packed[i >> 3] |= (unsigned char)(1u << (i & 7));
            static const char HEX[] = "0123456789abcdef";
            for (unsigned char b : packed) { m_stream << HEX[b >> 4] << HEX[b & 15]; }
        } else {
            for (size_t i = 0; i < NUM_ACC; i++) {
                kmer::Count freq = row[i];
                if (!show_count && freq > 0) freq = 1;
                if (i) m_stream << delimiter;
                m_stream << freq;
            }
        }
        m_stream << "\n";
    };

    // -- Open one cursor per accession -------------------------------------
    vector<unique_ptr<kmer::BinCursor>> cur;
    cur.reserve(NUM_ACC);
    for (size_t i = 0; i < NUM_ACC; i++)
        cur.push_back(make_unique<kmer::BinCursor>(paths[i], file_index));

    // -- Min-heap over the cursors' current keys ---------------------------
    // Greater-than comparison so priority_queue yields the SMALLEST key.
    struct HeapItem { Key key; uint32_t idx; };
    auto worse = [](const HeapItem& a, const HeapItem& b) {
        if (a.key == b.key) return a.idx > b.idx;
        return b.key < a.key;
    };
    priority_queue<HeapItem, vector<HeapItem>, decltype(worse)> heap(worse);

    for (size_t i = 0; i < NUM_ACC; i++)
        if (cur[i]->valid()) heap.push({cur[i]->current().key, (uint32_t)i});

    size_t rows = 0;
    while (!heap.empty())
    {
        Key key = heap.top().key;
        touched.clear();

        // Drain every accession positioned on this key.
        while (!heap.empty() && heap.top().key == key) {
            uint32_t i = heap.top().idx;
            heap.pop();
            row[i] = cur[i]->current().count;
            touched.push_back(i);
            cur[i]->advance();
            if (cur[i]->valid()) heap.push({cur[i]->current().key, i});
        }

        write_row(key, touched.size());
        rows++;
        for (uint32_t i : touched) row[i] = 0;    // clear only what was set
    }

    cout << "  merged " << rows << " distinct k-mers across "
         << NUM_ACC << " accessions" << endl;

    m_stream.close();
    if (write_core) ck_stream.close();
}

int main(int argc, char *argv[])
{
	string input_path, accessions_path;
	uint file_index = 0;
	bool have_index = false;      // --index is required; see the guard below
	uint min_occur = 0;
	std::string delimiter = "\t";  // Default value: tab
	bool show_count = false;  // Default is to show presence/absence
	bool write_core = false;  // Default is to skip core k-mers file
	bool bit_packed = false;  // --delimiter bits: 1 bit per accession
	uint num_bins = 0;
	uint available_threads;
	if (const char* s = getenv("SLURM_CPUS_PER_TASK")) {
	    uint cpus = (uint)atoi(s);
	    available_threads = cpus > 1 ? cpus : max(1u, thread::hardware_concurrency());
	} else {
	    available_threads = max(1u, thread::hardware_concurrency());
	}
	uint num_threads = available_threads;

	// Define the long options
	static struct option long_options[] = {
		{"input",    required_argument, 0, 'i'},
		{"accessions", required_argument, 0, 'a'},
		{"index",    required_argument, 0, 'f'},
		{"threshold", required_argument, 0, 't'},
		{"delimiter", required_argument, 0, 'd'},
		{"count",    required_argument, 0, 'c'},
		{"core",     required_argument, 0, 'r'},
		{"bins",     required_argument, 0, 'b'},
		{"threads",  required_argument, 0, 'T'},
		{0, 0, 0, 0}
	};
    int option_index = 0;
    int c;

    while ((c = getopt_long(argc, argv, "i:a:f:t:d:c:r:b:T:", long_options, &option_index)) != -1)
    {
        switch (c) 
        {
            case 'i':
                input_path = optarg;
                break;
            case 'a':
                accessions_path = optarg;
                break;
            case 'f':
            	if (*optarg == '\0' || !all_of(optarg, optarg + strlen(optarg), ::isdigit)) {
					cerr << "Error: --index requires a valid non-negative integer argument." << endl;
					return -1;
				}
				file_index = (uint)stoul(optarg);
				have_index = true;
				break;
            case 't':
                min_occur = stoi(optarg);
                break;
            case '?':
                // getopt_long will print an error message
                break;
            case 'd':
				if (string(optarg) == "tab")
				{
				    delimiter = "\t";
				}
				else if (string(optarg) == "none")
				{
				    delimiter = "";
				}
				else if (string(optarg) == "bits")
				{
				    // Bit-packed presence/absence: 1 bit per accession instead of
				    // 1-2 bytes. See the note in write_row for the layout.
				    delimiter = "";
				    bit_packed = true;
				}
				else
				{
				    cerr << "Invalid delimiter option. Use 'tab', 'none' or 'bits'." << endl;
				    return -1;
				}
				break;
            case 'c':
				if (string(optarg) == "y") 
				{
				    show_count = true;
				} 
				else if (string(optarg) == "n") 
				{
				    show_count = false;
				} 
				else 
				{
				    cerr << "Invalid display option. Use 'y' or 'n'." << endl;
				    return -1;
				}
				break;
            case 'r':
				if (string(optarg) == "y")
				    write_core = true;
				else if (string(optarg) == "n")
				    write_core = false;
				else {
				    cerr << "Invalid core option. Use 'y' or 'n'." << endl;
				    return -1;
				}
				break;
            case 'b':
                num_bins = stoi(optarg);
                break;
            case 'T':
                num_threads = stoi(optarg);
                break;
            default:
                abort();
        }
    }
    // Cap threads at available resources
    if (num_threads > available_threads) {
        cerr << "Warning: --threads " << num_threads << " exceeds available CPUs ("
             << available_threads << "), capping at " << available_threads << "." << endl;
        num_threads = available_threads;
    }

    // `--delimiter none` packs values with no separator, which is only
    // decodable when every value is a single character. With --count y the
    // counts are multi-digit and the row becomes unparseable (task #5).
    if (bit_packed && show_count) {
        cerr << "Error: --delimiter bits is presence/absence only and cannot be "
                "combined with --count y. Use --delimiter tab for counts." << endl;
        return -1;
    }
    if (delimiter.empty() && !bit_packed && show_count) {
        cerr << "Error: --delimiter none cannot be combined with --count y — "
                "multi-digit counts would be concatenated with no separator and "
                "the matrix could not be parsed. Use --delimiter tab for counts, "
                "or --count n for presence/absence." << endl;
        return -1;
    }

    // Check if all required options are provided.
    // --index must be included: it was previously read from an uninitialised
    // variable when omitted, silently merging an arbitrary bin (task #10).
    if (input_path.empty() || accessions_path.empty() || !have_index)
    /*{
        cout << "usage: " << argv[0] << " --input <input path> --accessions <accessions path> --index <file index> --threshold <min occurence threshold>\n";
        return -1;
    }*/
    {
		cout << "\nusage: " << argv[0] << "\n"
		     << "\t\t--input <input path> \n"
		     << "\t\t--accessions <accessions path> \n"
		     << "\t\t--index <file index which corresponds to bin> \n"
		     << "\t\t--threshold <value> (two-sided MAF filter, default: " << min_occur << ")\n"
		     << "\t\t            Keeps k-mers carried by [value, n_accessions - value] accessions.\n"
		     << "\t\t            Filters on accession occurrence, NOT on per-accession count. 0 = off.\n"
		     << "\t\t            Note: Use --threshold 20 or --threshold=20 (both work)\n"
		     << "\t\t--delimiter <tab|none|bits> (default: " << (delimiter == "\t" ? "tab" : (delimiter == " " ? "space" : "none")) << ")\n"
		     << "\t\t--count <print matrix as absence/presence or actual k-mer counts; type: y|n> (default: " << (show_count ? "y" : "n") << ")\n"
		     << "\t\t--core    <write core k-mers file (_core.txt); type: y|n> (default: n)\n"
		     << "\t\t            Core = present in ALL accessions; these are excluded from the matrix.\n"
		     << "\t\t            'bits' packs presence/absence 1 bit per accession as hex:\n"
		     << "\t\t            ~8x smaller than 'tab' at large cohorts. Excludes --count y.\n"
		     << "\t\t--bins    <number of bins (used in output folder name)> (default: 0)\n"
		     << "\t\t--threads <parallel threads for accession reading> (default: SLURM_CPUS_PER_TASK if set, else hardware concurrency)\n\n";
		return -1;
	}

    /*if (argc != 5)
    {
        cout << "usage: " << argv[0] << " <input path> <accessions path> <file index> <min occurence threshold (across panel)>\n";
        return -1;
    }
	*/
    try
    {
      	/*string input_path = argv[1];
		string accessions_path = argv[2];
        uint file_index = stoi(argv[3]);
		uint min_occur = stoi(argv[4]);*/
		
      	if (input_path[input_path.size()-1] != '/')
		input_path += '/';
    	auto accessions = get_accessions(accessions_path);
   		size_t NUM_ACC = accessions.size();

   		std::string ouput_dir = "";
   		std::string delim = bit_packed ? "bits" : (delimiter == "\t" ? "tab" : (delimiter == " " ? "space" : "none"));
		if (show_count)
		{
			ouput_dir = "matrix_acc"+to_string(NUM_ACC)+"_k"+to_string(kmer::K)+"_bins"+to_string(num_bins)+"_minOcc"+to_string(min_occur)+"_count_delim-"+delim+"/";
		}
		else
		{
			ouput_dir = "matrix_acc"+to_string(NUM_ACC)+"_k"+to_string(kmer::K)+"_bins"+to_string(num_bins)+"_minOcc"+to_string(min_occur)+"_pres-abs_delim-"+delim+"/";
		}
   		
		if (!filesystem::exists(ouput_dir)){
			filesystem::create_directory(ouput_dir);
		}

        cout << "***************************** " << endl;
        cout << "PROCESSING MATRIX CHUNK: " << file_index << endl;
        cout << "Number of accessions: " << NUM_ACC << std::endl;
        cout << "Output Folder: " << ouput_dir << std::endl;
        
        auto start = chrono::steady_clock::now();

        merge_chunk(file_index, min_occur, input_path, accessions_path, delimiter, show_count, ouput_dir, write_core, bit_packed, num_threads);

        auto end = chrono::steady_clock::now();
        cout << "processing index: " << file_index << " took "
             << chrono::duration_cast<chrono::seconds>(end - start).count()
             << " sec" << endl;

        cout << "FINISHED MATRIX CHUNK: " << file_index << endl;
        cout << "***************************** " << endl;
    }

    catch (exception const &e)
    {
        cerr << e.what() << endl;
        return EXIT_FAILURE; // This ensures a non-zero exit code
    }
    return EXIT_SUCCESS; // This ensures a zero exit code for success
}

