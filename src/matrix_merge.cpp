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
#include "mmap_io.hpp"
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <filesystem>
#include <getopt.h>
#include <unistd.h>


using namespace std;

const int k = 51;
const map<string, char> bases = {
    {"00", 'A'}, {"01", 'C'}, {"10", 'G'}, {"11", 'T'}
};
string bit_decode(bitset<2*k> code) {
    string kmer = "";
    string bcode = code.to_string();
    for (size_t i = 0; i < bcode.size(); i += 2)
        kmer += bases.at(bcode.substr(i, 2));
    return kmer;
}

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

void merge_chunk(const uint file_index, const uint min_occur, string input_path,
                 string accessions_path, string delimiter, bool show_count,
                 string ouput_dir, bool write_core, uint num_threads = 1)
{
    const uint   MAX_THREADS = 128;
    const size_t NUM_SHARDS  = 1024;  // must be power of 2

    auto accessions = get_accessions(accessions_path);
    size_t NUM_ACC = accessions.size();

    // Check all bin files exist before reading any
    for (auto &acc : accessions)
    {
        string file_path = input_path + acc + "/" + to_string(file_index) + "_nr.bin";
        if (!filesystem::exists(file_path))
            throw runtime_error("File: " + file_path + " does not exist!");
    }

    // Map a kmer to a shard index using its lower 10 bits (0–1023)
    auto shard_of = [](const bitset<2*k>& kmer) -> size_t {
        size_t id = 0;
        for (size_t i = 0; i < 10; i++)
            id |= ((size_t)kmer[i] << i);
        return id;
    };

    // Unified output: write one kmer row to the matrix and core files
    ofstream m_stream(ouput_dir + to_string(file_index) + "_matrix.tsv");
    ofstream ck_stream;
    if (write_core) ck_stream.open(ouput_dir + to_string(file_index) + "_core.txt");

    auto write_row = [&](const bitset<2*k>& kmer, const vector<ushort>& row) {
        if (write_core && row[NUM_ACC] == NUM_ACC)
            ck_stream << bit_decode(kmer) << "\n";
        if (row[NUM_ACC] < min_occur || row[NUM_ACC] > NUM_ACC - min_occur) return;
        m_stream << bit_decode(kmer) << "\t";
        for (int i = 0; i < (int)NUM_ACC; i++) {
            auto freq = row[i];
            if (!show_count && freq > 0) freq = 1;
            m_stream << delimiter << freq;
        }
        m_stream << "\n";
    };

    if (num_threads > 1)
    {
        // --- Sharded map parallel path ---
        // Each thread reads its accessions and inserts directly into the
        // appropriate shard. Threads only contend when hitting the same shard
        // simultaneously (probability ~1/NUM_SHARDS per insert).
        uint actual_threads = min({num_threads, (uint)NUM_ACC, MAX_THREADS});

        using ShardMap = unordered_map<bitset<2*k>, vector<ushort>>;
        vector<ShardMap>  shards(NUM_SHARDS);
        vector<mutex>     shard_mutexes(NUM_SHARDS);

        vector<thread>       threads;
        vector<exception_ptr> thread_errors(actual_threads, nullptr);
        size_t chunk = (NUM_ACC + actual_threads - 1) / actual_threads;

        for (uint t = 0; t < actual_threads; t++)
        {
            size_t t_start = t * chunk;
            size_t t_end   = min(t_start + chunk, NUM_ACC);
            if (t_start >= NUM_ACC) break;

            threads.emplace_back([&, t, t_start, t_end]()
            {
                try {
                    for (size_t i = t_start; i < t_end; i++)
                    {
                        string file_path = input_path + accessions[i] + "/"
                                           + to_string(file_index) + "_nr.bin";
                        ifstream stream(file_path, ios::binary);
                        if (!stream)
                            throw runtime_error("Failed to open file: " + file_path);

                        bitset<2*k> key;
                        ushort value;
                        while (stream.read(reinterpret_cast<char*>(&key), sizeof(key))) {
                            stream.read(reinterpret_cast<char*>(&value), sizeof(value));
                            size_t s = shard_of(key);
                            lock_guard<mutex> lock(shard_mutexes[s]);
                            if (shards[s].find(key) == shards[s].end())
                                shards[s][key] = vector<ushort>(NUM_ACC+1, 0);
                            shards[s][key][i] = value;
                            shards[s][key][NUM_ACC]++;
                        }
                    }
                } catch (...) {
                    thread_errors[t] = current_exception();
                }
            });
        }

        for (auto &th : threads) th.join();
        for (auto &ep : thread_errors)
            if (ep) rethrow_exception(ep);

        // Output: iterate all shards sequentially
        for (size_t s = 0; s < NUM_SHARDS; s++)
            for (auto &pair_ : shards[s])
                write_row(pair_.first, pair_.second);
    }
    else
    {
        // --- Single-threaded path (original behaviour) ---
        unordered_map<bitset<2*k>, vector<ushort>> matrix_;
        size_t acc_index = 0;

        for (auto &acc : accessions)
        {
            string file_path = input_path + acc + "/" + to_string(file_index) + "_nr.bin";
            ifstream stream(file_path, ios::binary);
            if (!stream)
                throw runtime_error("Failed to open file: " + file_path);

            bitset<2*k> key;
            ushort value;
            while (stream.read(reinterpret_cast<char*>(&key), sizeof(key))) {
                stream.read(reinterpret_cast<char*>(&value), sizeof(value));
                if (matrix_.find(key) == matrix_.end())
                    matrix_[key] = vector<ushort>(NUM_ACC+1, 0);
                matrix_[key][acc_index] = value;
                matrix_[key][NUM_ACC]++;
            }
            stream.close();
            acc_index++;
        }

        for (auto &pair_ : matrix_)
            write_row(pair_.first, pair_.second);
    }

    m_stream.close();
    if (write_core) ck_stream.close();
}

int main(int argc, char *argv[])
{
	string input_path, accessions_path;
	uint file_index;
	uint min_occur = 0;
	std::string delimiter = "\t";  // Default value: tab
	bool show_count = false;  // Default is to show presence/absence
	bool write_core = false;  // Default is to skip core k-mers file
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
            	if (!all_of(optarg, optarg + strlen(optarg), ::isdigit)) {
					cerr << "Error: --index requires a valid integer argument." << endl;
					return -1;
				}
				file_index = stoi(optarg);
				break;
                /*file_index = stoi(optarg);
                break;*/
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
				else 
				{
				    cerr << "Invalid delimiter option. Use 'tab' or 'space'." << endl;
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

    // Check if all required options are provided
    if (input_path.empty() || accessions_path.empty()) 
    /*{
        cout << "usage: " << argv[0] << " --input <input path> --accessions <accessions path> --index <file index> --threshold <min occurence threshold>\n";
        return -1;
    }*/
    {
		cout << "\nusage: " << argv[0] << "\n"
		     << "\t\t--input <input path> \n"
		     << "\t\t--accessions <accessions path> \n"
		     << "\t\t--index <file index which corresponds to bin> \n"
		     << "\t\t--threshold <value> (min/max occurence threshold, default: " << min_occur << ")\n"
		     << "\t\t            Note: Use --threshold 20 or --threshold=20 (both work)\n"
		     << "\t\t--delimiter <delimiter type: tab|none> (default: " << (delimiter == "\t" ? "tab" : (delimiter == " " ? "space" : "none")) << ")\n"
		     << "\t\t--count <print matrix as absence/presence or actual k-mer counts; type: y|n> (default: " << (show_count ? "y" : "n") << ")\n"
		     << "\t\t--core    <write core k-mers file (_core.txt); type: y|n> (default: n)\n"
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
   		std::string delim = (delimiter == "\t" ? "tab" : (delimiter == " " ? "space" : "none"));
		if (show_count)
		{
			ouput_dir = "matrix_acc"+to_string(NUM_ACC)+"_bins"+to_string(num_bins)+"_minOcc"+to_string(min_occur)+"_count_delim-"+delim+"/";
		}
		else
		{
			ouput_dir = "matrix_acc"+to_string(NUM_ACC)+"_bins"+to_string(num_bins)+"_minOcc"+to_string(min_occur)+"_pres-abs_delim-"+delim+"/";
		}
   		
		if (!filesystem::exists(ouput_dir)){
			filesystem::create_directory(ouput_dir);
		}

        cout << "***************************** " << endl;
        cout << "PROCESSING MATRIX CHUNK: " << file_index << endl;
        cout << "Number of accessions: " << NUM_ACC << std::endl;
        cout << "Output Folder: " << ouput_dir << std::endl;
        
        auto start = chrono::steady_clock::now();

        merge_chunk(file_index, min_occur, input_path, accessions_path, delimiter, show_count, ouput_dir, write_core, num_threads);

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

