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
#include <string>
#include <cstring>
#include <chrono>
#include "mmap_io.hpp"
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <filesystem>
#include <getopt.h>
#include <thread>
#include <atomic>
#include <memory>

using namespace std;

// Sparse matrix entry: stores only non-zero values
struct SparseEntry {
    ushort acc_index;
    ushort count;

    SparseEntry(ushort idx, ushort cnt) : acc_index(idx), count(cnt) {}
};

// Thread-safe sparse matrix
class SparseMatrix {
private:
    unordered_map<string, vector<SparseEntry>> data_;
    vector<unique_ptr<mutex>> mutexes_;  // Fine-grained locking by hash
    static const size_t NUM_LOCKS = 256;  // Reduced to avoid resource exhaustion

    size_t get_lock_index(const string& key) const {
        return hash<string>{}(key) % NUM_LOCKS;
    }

public:
    SparseMatrix() {
        mutexes_.reserve(NUM_LOCKS);
        for (size_t i = 0; i < NUM_LOCKS; i++) {
            mutexes_.emplace_back(make_unique<mutex>());
        }
    }

    void insert(const string& kmer, ushort acc_index, ushort count) {
        if (kmer.empty()) return;  // Safety check

        size_t lock_idx = get_lock_index(kmer);
        lock_guard<mutex> lock(*mutexes_[lock_idx]);

        // More efficient: try to find first, only iterate if key exists
        auto it = data_.find(kmer);
        if (it == data_.end()) {
            // New kmer - just insert
            data_[kmer] = vector<SparseEntry>{SparseEntry(acc_index, count)};
        } else {
            // Existing kmer - check if this accession already has an entry
            auto& entries = it->second;
            for (auto& entry : entries) {
                if (entry.acc_index == acc_index) {
                    entry.count = count;
                    return;
                }
            }
            entries.emplace_back(acc_index, count);
        }
    }

    const unordered_map<string, vector<SparseEntry>>& get_data() const {
        return data_;
    }

    size_t size() const {
        return data_.size();
    }
};

vector<string> split(const string &str, const char &sep)
{
    string segment;
    vector<string> ret;
    istringstream ss(str);
    while (getline(ss, segment, sep))
        ret.push_back(segment);
    return ret;
}

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

// Worker function for parallel file reading
void process_accession_file(
    const string& file_path,
    size_t acc_index,
    SparseMatrix& matrix,
    atomic<size_t>& files_processed,
    atomic<size_t>& total_kmers_read
) {
    try {
        if (!filesystem::exists(file_path)) {
            throw runtime_error("File: " + file_path + " does not exist!");
        }

        ifstream stream(file_path);
        if (!stream) {
            throw runtime_error("Failed to open file: " + file_path);
        }

        size_t local_kmer_count = 0;
        string line;
        while (getline(stream, line)) {
            if (line.empty()) continue;  // Skip empty lines

            auto record = split(line, '\t');
            if (record.size() != 2) {
                throw runtime_error("Invalid format in file: " + file_path + " (line: " + line + ")");
            }

            const string& key = record[0];
            if (key.empty()) continue;  // Skip empty keys

            ushort value = stoi(record[1]);
            if (value == 0) continue;  // Skip zero counts (optimization)

            matrix.insert(key, acc_index, value);
            local_kmer_count++;
        }

        stream.close();
        total_kmers_read += local_kmer_count;
        files_processed++;

    } catch (const exception& e) {
        cerr << "Error processing " << file_path << ": " << e.what() << endl;
        throw;
    }
}

void merge_chunk(
    const uint file_index,
    const uint min_occur,
    string input_path,
    string accessions_path,
    string delimiter,
    bool show_count,
    string output_dir,
    size_t num_threads
) {
    auto accessions = get_accessions(accessions_path);
    size_t NUM_ACC = accessions.size();

    cout << "Using " << num_threads << " threads for parallel processing" << endl;

    // Validate all files exist before processing
    cout << "Validating file existence..." << endl;
    for (size_t i = 0; i < accessions.size(); i++) {
        string file_path = input_path + accessions[i] + "/" + to_string(file_index) + "_nr.tsv";
        if (!filesystem::exists(file_path)) {
            throw runtime_error("File: " + file_path + " does not exist!");
        }
    }
    cout << "All files validated successfully" << endl;

    // Sparse matrix for efficient memory usage
    SparseMatrix matrix;

    // Thread pool for parallel file reading
    vector<thread> threads;
    atomic<size_t> files_processed(0);
    atomic<size_t> total_kmers_read(0);
    atomic<size_t> next_acc_index(0);
    atomic<bool> error_occurred(false);
    string error_message;
    mutex error_mutex;

    auto start_reading = chrono::steady_clock::now();

    // Launch worker threads
    auto worker = [&]() {
        try {
            while (true) {
                if (error_occurred.load()) break;  // Stop on error

                size_t acc_index = next_acc_index.fetch_add(1);
                if (acc_index >= accessions.size()) break;

                string file_path = input_path + accessions[acc_index] + "/" +
                                 to_string(file_index) + "_nr.tsv";

                process_accession_file(file_path, acc_index, matrix, files_processed, total_kmers_read);

                // Progress update (thread-safe)
                size_t processed = files_processed.load();
                if (processed % 100 == 0 || processed == accessions.size()) {
                    lock_guard<mutex> lock(error_mutex);
                    cout << "Processed " << processed << "/" << accessions.size()
                         << " accessions..." << endl;
                }
            }
        } catch (const exception& e) {
            error_occurred.store(true);
            lock_guard<mutex> lock(error_mutex);
            error_message = e.what();
        }
    };

    threads.reserve(num_threads);
    for (size_t i = 0; i < num_threads; i++) {
        try {
            threads.emplace_back(worker);
        } catch (const system_error& e) {
            cerr << "Failed to create thread " << i << ": " << e.what() << endl;
            cerr << "Continuing with " << threads.size() << " threads..." << endl;
            break;
        }
    }

    if (threads.empty()) {
        throw runtime_error("Failed to create any worker threads");
    }

    cout << "Successfully launched " << threads.size() << " worker threads" << endl;

    // Wait for all threads
    for (auto& t : threads) {
        t.join();
    }

    // Check if any errors occurred during processing
    if (error_occurred.load()) {
        throw runtime_error("Error during parallel processing: " + error_message);
    }

    auto end_reading = chrono::steady_clock::now();
    cout << "File reading took "
         << chrono::duration_cast<chrono::seconds>(end_reading - start_reading).count()
         << " seconds" << endl;
    cout << "Total kmers read: " << total_kmers_read.load() << endl;
    cout << "Unique kmers: " << matrix.size() << endl;

    // Write output files
    auto start_writing = chrono::steady_clock::now();

    ofstream m_stream(output_dir + to_string(file_index) + "_matrix.tsv");
    ofstream ck_stream(output_dir + to_string(file_index) + "_core.txt");

    if (!m_stream || !ck_stream) {
        throw runtime_error("Failed to open output files");
    }

    size_t core_kmers = 0;
    size_t output_kmers = 0;
    size_t filtered_kmers = 0;

    const auto& data = matrix.get_data();

    for (const auto& pair_ : data) {
        const string& kmer = pair_.first;
        const vector<SparseEntry>& entries = pair_.second;
        size_t num_accessions_with_kmer = entries.size();

        // Check if it's a core kmer (present in all accessions)
        if (num_accessions_with_kmer == NUM_ACC) {
            ck_stream << kmer << "\n";
            core_kmers++;
        }

        // Apply threshold filter: skip kmers that are too rare or too common
        if (min_occur > 0) {
            if (num_accessions_with_kmer < min_occur || num_accessions_with_kmer > NUM_ACC - min_occur) {
                filtered_kmers++;
                continue;  // Skip this kmer
            }
        }

        // Build dense representation for output
        m_stream << kmer << "\t";

        // Create dense vector from sparse entries
        vector<ushort> dense(NUM_ACC, 0);
        for (const auto& entry : entries) {
            dense[entry.acc_index] = entry.count;
        }

        // Write output
        for (size_t i = 0; i < NUM_ACC; i++) {
            ushort freq = dense[i];
            if (!show_count && freq > 0) freq = 1;
            m_stream << delimiter << freq;
        }
        m_stream << "\n";
        output_kmers++;
    }

    m_stream.close();
    ck_stream.close();

    auto end_writing = chrono::steady_clock::now();
    cout << "File writing took "
         << chrono::duration_cast<chrono::seconds>(end_writing - start_writing).count()
         << " seconds" << endl;
    cout << "Core kmers (present in all accessions): " << core_kmers << endl;
    if (min_occur > 0) {
        cout << "Filtered kmers (threshold=" << min_occur << "): " << filtered_kmers << endl;
        cout << "Kmers written after filtering: " << output_kmers << endl;
    } else {
        cout << "Total kmers written: " << output_kmers << endl;
    }
}

int main(int argc, char *argv[])
{
    string input_path, accessions_path;
    uint file_index;
    uint min_occur = 0;
    std::string delimiter = "\t";  // Default value: tab
    bool show_count = true;  // Default is to show counts
    size_t num_threads = thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 4;  // Fallback
    // Conservative defaults for stability
    if (num_threads > 32) num_threads = 32;  // Lower default cap for memory safety

    // Define the long options
    static struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"accessions", required_argument, 0, 'a'},
        {"index", required_argument, 0, 'f'},
        {"threshold", required_argument, 0, 't'},
        {"delimiter", required_argument, 0, 'd'},
        {"count", required_argument, 0, 'c'},
        {"threads", required_argument, 0, 'j'},
        {0, 0, 0, 0}
    };
    int option_index = 0;
    int c;

    while ((c = getopt_long(argc, argv, "i:a:f:t:d:c:j:", long_options, &option_index)) != -1)
    {
        switch (c)
        {
            case 'i':
                if (optarg == nullptr) {
                    cerr << "Error: --input requires an argument." << endl;
                    return -1;
                }
                input_path = optarg;
                break;
            case 'a':
                if (optarg == nullptr) {
                    cerr << "Error: --accessions requires an argument." << endl;
                    return -1;
                }
                accessions_path = optarg;
                break;
            case 'f':
                if (optarg == nullptr) {
                    cerr << "Error: --index requires an argument." << endl;
                    return -1;
                }
                if (!all_of(optarg, optarg + strlen(optarg), ::isdigit)) {
                    cerr << "Error: --index requires a valid integer argument." << endl;
                    return -1;
                }
                file_index = stoi(optarg);
                break;
            case 't':
                if (optarg == nullptr) {
                    cerr << "Error: --threshold requires an argument." << endl;
                    return -1;
                }
                min_occur = stoi(optarg);
                break;
            case '?':
                // getopt_long will print an error message
                break;
            case 'd':
                if (optarg == nullptr) {
                    cerr << "Error: --delimiter requires an argument." << endl;
                    return -1;
                }
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
                    cerr << "Invalid delimiter option. Use 'tab' or 'none'." << endl;
                    return -1;
                }
                break;
            case 'c':
                if (optarg == nullptr) {
                    cerr << "Error: --count requires an argument." << endl;
                    return -1;
                }
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
            case 'j':
                if (optarg == nullptr) {
                    cerr << "Error: --threads requires an argument." << endl;
                    return -1;
                }
                if (!all_of(optarg, optarg + strlen(optarg), ::isdigit)) {
                    cerr << "Error: --threads requires a valid integer argument." << endl;
                    return -1;
                }
                num_threads = stoi(optarg);
                if (num_threads == 0) {
                    cerr << "Error: --threads must be greater than 0." << endl;
                    return -1;
                }
                if (num_threads > 48) {
                    cout << "Warning: Capping threads at 48 (requested: " << num_threads << ")" << endl;
                    cout << "Note: For very large datasets, start with 8-16 threads to avoid memory issues" << endl;
                    num_threads = 48;
                }
                break;
            default:
                abort();
        }
    }
    // Check if all required options are provided
    if (input_path.empty() || accessions_path.empty())
    {
        cout << "\nOptimized Matrix Merge (Sparse + Parallel)\n"
             << "usage: " << argv[0] << "\n"
             << "\t\t--input <input path> \n"
             << "\t\t--accessions <accessions path> \n"
             << "\t\t--index <file index which corresponds to bin> \n"
             << "\t\t--threshold <min/max occurrence threshold> (default: " << min_occur << ")\n"
             << "\t\t            Filters out kmers appearing in < threshold or > (N-threshold) accessions\n"
             << "\t\t            Set to 0 to disable filtering (keep all kmers)\n"
             << "\t\t--delimiter <delimiter type: tab|none> (default: " << (delimiter == "\t" ? "tab" : (delimiter == " " ? "space" : "none")) << ")\n"
             << "\t\t--count <print matrix as absence/presence or actual k-mer counts; type: y|n> (default: " << (show_count ? "y" : "n") << ")\n"
             << "\t\t--threads <number of parallel threads> (default: " << num_threads << ", max: 48)\n"
             << "\n"
             << "Performance tips:\n"
             << "  - For large datasets (1000+ accessions), start with --threads 8 or 16\n"
             << "  - Monitor memory usage and increase threads if memory allows\n"
             << "  - Optimal thread count is usually 2-4x the number of physical cores\n\n";
        return -1;
    }

    try
    {
        if (input_path[input_path.size()-1] != '/')
            input_path += '/';
        auto accessions = get_accessions(accessions_path);
        size_t NUM_ACC = accessions.size();

        std::string output_dir = "";
        std::string delim = (delimiter == "\t" ? "tab" : (delimiter == " " ? "space" : "none"));
        if (show_count)
        {
            output_dir = "matrix_acc"+to_string(NUM_ACC)+"_count_delim-"+delim+"/";
        }
        else
        {
            output_dir = "matrix_acc"+to_string(NUM_ACC)+"_pres-abs_delim-"+delim+"/";
        }

        if (!filesystem::exists(output_dir)){
            filesystem::create_directory(output_dir);
        }

        cout << "***************************** " << endl;
        cout << "OPTIMIZED MATRIX MERGE (Sparse + Parallel)" << endl;
        cout << "PROCESSING MATRIX CHUNK: " << file_index << endl;
        cout << "Number of accessions: " << NUM_ACC << std::endl;
        cout << "Output Folder: " << output_dir << std::endl;
        cout << "Number of threads: " << num_threads << std::endl;
        if (min_occur > 0) {
            cout << "Threshold filtering: ENABLED (min_occur=" << min_occur
                 << ", keeping kmers in range [" << min_occur << ", " << (NUM_ACC - min_occur) << "])" << endl;
        } else {
            cout << "Threshold filtering: DISABLED (keeping all kmers)" << endl;
        }

        auto start = chrono::steady_clock::now();

        merge_chunk(file_index, min_occur, input_path, accessions_path, delimiter, show_count, output_dir, num_threads);

        auto end = chrono::steady_clock::now();
        cout << "Total processing time for index " << file_index << ": "
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
