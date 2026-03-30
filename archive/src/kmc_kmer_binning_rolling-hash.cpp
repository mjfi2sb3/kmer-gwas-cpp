#include "thread_pool.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <mutex>
#include <string>
#include <cstring>
#include <chrono>
#include <filesystem>
#include <getopt.h>

// KMC API includes
#include "/ibex/sw/rl9c/kmc/3.2.1/rl9_conda3/KMC/kmc_api/kmer_api.h"
#include "/ibex/sw/rl9c/kmc/3.2.1/rl9_conda3/KMC/kmc_api/kmc_file.h"

using namespace std;

// Function to split a string by a delimiter
vector<string> split(const string &str, const char &sep) {
    string segment;
    vector<string> ret;
    istringstream ss(str);
    while (getline(ss, segment, sep))
        ret.push_back(segment);
    return ret;
}

// Function to process a KMC database and distribute k-mers across bins
void process_kmc_database(const string& kmc_db_path, const string& output_path, 
                         const string& accession, size_t num_bins) {
    cout << "Processing KMC database: " << kmc_db_path << endl;
    
    // Create output directories
    filesystem::create_directory(output_path);
    filesystem::create_directory(output_path + "/" + accession);
    
    // Prepare output file streams for each bin
    vector<ofstream> bin_streams(num_bins);
    for (size_t i = 0; i < num_bins; i++) {
        string bin_path = output_path + "/" + accession + "/" + to_string(i) + "_nr.tsv";
        bin_streams[i].open(bin_path);
        if (!bin_streams[i].is_open()) {
            cerr << "Error: Could not open output file: " << bin_path << endl;
            exit(1);
        }
    }
    
    // Open the KMC database
    CKMCFile kmc_db;
    if (!kmc_db.OpenForListing(kmc_db_path)) {
        cerr << "Error: Cannot open KMC database: " << kmc_db_path << endl;
        exit(1);
    }
    
    // Get database info
    CKMCFileInfo db_info;
    kmc_db.Info(db_info);
    
    cout << "KMC database info:" << endl;
    cout << "  k-mer length: " << db_info.kmer_length << endl;
    cout << "  total k-mers: " << db_info.total_kmers << endl;
    cout << "  min count: " << db_info.min_count << endl;
    cout << "  max count: " << db_info.max_count << endl;
    
    // Create a k-mer object for iteration
    CKmerAPI kmer(db_info.kmer_length);
    uint32 count;
    
    // Variables for progress reporting
    uint64 processed = 0;
    uint64 total = db_info.total_kmers;
    uint64 report_step = total / 20; // Report progress 20 times
    if (report_step == 0) report_step = 1;
    
    auto start_time = chrono::steady_clock::now();
    
    // Iterate through all k-mers and distribute them to bins
    while (kmc_db.ReadNextKmer(kmer, count)) {
        // Get k-mer as string
        string kmer_str = kmer.to_string();
        
        // Generate a hash to determine which bin this k-mer belongs to
        // Using a simple hash function based on the first few characters of the k-mer
        size_t hash_val = 0;
        for (size_t i = 0; i < min(size_t(8), kmer_str.length()); i++) {
            hash_val = hash_val * 131 + kmer_str[i];
        }
        size_t bin_idx = hash_val % num_bins;
        
        // Write k-mer and count to the appropriate bin file
        bin_streams[bin_idx] << kmer_str << "\t" << count << "\n";
        
        // Report progress
        processed++;
        if (processed % report_step == 0) {
            auto current_time = chrono::steady_clock::now();
            auto elapsed = chrono::duration_cast<chrono::seconds>(current_time - start_time).count();
            cout << "Progress: " << (processed * 100 / total) << "% (" 
                 << processed << "/" << total << " k-mers) - " 
                 << elapsed << " seconds elapsed" << endl;
        }
    }
    
    // Close all bin files
    for (auto& stream : bin_streams) {
        stream.close();
    }
    
    // Close KMC database
    kmc_db.Close();
    
    auto end_time = chrono::steady_clock::now();
    auto elapsed = chrono::duration_cast<chrono::seconds>(end_time - start_time).count();
    cout << "Completed processing of " << processed << " k-mers in " << elapsed << " seconds" << endl;
}

// Process multiple KMC databases in parallel using a thread pool
void process_multiple_databases(const vector<string>& kmc_db_paths, 
                               const vector<string>& accessions,
                               const string& output_path, 
                               size_t num_bins,
                               size_t num_threads) {
    thread_pool pool(num_threads);
    
    if (kmc_db_paths.size() != accessions.size()) {
        cerr << "Error: Number of KMC databases and accessions must match" << endl;
        exit(1);
    }
    
    cout << "Processing " << kmc_db_paths.size() << " KMC databases using " 
         << num_threads << " threads" << endl;
    
    auto start_time = chrono::steady_clock::now();
    
    // Submit each database for processing
    for (size_t i = 0; i < kmc_db_paths.size(); i++) {
        const string& kmc_db_path = kmc_db_paths[i];
        const string& accession = accessions[i];
        
        pool.push_task([kmc_db_path, output_path, accession, num_bins]() {
            process_kmc_database(kmc_db_path, output_path, accession, num_bins);
        });
    }
    
    // Wait for all tasks to complete
    pool.wait_for_tasks();
    
    auto end_time = chrono::steady_clock::now();
    auto elapsed = chrono::duration_cast<chrono::seconds>(end_time - start_time).count();
    cout << "All databases processed in " << elapsed << " seconds" << endl;
}

// Function to read accessions from a file
vector<string> read_accessions_file(const string& file_path) {
    vector<string> accessions;
    ifstream file(file_path);
    
    if (!file.is_open()) {
        cerr << "Error: Could not open accessions file: " << file_path << endl;
        exit(1);
    }
    
    string line;
    while (getline(file, line)) {
        // Skip empty lines
        if (!line.empty()) {
            accessions.push_back(line);
        }
    }
    
    if (accessions.empty()) {
        cerr << "Warning: No accessions found in " << file_path << endl;
    } else {
        cout << "Read " << accessions.size() << " accessions from " << file_path << endl;
    }
    
    return accessions;
}

int main(int argc, char *argv[]) {
    string kmc_input_dir = "";
    string output_dir = "";
    string accessions_file = "";
    size_t num_bins = 50;
    size_t num_threads = 0; // 0 means use all available
    bool use_accession_file = false;
    vector<string> direct_accessions;
    
    // Define command line options
    static struct option long_options[] = {
        {"input", required_argument, NULL, 'i'},
        {"output", required_argument, NULL, 'o'},
        {"accessions-file", required_argument, NULL, 'a'},
        {"accessions", required_argument, NULL, 'A'},
        {"bins", required_argument, NULL, 'b'},
        {"threads", required_argument, NULL, 't'},
        {"help", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0}
    };
    
    int opt;
    int option_index = 0;
    
    // Parse command line options
    while ((opt = getopt_long(argc, argv, "i:o:a:A:b:t:h", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'i':
                kmc_input_dir = optarg;
                break;
            case 'o':
                output_dir = optarg;
                break;
            case 'a':
                accessions_file = optarg;
                use_accession_file = true;
                break;
            case 'A':
                direct_accessions.push_back(optarg);
                break;
            case 'b':
                num_bins = stoi(optarg);
                break;
            case 't':
                num_threads = stoi(optarg);
                break;
            case 'h':
            case '?':
                cout << "Usage: " << argv[0] << " [OPTIONS]" << endl;
                cout << "Options:" << endl;
                cout << "  -i, --input=DIR         Directory containing KMC databases" << endl;
                cout << "  -o, --output=DIR        Output directory for binned k-mers" << endl;
                cout << "  -a, --accessions-file=FILE File containing accession names, one per line" << endl;
                cout << "  -A, --accessions=ACC    Directly specify an accession (can be used multiple times)" << endl;
                cout << "  -b, --bins=N            Number of bins (default: 50)" << endl;
                cout << "  -t, --threads=N         Number of threads (default: all available)" << endl;
                cout << "  -h, --help              Show this help message" << endl;
                return opt == 'h' ? 0 : 1;
        }
    }
    
    // Check required options
    if (kmc_input_dir.empty() || output_dir.empty() || 
        (!use_accession_file && direct_accessions.empty())) {
        cerr << "Error: Missing required options" << endl;
        cerr << "Use --help for usage information" << endl;
        return 1;
    }
    
    // Add trailing slash to directories if not present
    if (kmc_input_dir.back() != '/') kmc_input_dir += '/';
    if (output_dir.back() != '/') output_dir += '/';
    
    // Create output directory if it doesn't exist
    if (!filesystem::exists(output_dir)) {
        cout << "Creating output directory: " << output_dir << endl;
        filesystem::create_directory(output_dir);
    }
    
    // Get accessions from file if specified
    vector<string> accessions;
    if (use_accession_file) {
        accessions = read_accessions_file(accessions_file);
    } else {
        accessions = direct_accessions;
    }
    
    // Prepare KMC database paths
    vector<string> kmc_db_paths;
    for (const auto& acc : accessions) {
        string db_path = kmc_input_dir + acc;
        kmc_db_paths.push_back(db_path);
    }
    
    // Use all available threads if not specified
    if (num_threads == 0) {
        num_threads = thread::hardware_concurrency();
    }
    
    try {
        cout << "Starting k-mer binning using KMC databases" << endl;
        cout << "Number of accessions: " << accessions.size() << endl;
        cout << "Number of bins: " << num_bins << endl;
        cout << "Number of threads: " << num_threads << endl;
        
        // Process the KMC databases in parallel
        process_multiple_databases(kmc_db_paths, accessions, output_dir, num_bins, num_threads);
        
        cout << "All done! Binned k-mers are available in " << output_dir << endl;
    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}

