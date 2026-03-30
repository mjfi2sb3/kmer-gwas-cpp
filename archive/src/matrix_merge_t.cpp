#include "thread_pool.hpp"
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
#include <condition_variable>
#include <thread>
#include <queue>

using namespace std;
std::mutex mtx;
std::condition_variable cv;
std::queue<std::function<void()>> tasks;

void worker_thread() {
    while (true) {
        std::function<void()> task;
        {
            std::unique_lock<std::mutex> lock(mtx);
            cv.wait(lock, []{ return !tasks.empty(); });
            task = std::move(tasks.front());
            tasks.pop();
        }
        task();
    }
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

vector<string> get_accessions(string accessions_path)
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
}

void process_file(const string& file_path, unordered_map<string, vector<ushort>>& matrix_, size_t acc_index, size_t NUM_ACC) {
    if (!filesystem::exists(file_path)) {
        std::lock_guard<std::mutex> lock(mtx);
        cout << "file: " << file_path << " doe not exist!" << endl;
    }
    ifstream stream(file_path);
    string line;
    while (getline(stream, line))
    {
        auto record = split(line, '\t');
        if (record.size() != 2)
            throw runtime_error("Could not open kmer count files!");
        string key = record[0];
        ushort value = stoi(record[1]);
        {
            std::lock_guard<std::mutex> lock(mtx);
            if(matrix_.find(key) == matrix_.end())
                matrix_[key] = vector<ushort>(NUM_ACC+1);
            matrix_[key][acc_index] = value;
            matrix_[key][NUM_ACC]++;
        }
    }
    stream.close();
}

void merge_chunk(const uint file_index, const uint min_occur, string input_path, string accessions_path)
{ 
    auto accessions = get_accessions(accessions_path);
    size_t NUM_ACC = accessions.size();
    unordered_map<string, vector<ushort>> matrix_;
    // unordered_map<string, ushort[NUM_ACC+1]> matrix_;
    size_t acc_index = 0;	


    std::vector<std::thread> threads;
    // size_t num_threads = std::min(accessions.size(), std::thread::hardware_concurrency());

    // size_t num_threads = std::min(accessions.size(), static_cast<size_t>(std::thread::hardware_concurrency()));
    size_t num_threads = 4;

    for (size_t i = 0; i < num_threads; ++i) {
        threads.emplace_back(std::thread(worker_thread));
    }
    for (auto &acc : accessions)
    {
        string file_path = input_path + acc + "/" + to_string(file_index) + "_nr.tsv";
        {
            std::lock_guard<std::mutex> lock(mtx);
            tasks.emplace([file_path, &matrix_, acc_index, NUM_ACC]{
                process_file(file_path, matrix_, acc_index, NUM_ACC);
            });
        }
        cv.notify_one();
        acc_index++;
    }
    for (auto& thread : threads) {
        thread.join();
    }

    // size_t num_accessions = accessions.size();
    std::cout << "Number of accessions: " << NUM_ACC << std::endl;
    ofstream m_stream("matrix_acc"+to_string(NUM_ACC)+"_"+to_string(min_occur)+"/" + to_string(file_index) + "_m.tsv");
    ofstream ck_stream("matrix_acc"+to_string(NUM_ACC)+"_"+to_string(min_occur)+"/" + to_string(file_index) + "_c.txt");
    for (auto & pair_ : matrix_) {
    if (pair_.second[NUM_ACC] == NUM_ACC ){
        // std::cout << "kmer: " << pair_.first << std::endl;
        ck_stream << pair_.first << "\n";
    }
	if (pair_.second[NUM_ACC] < min_occur || pair_.second[NUM_ACC] > NUM_ACC - min_occur ) continue;
        m_stream << pair_.first << "\t";
        auto freqs = pair_.second;
        for (int i = 0; i < NUM_ACC ; i++) {
            auto freq = freqs[i];
            if (freq > 0) freq = 1;
            m_stream << freq;
        }
            
        m_stream << "\n";    
    }
     
    m_stream.close();
    ck_stream.close();
}



// void merge_chunk(const uint file_index, const uint min_occur, string input_path, string accessions_path)
// { 

//     auto accessions = get_accessions(accessions_path);
//     size_t NUM_ACC = accessions.size();
//     unordered_map<string, vector<ushort>> matrix_;
//     // unordered_map<string, ushort[NUM_ACC+1]> matrix_;
//     size_t acc_index = 0;	

//     for (auto &acc : accessions)
//     {
//         string file_path = input_path + acc + "/" + to_string(file_index) + "_nr.tsv";
// 	if (!filesystem::exists(file_path)) {cout << "file: " << file_path << " doe not exist!" << endl;}
//         ifstream stream(file_path);
//         string line;
//         while (getline(stream, line))
//         {
//             auto record = split(line, '\t');
//             if (record.size() != 2)
//                 throw runtime_error("Could not open kmer count files!");
//             string key = record[0];
//             ushort value = stoi(record[1]);
//             if(matrix_.find(key) == matrix_.end())
//                 matrix_[key] = vector<ushort>(NUM_ACC+1);
//             matrix_[key][acc_index] = value;
//             matrix_[key][NUM_ACC]++;
	    	
//         }
//         stream.close();
//         acc_index++;
//     }

//     /* --frequency version - enable when needed!
//     ofstream m_stream("matrix/" + to_string(file_index) + "_m.tsv");
//     for (auto & pair_ : matrix_) {
// 	if (pair_.second[NUM_ACC] < min_occur || pair_.second[NUM_ACC] > NUM_ACC - min_occur ) continue;
//         m_stream << pair_.first;
//         for (auto & freq : pair_.second)
//             m_stream << "\t" << freq;
//         m_stream << "\n";    
//     }*/

//     // size_t num_accessions = accessions.size();
//     std::cout << "Number of accessions: " << NUM_ACC << std::endl;
//     ofstream m_stream("matrix_acc"+to_string(NUM_ACC)+"_"+to_string(min_occur)+"/" + to_string(file_index) + "_m.tsv");
//     ofstream ck_stream("matrix_acc"+to_string(NUM_ACC)+"_"+to_string(min_occur)+"/" + to_string(file_index) + "_c.txt");
//     for (auto & pair_ : matrix_) {
//     if (pair_.second[NUM_ACC] == NUM_ACC ){
//         // std::cout << "kmer: " << pair_.first << std::endl;
//         ck_stream << pair_.first << "\n";
//     }
// 	if (pair_.second[NUM_ACC] < min_occur || pair_.second[NUM_ACC] > NUM_ACC - min_occur ) continue;
//         m_stream << pair_.first << "\t";
//         auto freqs = pair_.second;
//         for (int i = 0; i < NUM_ACC ; i++) {
//             auto freq = freqs[i];
//             if (freq > 0) freq = 1;
//             m_stream << freq;
//         }
            
//         m_stream << "\n";    
//     }
     
//     m_stream.close();
//     ck_stream.close();
// }

int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        cout << "usage: " << argv[0] << " <input path> <accessions path> <file index> <min occurence threshold (across panel)>\n";
        return -1;
    }

    try
    {
      	string input_path = argv[1];
	string accessions_path = argv[2];
        uint file_index = stoi(argv[3]);
	uint min_occur = stoi(argv[4]);
      	if (input_path[input_path.size()-1] != '/')
                input_path += '/';
    auto accessions = get_accessions(accessions_path);
    size_t NUM_ACC = accessions.size();
	if (!filesystem::exists("matrix_acc"+to_string(NUM_ACC)+"_"+to_string(min_occur)))
      		filesystem::create_directory("./matrix_acc"+to_string(NUM_ACC)+"_"+to_string(min_occur)+"/");

        cout << "***************************** " << endl;
        cout << "PROCESSING MATRIX CHUNK: " << file_index << endl;
        auto start = chrono::steady_clock::now();

        std::vector<std::thread> threads;
        

        merge_chunk(file_index, min_occur, input_path, accessions_path);
        // merge_chunk(file_index, min_occur, input_path, accessions_path);

        auto end = chrono::steady_clock::now();
        cout << "processing index: " << file_index << "took "
             << chrono::duration_cast<chrono::seconds>(end - start).count()
             << " sec" << endl;

        cout << "FINISHED MATRIX CHUNK: " << file_index << endl;
        cout << "***************************** " << endl;
    }

    catch (exception const &e)
    {
        cout << e.what() << endl;
    }
}

