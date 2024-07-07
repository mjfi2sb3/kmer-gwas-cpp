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


using namespace std;

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

void merge_chunk(const uint file_index, const uint min_occur, string input_path, string accessions_path, string delimiter, bool show_count, string ouput_dir)
{ 

    auto accessions = get_accessions(accessions_path);
    size_t NUM_ACC = accessions.size();
    unordered_map<string, vector<ushort>> matrix_;
    // unordered_map<string, ushort[NUM_ACC+1]> matrix_;
    size_t acc_index = 0;	

    /*  - CHECK THAT ALL BINS EXIST BEFORE COMMITING TO ANY FILE READING
        - EXIT WITH SYSTEM ERROR IF ANY OF THE FILES DO NOT EXIT
    */
    for (auto &acc : accessions)
    {
        string file_path = input_path + acc + "/" + to_string(file_index) + "_nr.tsv";
	    /*if (!filesystem::exists(file_path)) {cout << "file: " << file_path << " doe not exist!" << endl;}*/
	    if (!filesystem::exists(file_path)) {
            throw runtime_error("File: " + file_path + " does not exist!");
        }
        
     }

    for (auto &acc : accessions)
    {
        string file_path = input_path + acc + "/" + to_string(file_index) + "_nr.tsv";
        
	    /*if (!filesystem::exists(file_path)) {cout << "file: " << file_path << " doe not exist!" << endl;}*/
	    /*if (!filesystem::exists(file_path)) {
            throw runtime_error("File: " + file_path + " does not exist!");
        }*/
	
        ifstream stream(file_path);
        if (!stream) {
            throw runtime_error("Failed to open file: " + file_path);
        }
        
        string line;
        while (getline(stream, line))
        {
            auto record = split(line, '\t');
            if (record.size() != 2)
                throw runtime_error("Could not open kmer count files!");
            string key = record[0];
            ushort value = stoi(record[1]);
            if(matrix_.find(key) == matrix_.end())
                matrix_[key] = vector<ushort>(NUM_ACC+1);
            matrix_[key][acc_index] = value;
            matrix_[key][NUM_ACC]++;
	    	
        }
        stream.close();
        acc_index++;
    }

    /* --frequency version - enable when needed!
    ofstream m_stream("matrix/" + to_string(file_index) + "_m.tsv");
    for (auto & pair_ : matrix_) {
	if (pair_.second[NUM_ACC] < min_occur || pair_.second[NUM_ACC] > NUM_ACC - min_occur ) continue;
        m_stream << pair_.first;
        for (auto & freq : pair_.second)
            m_stream << "\t" << freq;
        m_stream << "\n";    
    }*/

    // size_t num_accessions = accessions.size();
    //std::cout << "Number of accessions: " << NUM_ACC << std::endl;
    /*std::string ouput_dir = "";
    if (show_count)
    {
		ouput_dir = "matrix_acc"+to_string(NUM_ACC)+"_minOcc"+to_string(min_occur)+"_count/";
    }
    else
    {
    	ouput_dir = "matrix_acc"+to_string(NUM_ACC)+"_minOcc"+to_string(min_occur)+"_pres-abs/";
    }
    std:cout << ouput_dir << "\n";
    if (!filesystem::exists("matrix_acc"+to_string(NUM_ACC)+"_"+to_string(min_occur)))
		  		filesystem::create_directory("./matrix_acc"+to_string(NUM_ACC)+"_"+to_string(min_occur)+"/");*/
		  		
    //std::string ouput_dir = "matrix_acc"+to_string(NUM_ACC)+"_"+to_string(min_occur)+"/";
    
    /*ofstream m_stream("matrix_acc"+to_string(NUM_ACC)+"_"+to_string(min_occur)+"/" + to_string(file_index) + "_matrix.tsv");
    ofstream ck_stream("matrix_acc"+to_string(NUM_ACC)+"_"+to_string(min_occur)+"/" + to_string(file_index) + "_core.txt");*/
    
    ofstream m_stream(ouput_dir + to_string(file_index) + "_matrix.tsv");
    ofstream ck_stream(ouput_dir + to_string(file_index) + "_core.txt");
    
    for (auto & pair_ : matrix_) {
      // check if a k-mer occurs in all accessions and therefore is flagged as core k-mer
    if (pair_.second[NUM_ACC] == NUM_ACC ){
        // std::cout << "kmer: " << pair_.first << std::endl;
        ck_stream << pair_.first << "\n";
    }
	// if (pair_.second[NUM_ACC] < min_occur || pair_.second[NUM_ACC] > NUM_ACC - min_occur ) continue;
        m_stream << pair_.first << "\t";
        auto freqs = pair_.second;
        for (int i = 0; i < NUM_ACC ; i++) {
            auto freq = freqs[i];
            if (!show_count && freq > 0) freq = 1;
            m_stream << delimiter << freq;
        }
            
        m_stream << "\n";    
    }
     
    m_stream.close();
    ck_stream.close();
}

int main(int argc, char *argv[])
{
	string input_path, accessions_path;
	uint file_index;
	uint min_occur = 0;
	std::string delimiter = "\t";  // Default value: tab
	bool show_count = true;  // Default is to show counts

	// Define the long options
	static struct option long_options[] = {
		{"input", required_argument, 0, 'i'},
		{"accessions", required_argument, 0, 'a'},
		{"index", required_argument, 0, 'f'},
		{"threshold", optional_argument, 0, 't'},
		{"delimiter", required_argument, 0, 'd'},
		{"count", required_argument, 0, 'c'},
		{0, 0, 0, 0}
	};
    int option_index = 0;
    int c;
    
    while ((c = getopt_long(argc, argv, "i:a:f:t:d:c:", long_options, &option_index)) != -1) 
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
            default:
                abort();
        }
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
		     << "\t\t--threshold OBSOLETE & INACTIVE <min/max occurence threshold> (default: " << min_occur << ")\n"
		     << "\t\t--delimiter <delimiter type: tab|none> (default: " << (delimiter == "\t" ? "tab" : (delimiter == " " ? "space" : "none")) << ")\n"
		     << "\t\t--count <print matrix as absence/presence or actual k-mer counts; type: y|n> (default: " << (show_count ? "y" : "n") << ")\n\n";
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
   		std:string delim = (delimiter == "\t" ? "tab" : (delimiter == " " ? "space" : "none"));
		if (show_count)
		{
			//ouput_dir = "matrix_acc"+to_string(NUM_ACC)+"_minOcc"+to_string(min_occur)+"_count_delim-"+delim+"/";
			ouput_dir = "matrix_acc"+to_string(NUM_ACC)+"_count_delim-"+delim+"/";
		}
		else
		{
			//ouput_dir = "matrix_acc"+to_string(NUM_ACC)+"_minOcc"+to_string(min_occur)+"_pres-abs_delim-"+delim+"/";
			ouput_dir = "matrix_acc"+to_string(NUM_ACC)+"_pres-abs_delim-"+delim+"/";
		}
   		
		if (!filesystem::exists(ouput_dir)){
			filesystem::create_directory(ouput_dir);
		}

        cout << "***************************** " << endl;
        cout << "PROCESSING MATRIX CHUNK: " << file_index << endl;
        cout << "Number of accessions: " << NUM_ACC << std::endl;
        cout << "Output Folder: " << ouput_dir << std::endl;
        
        auto start = chrono::steady_clock::now();

        merge_chunk(file_index, min_occur, input_path, accessions_path, delimiter, show_count, ouput_dir);

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

