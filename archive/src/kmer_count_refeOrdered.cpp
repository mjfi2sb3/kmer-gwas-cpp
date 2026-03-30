#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>  // Added this header
#include <map>
#include <string>
#include <bitset>
#include <algorithm>
#include <mutex>
#include <thread>
#include <filesystem>
#include <fstream>
#include <functional>

using namespace std;

// Compile-time k-mer length
const int k = 51;

// Bit encoding for DNA bases
const map<char, string> bits = {
    {'A', "00"}, {'C', "01"}, {'G', "10"}, {'T', "11"}
};

const map<string, char> bases = {
    {"00", 'A'}, {"01", 'C'}, {"10", 'G'}, {"11", 'T'}
};

// Existing bit encoding functions remain the same
bitset<2*k> bit_encode(const string& kmer) {
    string code;
    for (char c : kmer) 
        code += bits.at(c);
    return bitset<2*k>(code);
}

string bit_decode(const bitset<2*k>& code) {
    string kmer;
    string bcode = code.to_string();
    for (size_t i = 0; i < bcode.size(); i += 2)
        kmer += bases.at(bcode.substr(i, 2));
    return kmer;
}

string canonical_kmer(const string& seq) {
    string complement;
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        switch (*it) {
            case 'A': complement += 'T'; break;
            case 'C': complement += 'G'; break;
            case 'G': complement += 'C'; break;
            case 'T': complement += 'A'; break;
            default: return ""; // Invalid base
        }
    }
    return (complement < seq) ? complement : seq;
}

class ReferenceGenomeKmerOrganizer {
public:
    // Enhanced genomic location structure
    struct KmerLocation {
        string chromosome;
        int64_t position;
        string accession;  // Added to track specific accession
        
        bool operator<(const KmerLocation& other) const {
            if (chromosome != other.chromosome)
                return chromosome < other.chromosome;
            return position < other.position;
        }
    };

    // Genomic interval structure
    struct GenomicInterval {
        string chromosome;
        int64_t start;
        int64_t end;
        vector<bitset<2*k>> kmers;
        unordered_set<string> accessions;

        bool contains(const KmerLocation& loc) const {
            return loc.chromosome == chromosome && 
                   loc.position >= start && 
                   loc.position < end;
        }
    };

private:
    // Map of k-mers to their genomic locations across all accessions
    unordered_map<bitset<2*k>, vector<KmerLocation>> reference_kmer_map;

public:
    // Build reference index with accession information
    void build_reference_index(const string& reference_path) {
        ifstream ref_file(reference_path);
        if (!ref_file) {
            throw runtime_error("Cannot open reference genome file");
        }

        string line, chromosome, sequence, current_accession;
        while (getline(ref_file, line)) {
            if (line[0] == '>') {
                // Extract accession from header, assuming FASTA format
                current_accession = line.substr(1);
                size_t first_space = current_accession.find(' ');
                if (first_space != string::npos) {
                    current_accession = current_accession.substr(0, first_space);
                }

                // Process previous chromosome
                if (!sequence.empty()) {
                    index_chromosome_sequence(chromosome, sequence, current_accession);
                }

                // Parse chromosome name
                chromosome = current_accession;
                sequence.clear();
            } else {
                sequence += line;
            }
        }

        // Index last chromosome
        if (!sequence.empty()) {
            index_chromosome_sequence(chromosome, sequence, current_accession);
        }
    }

private:
    void index_chromosome_sequence(const string& chromosome, const string& sequence, const string& accession) {
        for (size_t i = 0; i <= sequence.length() - k; ++i) {
            string kmer = sequence.substr(i, k);
            string canon_kmer = canonical_kmer(kmer);
            
            KmerLocation location{chromosome, static_cast<int64_t>(i), accession};
            
            bitset<2*k> encoded_kmer = bit_encode(canon_kmer);
            reference_kmer_map[encoded_kmer].push_back(location);
        }
    }

public:
    // Generate genomically ordered bins with contiguous intervals
    vector<GenomicInterval> generate_genomically_ordered_bins(int num_bins) {
        vector<GenomicInterval> intervals;
        
        // Group k-mers by chromosome
        map<string, vector<KmerLocation>> chromosome_locations;
        for (const auto& [kmer, locations] : reference_kmer_map) {
            for (const auto& loc : locations) {
                chromosome_locations[loc.chromosome].push_back(loc);
            }
        }

        // Process each chromosome separately
        for (auto& [chromosome, locations] : chromosome_locations) {
            // Sort locations within chromosome
            sort(locations.begin(), locations.end());

            // Calculate interval size for this chromosome
            int64_t chromosome_length = locations.back().position;
            int64_t interval_size = chromosome_length / num_bins;

            // Create bins for this chromosome
            for (int bin = 0; bin < num_bins; ++bin) {
                int64_t start = bin * interval_size;
                int64_t end = (bin == num_bins - 1) ? chromosome_length : (bin + 1) * interval_size;

                GenomicInterval interval;
                interval.chromosome = chromosome;
                interval.start = start;
                interval.end = end;

                // Add k-mers and accessions in this interval
                for (const auto& [kmer, locs] : reference_kmer_map) {
                    for (const auto& loc : locs) {
                        if (loc.chromosome == chromosome && 
                            loc.position >= start && 
                            loc.position < end) {
                            interval.kmers.push_back(kmer);
                            interval.accessions.insert(loc.accession);
                        }
                    }
                }

                if (!interval.kmers.empty()) {
                    intervals.push_back(interval);
                }
            }
        }

        return intervals;
    }

    // Utility method to output bin information
    void output_bins(const vector<GenomicInterval>& intervals, const string& output_path) {
        for (size_t i = 0; i < intervals.size(); ++i) {
            ofstream bin_file(output_path + "/interval_" + to_string(i) + ".txt");
            
            bin_file << "Chromosome: " << intervals[i].chromosome << "\n";
            bin_file << "Start: " << intervals[i].start << "\n";
            bin_file << "End: " << intervals[i].end << "\n";
            bin_file << "Accessions: ";
            for (const auto& acc : intervals[i].accessions) {
                bin_file << acc << " ";
            }
            bin_file << "\n\nK-mers:\n";

            for (const auto& kmer : intervals[i].kmers) {
                bin_file << bit_decode(kmer) << "\n";
            }
        }
    }
};

int main(int argc, char* argv[]) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] 
             << " <reference_genome> <num_bins> <output_path>" << endl;
        return 1;
    }

    try {
        string reference_genome_path = argv[1];
        int num_bins = stoi(argv[2]);
        string output_path = argv[3];

        // Create output directory
        filesystem::create_directories(output_path);

        // Initialize reference genome organizer
        ReferenceGenomeKmerOrganizer ref_organizer;

        // Build reference genome index
        ref_organizer.build_reference_index(reference_genome_path);

        // Generate genomically ordered intervals
        auto ordered_intervals = ref_organizer.generate_genomically_ordered_bins(num_bins);

        // Output intervals to files
        ref_organizer.output_bins(ordered_intervals, output_path);

    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}
