#!/bin/bash
#cmake -DCMAKE_BUILD_TYPE=Release -S/scratch/bougous/kmer-gwas-cpp/src -B/scratch/bougous/kmer-gwas-cpp/build_sh3
#cmake --build /scratch/bougous/kmer-gwas-cpp/build_sh3 --config Release --target all


g++ -std=c++17 -O3 -march=native -pthread -o kmer_counter_k51 kmer_count.cpp mmap_io.cpp

g++ -std=c++17 -O3 -march=native -pthread -o matrix_merge matrix_merge.cpp mmap_io.cpp
