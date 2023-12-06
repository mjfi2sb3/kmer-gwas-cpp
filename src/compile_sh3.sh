#!/bin/bash
cmake -DCMAKE_BUILD_TYPE=Release -S/scratch/bougous/kmer-gwas-cpp/src -B/scratch/bougous/kmer-gwas-cpp/build_sh3
cmake --build /scratch/bougous/kmer-gwas-cpp/build_sh3 --config Release --target all
