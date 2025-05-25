#!/bin/bash
#/usr/bin/cmake --no-warn-unused-cli -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE -DCMAKE_BUILD_TYPE:STRING=Release -DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc-10 -DCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++-10 -S/project/src -B/project/build -G Ninja
#/usr/bin/cmake --build /project/build --config Release --target all
## compile matrix_merge.cpp alone

# compile matrix_merge
g++ -std=c++17 -o matrix_merge_noThread_v2 matrix_merge.cpp -lstdc++fs

# compile kmc_kmer_binning_ori
g++ -std=c++17 -O3 -o kmc_kmer_binning_ori kmc_kmer_binning_ori.cpp  /ibex/sw/rl9c/kmc/3.2.1/rl9_conda3/KMC/kmc_api/kmer_api.cpp /ibex/sw/rl9c/kmc/3.2.1/rl9_conda3/KMC/kmc_api/kmc_file.cpp /ibex/sw/rl9c/kmc/3.2.1/rl9_conda3/KMC/kmc_api/mmer.cpp -I. -pthread

