#!/bin/bash

mkdir -p temp/lib temp/bin
cp src/*.cpp src/*.hpp temp/
g++-10 -std=c++20 -O3 -DNDEBUG   -pthread -o temp/lib/mmap_io.cpp.o -c temp/mmap_io.cpp
ar rcs temp/lib/libgwas_tools.a temp/lib/mmap_io.cpp.o
g++-10 -std=c++20 -O3 -DNDEBUG   -pthread -o temp/kmer_count.cpp.o -c temp/kmer_count.cpp
g++-10 temp/kmer_count.cpp.o -Ltemp/lib -lgwas_tools -pthread -o temp/bin/kmer_count

#remove -O3 and replace -DNDEBUG for a debug version ! thanks and good night.

