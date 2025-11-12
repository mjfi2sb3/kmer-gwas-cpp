#!/bin/bash

# Load gcc module (version 12.2.0 for C++17 filesystem support)
module load gcc/12.2.0

# Compile matrix_merge.cpp
g++ -std=c++17 -o matrix_merge matrix_merge.cpp

# Copy executable to build_sh3
cp matrix_merge ../build_sh3/

echo "Compilation complete! matrix_merge executable created and copied to build_sh3/"
