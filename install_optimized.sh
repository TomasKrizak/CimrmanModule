#!/bin/bash
set -e

# Clean build directory
rm -rf ./build
mkdir build
cd build

# Configure with optimizations and debug symbols for profiling
cmake .. \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_CXX_FLAGS="-march=native -funroll-loops -fomit-frame-pointer"

# Build using all available cores
make -j$(nproc)
