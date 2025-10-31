#!/bin/bash
set -e

rm -rf ./build
mkdir build
cd build

# Configure with optimizations + debug symbols (good for profiling)
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo

# Build
make -j 4
