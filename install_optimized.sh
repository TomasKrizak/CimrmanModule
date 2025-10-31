#!/bin/bash
set -euo pipefail

# Clean build directory
rm -rf ./build
mkdir build
cd build

cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-O3 -fno-math-errno"

# Build using many jobs
make -j 4
