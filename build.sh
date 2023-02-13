#!/bin/zsh

# This script is used to build the project.
rm -r build
mkdir build
cd build
cmake .. -DPYBIND11_FINDPYTHON=ON
make
rm ../examples/*.so
cp *.so ../examples/
