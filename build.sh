#!/bin/bash

clean=0
for i in $*; do
    if [[ $i == "-clean" ]]; then
        clean=1
    fi
done

# This script is used to build the project.
if [[ $clean -eq 1 ]]; then
    rm -r build
    mkdir build
    cd build
    cmake .. -DPYBIND11_FINDPYTHON=ON
    make
    rm ../examples/*.so
    cp *.so ../examples/
else
    cd build
    make
    rm ../examples/*.so
    cp *.so ../examples/
fi