#!/bin/bash

clean=0
for i in $*; do
    if [[ $i == "-clean" ]]; then
        clean=1
    fi
done

# This script is used to build the project.
if [[ $clean -eq 1 ]]; then
    rm -rf build
    mkdir build
    cd build
    rm ../grss/*.so
    cmake ..
    make
    cp *.so ../grss/
else
    # if build directory does not exist, create it
    if [[ ! -d "build" ]]; then
        mkdir build
    fi
    cd build
    rm ../grss/*.so
    cmake ..
    make
    cp *.so ../grss/
fi
