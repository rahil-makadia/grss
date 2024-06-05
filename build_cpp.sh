#!/bin/bash

clean=0
docs=0
for i in $*; do
    if [[ $i == "-clean" ]]; then
        clean=1
    fi
    if [[ $i == "-docs" ]]; then
        docs=1
    fi
done

# This script is used to build the project.
if [[ $clean -eq 1 ]]; then
    rm -rf build
fi
mkdir -p build
cd build
if [[ $docs -eq 1 ]]; then
    cmake -DBUILD_DOCS=ON ..
else
    cmake -DBUILD_DOCS=OFF ..
fi
make
cd ..
