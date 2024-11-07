#!/bin/bash
# make error stop the script
set -e

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
if [[ $docs -eq 1 ]]; then
    DOCFLAG="-DBUILD_DOCS=ON"
else
    DOCFLAG="-DBUILD_DOCS=OFF"
fi
cmake $DOCFLAG -Dpybind11_DIR=$(pybind11-config --cmakedir) -DPYTHON_EXECUTABLE=$(which python) -B./build -S.
cd build
make
cd ..
