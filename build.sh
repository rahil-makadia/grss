#!/bin/sh

# This script is used to build the project.
rm -r build
mkdir build
cd build
cmake ..
make
