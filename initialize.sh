#!/bin/bash

# initialize and update git submodules
cd ./extern/pybind11/
git submodule init
cd ../../
git submodule update
