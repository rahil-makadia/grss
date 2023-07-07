#!/bin/bash

# get cspice
cd ./extern/
python3 get_cspice.py

# initialize and update git submodules
cd ./pybind11/
git submodule init
cd ../../
git submodule update

# get debiasing data and kernels
cd ./grss/debias/
python3 get_debiasing_data.py
cd ../kernels/
python3 get_kernels.py

# return to root
cd ../..
