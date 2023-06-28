#!/bin/bash

# initialize and update git submodules
cd ./extern/pybind11/
git submodule init
cd ../../
git submodule update

cd ./grss/debias/
python get_debiasing_data.py

cd ../kernels/
python get_kernels.py
