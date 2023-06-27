#!/bin/bash

conda env create -f environment.yml
conda activate grss

# initialize and update git submodules
cd ../pybind11/
git submodule init
cd ../../
git submodule update

# download debiasing data
cd ./grss/debias/
./get_debiasing_data.sh

# download/update earth orientation kernels
cd ./kernels/
./update_kernels.sh

# build the python wrapper
cd ..
./build.sh -clean
