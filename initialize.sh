#!/bin/bash

conda env create -f environment.yml
conda activate grss

# initialize and update git submodules
cd ./extern/astrocat_debiasing/
git submodule init
cd ../pybind11/
git submodule init
cd ../../
git submodule update

# TODO: download debiasing data, add a downloader to the astrocat_debiasing repo

# download/update earth orientation kernels
cd ./kernels/
./update_kernels.sh

# build the python wrapper
cd ..
./build.sh -clean
