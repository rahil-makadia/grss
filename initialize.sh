#!/bin/bash

# get cspice
cd ./extern/
python3 get_cspice.py

# initialize and update git submodules
cd ../
git submodule init
git submodule update

# get debiasing data and kernels
cd ./grss/debias/
python3 get_debiasing_data.py
cd ../kernels/
if [ "$1" = "--no-tm-overwrite" ]; then
    python3 get_kernels.py --no-tm-overwrite
else
    python3 get_kernels.py
fi

# return to root
cd ../..
