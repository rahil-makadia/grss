#!/bin/bash

# get cspice and install pybind11
cd ./extern/
python3 get_cspice.py
pip3 install "pybind11[global]>=2.10.0"

# get debiasing data and kernels
cd ../grss/debias/
python3 get_debiasing_data.py
cd ../kernels/
if [ "$1" = "--tm-overwrite" ]; then
    python3 get_kernels.py
else
    python3 get_kernels.py --no-tm-overwrite
fi

# return to root
cd ../..
