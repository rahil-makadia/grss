#!/bin/bash

# get pybind11
pip3 install "pybind11[global]>=2.10.0"

# get debiasing data and kernels
cd ../grss/debias/
python3 get_debiasing_data.py
cd ../kernels/
python3 get_kernels.py

# return to root
cd ../..
