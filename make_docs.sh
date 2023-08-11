#!/bin/bash

# This script is used to generate the documentation for the project.
pip3 install .
cd docs
rm -rf build
rm -rf source/_autosummary
make clean html
make latexpdf
cd ..
