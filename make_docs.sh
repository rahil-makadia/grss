#!/bin/bash

# This script is used to generate the documentation for the project.
pip install .
cd docs
rm -rf build
rm -rf source/_autosummary
make clean html
make epub
make latexpdf
cd ..
