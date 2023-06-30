#!/bin/bash

# This script is used to generate the documentation for the project.
pip install .
cd docs
rm -rf source/generated
sphinx-apidoc -o source ../grss
make clean html
# make epub
# make latexpdf
cd ..
