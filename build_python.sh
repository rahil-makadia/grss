#!/bin/bash
rm -rf build dist grss.egg-info

python3 -m pip install --upgrade pip build twine cibuildwheel

python3 -m build --sdist --outdir dist
python3 -m cibuildwheel --output-dir dist

twine check dist/*
