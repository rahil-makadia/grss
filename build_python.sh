#!/bin/bash
rm -rf build dist grss.egg-info

python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build twine cibuildwheel

python3 -m build --sdist --outdir dist

twine check dist/*
