#!/bin/bash
rm -rf build dist grss.egg-info

python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build twine cibuildwheel

python3 -m build --sdist --outdir dist
python3 -m cibuildwheel --output-dir dist

# if file run with -pypi then upload to pypi, else upload to testpypi
twine check dist/*
