#!/bin/bash
rm -rf build dist grss.egg-info

python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build twine

python3 -m build --sdist --wheel

# if file run with -pypi then upload to pypi, else upload to testpypi
twine check dist/*
