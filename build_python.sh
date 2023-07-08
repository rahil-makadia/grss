#!/bin/bash
rm -rf build dist grss.egg-info

python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build twine cibuildwheel

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=Mac;;
    CYGWIN*)    machine=Cygwin;;
    MINGW*)     machine=MinGw;;
    MSYS_NT*)   machine=Git;;
    *)          machine="UNKNOWN:${unameOut}"
esac

python3 -m build --sdist --outdir dist
if [ $machine = "Linux" ]; then
    python3 -m cibuildwheel --output-dir dist
fi
if [ $machine = "Mac" ]; then
    python3 -m build --wheel --outdir dist
fi

# if file run with -pypi then upload to pypi, else upload to testpypi
twine check dist/*
