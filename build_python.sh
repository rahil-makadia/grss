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

if [ $machine = "Mac" ]; then
    # force download arm64 version of cspice
    export ARCHFLAGS="-arch arm64"
    cd ./extern
    python3 get_cspice.py
    cd ..
fi

python3 -m build --sdist --outdir dist
python3 -m cibuildwheel --output-dir dist

# if file run with -pypi then upload to pypi, else upload to testpypi
twine check dist/*
