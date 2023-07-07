#!/bin/bash
rm -rf build dist grss.egg-info
pip3 uninstall grss -y

python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build twine

python3 -m build > release.log

# if file run with -pypi then upload to pypi, else upload to testpypi
server=$1
if [ "$server" = "-pypi" ]; then
    python3 -m twine upload dist/*
else
    python3 -m twine upload --repository testpypi dist/*
fi

# install from pypi or testpypi based on -pypi flag
if [ "$server" = "-pypi" ]; then
    python3 -m pip install --upgrade grss
else
    python3 -m pip install --upgrade --index-url https://test.pypi.org/simple/ --no-deps grss
fi
