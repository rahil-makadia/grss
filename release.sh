#!/bin/bash
rm -rf build dist grss.egg-info
pip uninstall grss -y

pip install .

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
