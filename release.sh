#!/bin/bash
rm -rf build dist grss.egg-info
pip uninstall grss -y

pip install .

python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build twine

python3 -m build > release.log

python3 -m twine upload --repository testpypi dist/*
