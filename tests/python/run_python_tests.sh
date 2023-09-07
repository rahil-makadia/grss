#!/bin/bash
rm -rf .pytest_cache
rm -rf fit/__pycache__
rm -rf prop/__pycache__
jupyter nbconvert --to script fit/*.ipynb
jupyter nbconvert --to script prop/*.ipynb
clear
pytest
rm -rf fit/*.py
rm -rf prop/*.py