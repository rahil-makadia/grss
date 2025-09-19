#!/bin/bash
# raise error on failure
set -e

# raise error if no arguments are provided
if [ $# -eq 0 ]; then
    echo "No arguments provided. Please provide the name of the directory containing the test files."
    exit 1
fi

# first argument is the path to the directory containing the test files
cd $1

# for each file in the directory, run the test
jupyter nbconvert --to script *.ipynb
for file in *.py
do
    python $file
done
rm *.py

# return to the original directory
cd ..
