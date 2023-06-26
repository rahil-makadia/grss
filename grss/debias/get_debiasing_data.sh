#!/bin/bash

# get debiasing data files from the JPL Solar System Dynamics Group's FTP server
wget --no-clobber ftp://ssd.jpl.nasa.gov/pub/ssd/debias/debias_2018.tgz -O lowres_data.tgz
wget --no-clobber ftp://ssd.jpl.nasa.gov/pub/ssd/debias/debias_hires2018.tgz -O hires_data.tgz

# extract the data files
# create the directories if they don't exist, delete the files if they do
mkdir -p lowres_data
mkdir -p hires_data
rm -f lowres_data/*
rm -f hires_data/*
tar -xzvf lowres_data.tgz -C lowres_data
tar -xzvf hires_data.tgz -C hires_data
