#!/bin/bash

# Launch a reading of the scattering_devel spectrum to plot a comparison to the lattice
# This assumes the data is located in 
# ../../spectrum/data/
# There is a file with lattice name in
# ../this_lattice.py
# There is a file with the location of Ecm_ini for this fit in 
# ./dataset
# And that the spectrum from scatdevel is in the folder
# FV_spec

# Check the number of arguments
if [ "$#" -lt 1s ]; then
    echo "Usage is: $0"
    exit 1
fi

echo "Assuming data in ../../spectrum/data/"
echo "Assuming lattice name in ../this_lattice.py"
echo "Assuming location of Ecm_ini in ./dataset"
echo "Assuming interacting spectrum from scatdevel in FV_spec"


cd "FV_spec"
scatfiles=`ls *spectrum`
cd ..

echo "Plotting: " $scatfiles

plot_data_fit.py $scatfiles
