#!/bin/bash

# Launch a reading of the scattering_devel spectrum to plot a comparison to the lattice
# This assumes the data is located in 
# ../../spectrum/data/
# This can be modified by adding a variable header into the dataset file
# There is a file with lattice name in
# ../this_lattice.py
# There is a file with the location of Ecm_ini for this fit in 
# ./dataset
# And that the spectrum from scatdevel is in the folder
# FV_spec

# Check the number of arguments
if [ "$#" -gt 0 ]; then
    echo "Usage is: $0"
    exit 1
fi

# cat dataset.py

echo "Assuming all reconfit dataset in ../../spectrum/data_useme/ unless modified in ./dataset"
echo "Assuming lattice name in ../this_lattice.py unless modified in ./dataset"
echo "Assuming interacting spectrum from scatdevel in FV_spec"


print_xpath="/Users/felipeortg/Documents/scattering/install/adat/bin/print_xpath"

if [ ! -f fit_finite_volume_spectrum.ini.xml ]; then
    echo "Run in folder with fit_finite_volume_spectrum.ini.xml"
    exit 1
fi

Ecm_xml=`$print_xpath fit_finite_volume_spectrum.ini.xml "//EcmData"`

# cd "FV_spec"
scatfiles=`ls FV_spec/*spectrum`
# cd ..

echo "Plotting: " $scatfiles

plot_data_fit.py $Ecm_xml $scatfiles
