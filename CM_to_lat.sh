#!/bin/bash

# Create a directory with the lattice energies from the amplitude
# This is a quick way to compare lattice prin_corr fits directly to amplitude

# Check the number of arguments
if [[ $# -lt 1 ]]; then
    echo "Usage is: $0 finite_volume_spectrum.output(s).spectrum"
    exit 1
fi


[[ ! -d Latt_spec ]] && mkdir Latt_spec

for file in $@; do
    awk '{print $1, $2,$3}' $file > /tmp/spec.tmp
    mom=`echo $file | sed 's/d//g' | sed 's/_/ /g' | awk '{print $2}'`
    CM_to_lat_spectrum.py /tmp/spec.tmp $mom | tail +2 > Latt_spec/${file}_lat
done

rm /tmp/spec.tmp

