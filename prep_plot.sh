#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "Usage is $0 C_n_m.dat"
    exit 1
fi

file=$1

filewo="${file%.*}"


sed "s/modify/${file}/g" ~/FVplotscripts/ini.gnu > ${filewo}.gnu

gnuplot -persist ${filewo}.gnu
