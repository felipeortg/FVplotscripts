#!/bin/bash

# Plot an ensemble jackknife file with errorbars
# Use calc or calcbc and dfplotter
# Check for dfplotter
if ! command -v dfplotter.py &> /dev/null; then
    echo "dfplotter.py not installed in path..."
    exit 1
fi

# Check for calc
if ! command -v calc &> /dev/null; then
    echo "calc not installed in path..."
    exit 1
fi

# Check for calcbc
if ! command -v calcbc &> /dev/null; then
    echo "calcbc not installed in path..."
    exit 1
fi

# Check the number of arguments
if [[ ($# -lt 1) || ($# -gt 2) ]]; then
    echo "Usage is: $0 corr.jack [imag_part]"
    exit 1
fi

correl=$1

if [[ $# -gt 1 ]]; then
    realpart=$2
else
    realpart=0
fi

comple=`head -1 $correl | awk '{print $3}'`

if [[ ( $comple -eq 0 ) && ($realpart -eq 0) ]]; then
    calc $correl | awk '{print $1,$2,$3}' > /tmp/mean.jack
elif [[ ($comple -eq 1 ) && ($realpart -eq 0) ]]; then
    echo "plot real part"
    calcbc "real( $correl )" | awk '{print $1,$2,$3}' > /tmp/mean.jack
elif [[ ($comple -eq 1 ) && ($realpart -eq 1) ]]; then
    echo "plot imag part"
    calcbc "imag( $correl )" | awk '{print $1,$2,$3}' > /tmp/mean.jack
else
    echo "Cannot plot the imaginary part of a real correlator ..."
    exit 1
fi

dfplotter.py /tmp/mean.jack 0 1 2

rm /tmp/mean.jack
