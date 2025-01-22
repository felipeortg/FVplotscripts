#!/bin/bash

# Calculate the effective mass of correlation
# Use ensbc for jackknife statistics

# Check for ensbc
if ! command -v ensbc &> /dev/null; then
    echo "ensbc not installed in path..."
    exit 1
fi


# Check the number of arguments
if [[ ($# -lt 2) || ($# -gt 3) ]]; then
    echo "Usage is: $0 corr.jack effmass.jack [dt=3]"
    exit 1
fi

correl=$1
effmass=$2

if [[ $# -gt 2 ]]; then
    dt=$3
else
    dt=3
fi


timetemp=`head -1 $correl | awk '{print $2}'`
tmax=`expr $timetemp - 1 `
tmaxdt=`expr $tmax - $dt `


#echo "$1 = - log( real( extract( $correl , $dt, $tmax ) ) / real( extract( $correl , 0, $tmaxdt ) ) ) / $dt "

ensbc "$effmass = - log( real( extract("'"'${correl}'"'", $dt + 1 , $tmax ) ) / real( extract( "'"'${correl}'"'" , 1 , $tmaxdt ) ) ) / $dt "


