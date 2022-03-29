#!/bin/bash

# Plot gnu_files in directory

for i in *gnu; do echo $i; done
echo "Plot? (y/n)"
read ans
if [[ $ans = y* ]]; then
    echo "Plotting..."
    for i in *gnu; do gnuplot $i; done
fi

exit 0
