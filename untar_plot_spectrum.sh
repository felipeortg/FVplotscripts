#!/bin/bash

# Untar and plot the spectrum for various t0 folders

# No arguments is for tarfiles
if [ "$#" -lt 1 ]; then
    echo "Untaring"

    for tarfile in `ls *.tar`
    do
        tar -xf $tarfile  &&  rm $tarfile
        folder=`echo $tarfile | awk -F. '$0=$1'`
	    echo $folder
        cd $folder
        plotting_spectrum.py out_fit_t0*
        cd ..
    done

    exit 0
fi

#any arguments assumes just plotting
echo "Just plotting"
for folder in `ls`
    do
        echo $folder
        cd $folder
        plotting_spectrum.py out_fit_t0*
        cd ..
    done

exit 0