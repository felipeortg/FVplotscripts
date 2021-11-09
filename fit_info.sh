#!/bin/bash

# Get info in each fit
# Print t0, size of matrix, size if SVD cuts, chi2 per dof
# t0, oi, of, x2/dof

echo "Run in the irrep with all the fits"
irrep=${PWD##*/}
echo "Loping over " `ls` 
echo " ------ "
echo "     oi of x2/dof"

for fit in `ls`
do

    cd $fit
    echo $fit
    for t0 in `ls`
    do
        cd $t0
        t0fol=${PWD##*/}
        t0num=`echo $t0fol | sed 's/\_/\ /g' | awk '{print $3}'`
        
        opsi=`grep "${irrep}P " fit_t0*| wc -l`
        let opsf=`wc -l out* | awk '{print $1}'`-15

        chi2=`grep "best" out* | awk '{print $5}'`

        echo "    "$opsi $opsf $chi2 "t0_$t0num"

        # echo "    t0_"$t0num": oi" $opsi" of" $opsf" x2/dof" $chi2
        cd ..
    done
    cd ..
done