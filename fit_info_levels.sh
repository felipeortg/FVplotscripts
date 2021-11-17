#!/bin/bash

# Get info in each fit
# Print t0, size of matrix, size if SVD cuts, chi2 per dof
# t0, oi, of, x2/dof
# print also the energy (and error) of the first n levels 
# Check the number of arguments
if [ "$#" -lt 1 ]; then
    echo "Usage is: $0 n-levels"
    exit 1
fi

echo "Run in the irrep with all the fits"
irrep=${PWD##*/}
echo "Loping over " `ls` 
echo " ------ "
printf "  %5s %2s %2s %6s" "t0" opsi opsf "chi2/dof"
for n in `seq $1`
do
    printf "%5s" $n
done 
printf "\n"


for fit in `ls`
do
    cd $fit

    echo $fit
    for t0 in `ls -F | grep /`
    do

	cd $t0

        t0fol=${PWD##*/}
        t0num=`echo $t0fol | sed 's/\_/\ /g' | awk '{print $3}'`
        
        opsi=`grep "${irrep}P " fit_t0*| wc -l`
        let opsf=`wc -l out* | awk '{print $1}'`-15

        chi2=`grep "best" out* | awk '{print $5}'`

        printf "  %5s %2s %2s %5.2f" "t0_$t0num" $opsi $opsf $chi2

	let nm1=$1-1
	for n in `seq 0 $nm1`
 	do
	    enerde=`awk '{if ($1 == "0|") print $3$4}' out*`
 	    printf "  %17s" $enerde
 	done
 	printf "\n"

        # echo "    t0_"$t0num": oi" $opsi" of" $opsf" x2/dof" $chi2
        cd ..
    done 
    
    cd ..
done

