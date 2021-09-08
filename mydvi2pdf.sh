#!/bin/bash -v

# Fix ghostcript 9.5 by hand since I use it for latexit + feynmp

# Check the number of arguments
if [ "$#" -lt 1 ]; then
    echo "Usage is: $0 dvifile"
    exit 1
fi

input=$1
inputwoext=${input%.dvi}

dvips $input

ps2pdf14 -dPDFSETTINGS=/prepress $inputwoext.ps $inputwoext.pdf

#rm $inputwoext.ps # Maybe is not a good idea to delete the old file in case something goes wrong