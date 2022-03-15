#!/bin/bash

# Create directories from a file

# Check the number of arguments
if [[ $# -ne 1 ]]; then
    echo "Usage is: $0 file_w_dirslist"
    exit 1
fi

for i in `cat $1`; do mkdir -p ${i}; done
