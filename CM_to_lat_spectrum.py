#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2022-03-28 20:31:19
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 1
# Calculate the lattice spectrum from an scattering_devel CM spectrum output


import numpy as np
import csv
import sys

# These are in the same git repo, and therefore no need to append PATH
try:
    import spectrum as spec
    import lattice
except:
    print('Libraries spectrum and/or lattice not found')
    raise Exception('This script is part of a repo and needs the other libraries in the repo to work...')

# Import the name of this lattice_channel
sys.path.append('../../')
try:
    import this_lattice as tl
except:
    raise Exception('Run in amplitudes/amplitude/FV_spec, with file amplitudes/this_lattice.py')

if len(sys.argv) < 3:
    raise Exception('Usage is {0} spectrum_file P_latt'.format(sys.argv[0]))

spectrum_file = str(sys.argv[1])
vec = spec.label2vec(sys.argv[2])

print("Lattice momentum: ", vec)

lattice.init_lattice_props(tl.lattice)
chi = lattice.chi


with open(spectrum_file, 'r') as f:
    file = csv.reader(f, delimiter=' ')
    for line in file:
        L = int(line[0])
        bp, ap = line[1].split('.')
        sd = len(ap)
        bp, ap = line[2].split('.')
        sde = len(ap)
        energy = float(line[1])
        unc = float(line[2])

        en, unc = spec.latt_prop_unc([energy,unc] ,vec , L*chi)

        print("{0} {1:.{3}f} {2:.{4}f}".format(L, en, unc, sd, sde))



