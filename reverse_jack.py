#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2023-05-01 09:50:50
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 1.0
# Reverse time order of jack file

import sys
import iminuitwJK as mJK

if len(sys.argv) < 3 :
    print(f"Usage is: {sys.argv[0]} old_jack new_jack")
    exit(1)


old_jack = str(sys.argv[1])
new_jack = str(sys.argv[2])

# cfgs, tl, comp
type_jack = mJK.type_Jack_file(old_jack)[2]

if type_jack == 0:
    correl = mJK.get_data(old_jack)[3]

    new_correl = correl[:,::-1]

    mJK.write_Jack_file(new_jack, new_correl)

elif type_jack == 1:
    print("Reading complex correlator, real and imaginary parts")
    correl_re = mJK.get_data(old_jack, 1)[3]
    correl_im = mJK.get_data(old_jack, 2)[3]

    new_correl = correl_re[:,::-1] +1j*correl_im[:,::-1]

    mJK.write_Jack_file(new_jack, new_correl, 1)
