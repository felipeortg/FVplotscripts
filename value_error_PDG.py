#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2023-04-28 16:14:17
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 1.0
# Pretty print num+error in PDG convention


import iminuitwJK as mJK
import sys

if len(sys.argv) < 3:
    print(f"Usage is {sys.argv[0]} value [+/-] error")
    exit(1)
if sys.argv[2] == '+/-':
    error = float(sys.argv[3])
else:
    error = float(sys.argv[2])

value = float(sys.argv[1])

ve_dict = mJK.value_error_rounding(value, error)

print(mJK.ve_dict2string(ve_dict))  
