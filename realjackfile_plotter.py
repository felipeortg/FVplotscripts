#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2024-12-11 17:40:11
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 0.1
# Plot [save plot] real part of correlator


import numpy as np
import matplotlib.pyplot as plt
import sys
import iminuitwJK as mJK
from cycler import cycler


if len(sys.argv) < 2:
    raise Exception("Not enough arguments:\n"+sys.argv[0]+" jackfile [filename_to_save]")

colors_def = [
    [2/255,103/255,253/255,1],
    [254/255,30/255,4/255,1]
    ]  
default_cycler = cycler(color=colors_def) 

if len(sys.argv) == 3:
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    plt.rc('font', family='serif')


plt.rc('font', size=11)
plt.rc('axes.formatter', useoffset = False)
plt.rc('lines', linewidth=1)
plt.rc('axes', prop_cycle=default_cycler)

filename = str(sys.argv[1])

# get the real part, imaginary should be zero
cfgs, npoints, xdata, ydata = mJK.get_data(filename)


f, ax = plt.subplots()#1,1,figsize=(4.5,3))


mJK.plot_ensemble_mean_err(ax, xdata, ydata)

plt.xlabel(r'$t/a_t$')

plt.xlim([min(xdata), max(xdata)])

plt.grid(ls=':')
plt.subplots_adjust(right=0.97)
plt.subplots_adjust(bottom=0.16)


if len(sys.argv) == 3:  
    plt.savefig(str(sys.argv[2])+'.pdf',transparent=True, bbox_inches='tight')
else:
    plt.show()
