#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2022-02-11 15:55:48
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 2

# Give the name of the file
# Give the numbers of columns of the data to plot
# Use pandas

import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
from cycler import cycler

colors_def = [
    [2/255,103/255,253/255,1],
    [254/255,30/255,4/255,1]
    ]  
default_cycler = cycler(color=colors_def) 


plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
plt.rc('font', family='serif')

plt.rc('font', size=11)
plt.rc('axes.formatter', useoffset = False)
plt.rc('lines', linewidth=1)
plt.rc('axes', prop_cycle=default_cycler)

if len(sys.argv) < 5:
    print("Usage is {0} filename xcol ycol yerr name".format(sys.argv[0]))
    exit()

filename = sys.argv[1]

df = pd.read_csv(filename, delim_whitespace=True, header=None, comment="#")

xcol = int(sys.argv[2])

ycol = int(sys.argv[3])

yerr = int(sys.argv[4])

f, ax = plt.subplots(1,1,figsize=(4.5,3))

df.plot(x=xcol,y=ycol,yerr=yerr,ax=ax, marker='s',ls='',mfc='w',capsize=5, ms=5,legend=False)


plt.grid(linestyle=':')

plt.xlabel(r'$t/a_t$')

plt.subplots_adjust(right=0.97)
plt.subplots_adjust(bottom=0.16)

plt.savefig(f"{sys.argv[5]}.pdf", transparent=True, bbox_inches='tight')
