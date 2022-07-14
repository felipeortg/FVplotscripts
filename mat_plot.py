#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2022-07-14 13:43:18
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 1
# Plot a matrix from a dat file


import sys
import pandas as pd
import iminuitwJK as mJK

if len(sys.argv) != 2:
    print("Usage is: " + sys.argv[0] + " mat_file")
    sys.exit(1)


# -----------------
mat_file = str(sys.argv[1])


df = pd.read_csv(mat_file, delim_whitespace=True, header=None, comment="#")

nparray = df.to_numpy()

xdata = [int(val) for val in nparray[:,0]]
corr = nparray[:,1:]

axs = mJK.plt.subplots()[1]

cmap = mJK.mpl.cm.RdBu_r
norm = mJK.mpl.colors.Normalize(vmin=-1, vmax=1)

mJK.matrix_plot(axs, xdata, corr, cmap=cmap, norm=norm)


mJK.plt.show()

