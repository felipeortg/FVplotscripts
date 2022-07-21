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

maxval = np.max(np.abs(mm))

# we normally plot correlations, but in case you plot things that are bigger
# smaller things would not work, but then the purpose of this is nor clear...
if maxval < 1:
    maxval = 1

cmap = mJK.mpl.cm.RdBu_r
norm = mJK.mpl.colors.Normalize(vmin=-maxval, vmax=maxval)

mJK.matrix_plot(axs, xdata, corr, cmap=cmap, norm=norm)


mJK.plt.show()

