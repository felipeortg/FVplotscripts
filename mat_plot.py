#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2022-07-14 13:43:18
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 1
# Plot a matrix from a dat file


import sys
import pandas as pd
import iminuitwJK as mJK

if len(sys.argv) < 2:
    print("Usage is: " + sys.argv[0] + " mat_file [save_file]")
    sys.exit(1)

if len(sys.argv) == 3:
    mJK.plt.rc('text', usetex=True)
    mJK.plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    mJK.plt.rc('font', family='serif')
    mJK.plt.rc('font', size=11)
    mJK.plt.rc('lines', linewidth=1)


# -----------------
mat_file = str(sys.argv[1])


df = pd.read_csv(mat_file, delim_whitespace=True, header=None, comment="#")

nparray = df.to_numpy()

xdata = [int(val) for val in nparray[:,0]]
corr = nparray[:,1:]

axs = mJK.plt.subplots()[1]

maxval = mJK.np.max(mJK.np.abs(corr))

# we normally plot correlations, but in case you plot things that are bigger
# smaller things would not work, but then the purpose of this is nor clear...
if maxval < 1:
    maxval = 1

cmap = mJK.mpl.cm.RdBu_r
norm = mJK.mpl.colors.Normalize(vmin=-maxval, vmax=maxval)

mJK.matrix_plot(axs, xdata, corr, cmap=cmap, norm=norm)


if len(sys.argv) > 2:
    mJK.plt.savefig(str(sys.argv[2])+'.pdf',transparent=True, bbox_inches='tight')
else:
    mJK.plt.show()

