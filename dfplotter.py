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

if len(sys.argv) < 3 or len(sys.argv) > 5:
    print("Usage is {0} filename [xcol] ycol [yerr]".format(sys.argv[0]))
    print("Need x-axis to have yerrors")
    exit()

filename = sys.argv[1]

df = pd.read_csv(filename, delim_whitespace=True, header=None, comment="#")

if len(sys.argv) < 4:
    ycol = int(sys.argv[2])
    df.plot(y=ycol,marker='o',ls='',legend=False)

else:
    xcol = int(sys.argv[2])

    ycol = int(sys.argv[3])

    if len(sys.argv) < 5:
        df.plot(x=xcol,y=ycol,marker='o',ls='',legend=False)

    else:
        print("printing error bars")
        yerr = int(sys.argv[4])

        df.plot(x=xcol,y=ycol,yerr=yerr,marker='s',ls='',mfc='w',capsize=5, ms=5,legend=False)


# print(df,"\n",df.dtypes)

plt.grid(linestyle=':')
plt.xlabel('')

plt.show()
