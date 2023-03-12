#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2022-07-13 16:31:10
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 1.0
# Create pdf with latex plot of data and correlation

import sys

if len(sys.argv) < 2:
    print("Usage is: " + sys.argv[0] + " jackfile [mask] [output]")
    sys.exit(1)

jackfile = str(sys.argv[1])

mask = []
save = False
outn = 2
if len(sys.argv) > 2:
    try:
        masksints = [int(n) for n in str(sys.argv[2]).split("-")]

        outn += 1

        # when two values given make a range
        if len(masksints) == 2:
            mask = [n for n in range(masksints[0],masksints[1]+1)]

        # otherwise use directly as a list
        else:
            mask = masksints

    except:
    # if that argument value is not numeric means mask is empty
    # this has the drawback that the out files cannot have purely numeric names...
        pass

    # print(len(sys.argv), outn)

    if len(sys.argv) > outn: 

        outfile = str(sys.argv[outn])
        save = True



import iminuitwJK as mJK

if save:
    mJK.plt.rc('text', usetex=True)
    mJK.plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    mJK.plt.rc('font', family='serif')
    mJK.plt.rc('font', size=14)
    mJK.plt.rc('axes.formatter', useoffset = False)





fig, axs = mJK.plt.subplots(num=1)
corr = mJK.corr_mat_from_file(jackfile, axs, mask=mask)



# data_size = mJK.np.shape(corr)[0]


eigvals = mJK.np.linalg.eigvalsh(corr)[::-1]

evals_norm = eigvals/eigvals[0]

fig1, axs1 = mJK.plt.subplots(num=2, figsize=(6,3))

axs1.plot(evals_norm, ls="", marker='s')

axs1.set_yscale('log')

axs1.grid(True, ls=':')

axs1.set_ylabel(r"$\lambda_n/\lambda_0$")
axs1.set_title("Correlation matrix normalized eigenvalues")

if save:
    fig1.savefig(outfile + '.pdf',transparent=True,bbox_inches='tight')
else:
    mJK.plt.show()


