#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2022-02-10 15:42:30
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 1.0
# Use the data from the ax file to plot the pc
# It requires a working installation of gnuplot and grep

import matplotlib.pyplot as plt
import pandas as pd
import subprocess as shell
import sys
from cycler import cycler

if len(sys.argv) < 3:
    raise Exception("Not enough arguments\n"+sys.argv[0]+" axfile state# [s:ave_w_tex]")

colors_def = [
    [2/255,103/255,253/255,1],
    [254/255,30/255,4/255,1]
    ]  
default_cycler = cycler(color=colors_def) 

if len(sys.argv) > 3 and sys.argv[3] == 's':
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    plt.rc('font', family='serif')

plt.rc('font', size=11)
plt.rc('axes.formatter', useoffset = False)
plt.rc('lines', linewidth=1)
plt.rc('axes', prop_cycle=default_cycler)

axfile = str(sys.argv[1])

GNUPLOT = "gnuplot"
GREP = "grep"

for ii in range(11):
    shell.run([GNUPLOT, "-e", "set table 'out{0}.dat'; plot '{1}' i {0} with table".format(ii,axfile)])

lines = []
for i in range(9):
    lines.append(pd.read_csv('out{0}.dat'.format(i), delim_whitespace=True, names=['t','y']))


data =[]
data.append(pd.read_csv('out9.dat', delim_whitespace=True, names=['t','y','e']))
data.append(pd.read_csv('out10.dat', delim_whitespace=True, names=['t','y','e']))

shell.run("rm out*.dat", shell=True)

plt.subplots(1,1,figsize=(4.5,3))
ax=plt.gca()

for ii in range(3):
    lines[ii].plot(x='t',y='y',color='C1',ax=ax, lw=0.7)

for ii in range(6,9):
    lines[ii].plot(x='t',y='y',color='C1',ax=ax, lw=0.7)

for ii in range(3,6):
    lines[ii].plot(x='t',y='y',color='C0',ax=ax, lw=0.7)


data[0].plot(x='t',y='y',yerr='e',color='C0',ax=ax,
    marker='s',ls='',fillstyle='none',capsize=3, ms=3)
data[1].plot(x='t',y='y',yerr='e',color='C1',ax=ax,
    marker='s',ls='',fillstyle='none',capsize=3, ms=3)


plt.xlabel(r'$t/a_t$')
plt.title('state ' + str(sys.argv[2]))

grepx = shell.run([GREP, "x ", axfile], capture_output=True).stdout.split()
xlims = [float(x) for x in grepx[1:]]

grepy = shell.run([GREP, "y ", axfile], capture_output=True).stdout.split()
ylims = [float(y) for y in grepy[1:]]

plt.xlim(xlims)
plt.ylim(ylims)

label = shell.run([GREP, "gx", axfile], capture_output=True).stdout.decode('UTF-8')

jj, chi, t0, m = label.split(sep='=')

chil = chi[:-6]
ml = m[:-2].split('\\+-')
t0l = t0[:-1]

fit_info = [
    r'$\chi^2/\mathrm{dof}=' + chil + r'$',
    r'$m= ' + ml[0] + "\\pm" + ml[1] + r'$',
    r'$t_0=' +t0l + r'$'
    ]


plt.legend('',title="\n".join(fit_info), loc='best')

plt.subplots_adjust(left=0.1)
plt.subplots_adjust(right=0.97)
plt.subplots_adjust(bottom=0.16)

if len(sys.argv) > 3 and sys.argv[3] == 's':
    plt.savefig('singleplot.pdf',transparent=True)
else:
    plt.show()
