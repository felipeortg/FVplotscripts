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

datalen = int(shell.run("grep '#e' " + axfile + " | wc -l", shell=True, capture_output=True).stdout)

if datalen == 12:
    print("axfile with data and unused data")
elif datalen == 11:
    print("axfile only with data, no unused data")
else:
    raise Exception('There is not 9 lines , and 1 or 2 data used(unused) sets, but instead {0} sets of points'.format(datalen-1))

for ii in range(datalen - 1):
    shell.run([GNUPLOT, "-e", "set table 'out{0}.dat'; plot '{1}' i {0} with table".format(ii,axfile)])


#out=shell.run('head -10 out*.dat', shell=True, capture_output=True)
#print('a\n',out.stdout.decode('UTF-8'))

lines = []
for i in range(9):
    lines.append(pd.read_csv('out{0}.dat'.format(i), delim_whitespace=True, names=['t','y']))


data =[]
data.append(pd.read_csv('out9.dat', delim_whitespace=True, names=['t','y','e']))
if datalen > 11:
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
if datalen > 11:
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


info0 = label.split(sep='"')
#print(info0)

info1 = info0[1].split(sep='\gx\sp2\ep/N\sbdof\eb')
info2 = info1[1].split(';')
chiinfo = info2[0]
if len(info2) > 1:
    rest = info2[1]
    restclean = rest.replace('\\+-', '\\pm')
else:
    restclean = "."

fit_info = [
    r'$\chi^2/\mathrm{dof}' + chiinfo + '$',
    r'$' + restclean +'$']


plt.legend('',title="\n".join(fit_info), loc='best')

#plt.subplots_adjust(left=0.15)
plt.subplots_adjust(right=0.97)
plt.subplots_adjust(bottom=0.16)

if len(sys.argv) > 3 and sys.argv[3] == 's':
    plt.savefig('prin_corr_t0'+str(t0l)+'_state'+str(sys.argv[2])+'.pdf',transparent=True)
else:
    plt.show()
