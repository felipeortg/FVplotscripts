#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2022-02-10 15:42:30
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 1.0
# Use the data from the ax file to plot the pc
# It requires a working installation of awk, sed and grep

import matplotlib.pyplot as plt
import pandas as pd
import io
import subprocess as shell
import sys
from cycler import cycler
sys.path.append("${HOME}/FVplotscripts")
import iminuitwJK as mJK

if len(sys.argv) < 2:
    raise Exception("Not enough arguments\n"+sys.argv[0]+" axfile state# [s:ave_w_tex]")

colors_def = [
    [2/255,103/255,253/255,1],
    [254/255,30/255,4/255,1]
    ]  
default_cycler = cycler(color=colors_def) 

#if len(sys.argv) > 2 and sys.argv[2] == 's':
if len(sys.argv) > 2:
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    plt.rc('font', family='serif')

plt.rc('font', size=11)
plt.rc('axes.formatter', useoffset = False)
plt.rc('lines', linewidth=1)
plt.rc('axes', prop_cycle=default_cycler)

axfile = str(sys.argv[1])

GREP = "grep"
AWK = "awk"
SED = "sed"

def get_df(index, names):
    sedcomm = SED + ' /#/d ' + axfile
    awkcomm =AWK + ' -v RS="" '+"'{if (NR == "+'"{0}"'.format(index + 1)+") print $0 }' "
    out = shell.run(sedcomm + " | " + awkcomm, 
        shell=True, capture_output=True).stdout
    
    
    file = io.StringIO(out.decode('UTF-8'))
    
    return pd.read_csv(file, delim_whitespace=True, names=names)
    

datalen = int(shell.run(GREP + " '#e' " + axfile + " | wc -l", shell=True, capture_output=True).stdout)

mode='Jackfitter'
lineslen = 9

if datalen == 12:
    print("axfile with data and unused data")
elif datalen == 11:
    print("axfile only with data, no unused data")
elif datalen == 6:
    mode='eigen'
    lineslen = 3
    print("axfile from eigen")
else:
    raise Exception('Jackfitter files have 9 lines and 1 or 2 data used(unused) sets,\
     eigen has 3 lines and 2 data sets, but {1} has instead {0} sets of points'.format(datalen-1, axfile))
    
lines = []
data =[]

for ii in range(lineslen):
    lines.append(get_df(ii, names=['t','y']))
    

data.append(get_df(lineslen, names=['t','y','e']))
if datalen > 11 or mode=='eigen':
    data.append(get_df(lineslen+1, names=['t','y','e']))


plt.subplots(1,1,figsize=(4.5,3))
ax=plt.gca()

color3first='C1'
if mode == 'eigen':
    color3first='C0'

for ii in range(3):
    lines[ii].plot(x='t',y='y',color=color3first,ax=ax, lw=0.7)

if mode == 'Jackfitter':
    for ii in range(6,9):
        lines[ii].plot(x='t',y='y',color='C1',ax=ax, lw=0.7)

    for ii in range(3,6):
        lines[ii].plot(x='t',y='y',color='C0',ax=ax, lw=0.7)


data[0].plot(x='t',y='y',yerr='e',color='C0',ax=ax,
    marker='s',ls='',fillstyle='none',capsize=3, ms=3)
if datalen > 11 or mode=='eigen':
    data[1].plot(x='t',y='y',yerr='e',color='C1',ax=ax,
        marker='s',ls='',fillstyle='none',capsize=3, ms=3)


plt.xlabel(r'$t/a_t$')

 
if len(sys.argv) > 2 and sys.argv[2] != 's':   
    plt.title(str(sys.argv[2]))

grepx = shell.run([GREP, "x ", axfile], capture_output=True).stdout.split()
xlims = [float(x) for x in grepx[1:]]

grepy = shell.run([GREP, "y ", axfile], capture_output=True).stdout.split()
ylims = [float(y) for y in grepy[1:]]

plt.xlim(xlims)
plt.ylim(ylims)

if mode == 'Jackfitter':
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

elif mode == 'eigen':

    linenum = int(shell.run("grep -n '#cs 1' {0}".format(axfile) + " | sed 's/:/ /g' | awk '{print $1}'", shell=True, capture_output=True).stdout)
    fit_info = []
    label = []
    chiinfo = ""
    with open(axfile,'r') as f:
        for nn, line in enumerate(f):
            if nn >= linenum:
                linewonl = line.rstrip('\n')

                if chiinfo == "":
                    chiinfo = linewonl.split(sep='\gx\sp2\ep/N\sbdof\eb')[1][:-2]
                    fit_info.append(r'$\chi^2/\mathrm{dof}' + chiinfo + '$',)
                else:
                    # fit_info.append(r'$'+ linewonl.split('"')[1].replace('\\+-', '\\pm') + r'$')
                    p, ve = linewonl.split('"')[1].split('=')
                    v, e = ve.split('\\+-')
                    v=float(v)
                    e=float(e)
                    fit_info.append(mJK.add_fit_info_ve(p,v,e))



plt.grid(ls=':')
plt.legend('',title="\n".join(fit_info), loc='best')

#plt.subplots_adjust(left=0.15)
plt.subplots_adjust(right=0.97)
plt.subplots_adjust(bottom=0.16)

if len(sys.argv) > 2 and sys.argv[2] == 's':
# if len(sys.argv) > 2:    
    plt.savefig(str(sys.argv[2])+'.pdf',transparent=True, bbox_inches='tight')
else:
    plt.show(block=True)

