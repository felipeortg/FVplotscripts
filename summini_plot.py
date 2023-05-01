#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2023-04-26 16:20:07
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : num
# Use the data from the SUMMARY_GNU_DATA to create plot
# It requires a working installation of awk, sed and grep

import matplotlib.pyplot as plt
import pandas as pd
import io
import subprocess as shell
import sys
from cycler import cycler
# sys.path.append("${HOME}/FVplotscripts")
import iminuitwJK as mJK

if len(sys.argv) < 2:
    raise Exception("Not enough arguments\n"+sys.argv[0]+" SUMMARY_GNU_DATA [filename_to_save] [title:will not save] ")

colors_def = [
    [2/255,103/255,253/255,1],
    [254/255,30/255,4/255,1]
    ]  
default_cycler = cycler(color=colors_def) 

#if len(sys.argv) > 2 and sys.argv[2] == 's':
if len(sys.argv) == 3:
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    plt.rc('font', family='serif')

plt.rc('font', size=11)
plt.rc('axes.formatter', useoffset = False)
plt.rc('lines', linewidth=1)
plt.rc('axes', prop_cycle=default_cycler)

summinifile = str(sys.argv[1])

GREP = "grep"
AWK = "awk"
SED = "sed"

def get_df(index, names):
    sedcomm = SED + ' /#/d ' + summinifile
    awkcomm =AWK + ' -v RS="" '+"'{if (NR == "+'"{0}"'.format(index + 1)+") print $0 }' "
    out = shell.run(sedcomm + " | " + awkcomm, 
        shell=True, capture_output=True).stdout
    
    
    file = io.StringIO(out.decode('UTF-8'))

    return pd.read_csv(file, delim_whitespace=True, names=names)
    


data =[]

# get active
data.append(get_df(0, names=['t','y','e']))

# get inactive
data.append(get_df(1, names=['t','y','e']))

# get function
line = get_df(2, names=['t','y','ym','yp'])




plt.subplots(1,1,figsize=(4.5,3))
ax=plt.gca()

line.plot(x='t',y='y',color='C0',ax=ax, lw=0.7, legend=False)
line.plot(x='t',y='yp',color='C0',ax=ax, lw=0.7, legend=False)
line.plot(x='t',y='ym',color='C0',ax=ax, lw=0.7, legend=False)




data[0].plot(x='t',y='y',yerr='e',color='C0',ax=ax,
    marker='s',ls='',fillstyle='none',capsize=3, ms=3, legend=False)

data[1].plot(x='t',y='y',yerr='e',color='C1',ax=ax,
        marker='s',ls='',fillstyle='none',capsize=3, ms=3, legend=False)


plt.xlabel(r'$t/a_t$')

 
if len(sys.argv) > 3:   
    plt.title(str(sys.argv[3]))

# grepx = shell.run([GREP, "x ", summinifile], capture_output=True).stdout.split()
# xlims = [float(x) for x in grepx[1:]]

# grepy = shell.run([GREP, "y ", summinifile], capture_output=True).stdout.split()
# ylims = [float(y) for y in grepy[1:]]

# plt.xlim(xlims)
# plt.ylim(ylims)

# if mode == 'Jackfitter':
#     label = shell.run([GREP, "gx", summinifile], capture_output=True).stdout.decode('UTF-8')


#     info0 = label.split(sep='"')
#     #print(info0)

#     info1 = info0[1].split(sep='\gx\sp2\ep/N\sbdof\eb')
#     info2 = info1[1].split(';')
#     chiinfo = info2[0]
#     if len(info2) > 1:
#         rest = info2[1]
#         restclean = rest.replace('\\+-', '\\pm')
#     else:
#         restclean = "."

#     fit_info = [
#         r'$\chi^2/\mathrm{dof}' + chiinfo + '$',
#         r'$' + restclean +'$']

# elif mode == 'eigen':

#     linenum = int(shell.run("grep -n '#cs 1' {0}".format(summinifile) + " | sed 's/:/ /g' | awk '{print $1}'", shell=True, capture_output=True).stdout)
#     fit_info = []
#     label = []
#     chiinfo = ""
#     with open(summinifile,'r') as f:
#         for nn, line in enumerate(f):
#             if nn >= linenum:
#                 linewonl = line.rstrip('\n')

#                 if chiinfo == "":
#                     chiinfo = linewonl.split(sep='\gx\sp2\ep/N\sbdof\eb')[1][:-2]
#                     fit_info.append(r'$\chi^2/\mathrm{dof}' + chiinfo + '$',)
#                 else:
#                     # fit_info.append(r'$'+ linewonl.split('"')[1].replace('\\+-', '\\pm') + r'$')
#                     p, ve = linewonl.split('"')[1].split('=')
#                     v, e = ve.split('\\+-')
#                     v=float(v)
#                     e=float(e)
#                     fit_info.append(mJK.add_fit_info_ve(p,v,e))



plt.grid(ls=':')
# plt.legend('',title="\n".join(fit_info), loc='best')

#plt.subplots_adjust(left=0.15)
plt.subplots_adjust(right=0.97)
plt.subplots_adjust(bottom=0.16)

if len(sys.argv) == 3:
# if len(sys.argv) > 2:    
    plt.savefig(str(sys.argv[2])+'.pdf',transparent=True, bbox_inches='tight')
else:
    plt.show()

