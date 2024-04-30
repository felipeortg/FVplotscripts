#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2022-11-08 11:52:01
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 1.0
# Create correlation plot of all files in a folder, ordering might be given by second set of arguments

import iminuitwJK as mJK
import subprocess as sp
import pandas as pd
import spectrum as spec
import sys

if( len(sys.argv) < 3):
    print(f"Usage is: {sys.argv[0]} list_file folder [Ecm.ini mass_folder]")
    exit(1)

datalist_file = str(sys.argv[1])
datafolder = str(sys.argv[2])

if( len(sys.argv) > 3):
    enelist_file = str(sys.argv[3])
    enefolder = str(sys.argv[4])
else:
    # if not supplied, then order with its own folder (case of energies)
    enelist_file = datalist_file
    enefolder = datafolder   


mJK.plt.rc('text', usetex=True)
mJK.plt.rc('text.latex', preamble=r'\usepackage{amsmath,amssymb}')
mJK.plt.rc('font', family='serif')
mJK.plt.rc('font', size=10)
mJK.plt.rc('axes.formatter', useoffset = False)

ecm_list = spec.read_ini_file(enelist_file)
data_list = spec.read_ini_file(datalist_file)
print(ecm_list)

ecm = []
data = []
irrep_list = []
lvl_num = []
vol=24



for irrep in ecm_list:
    for level in ecm_list[irrep]:
        label = spec.label_state(level)

        print(irrep)
        jackfile = level['file_name']
        tmp = mJK.get_data(f"{enefolder}/" + jackfile)
        ecm.append(tmp[3])
        irrep_list.append(irrep)
        lvl_num.append(level['lvl_num'])

        data_file = None
        for data_level in data_list[irrep]:
            if spec.label_state(data_level) == label:
                data_file = data_level["file_name"]

        tmp = mJK.get_data(f"{datafolder}/" + data_file)
        data.append(tmp[3])

data_df = pd.DataFrame(irrep_list,columns=["irrep"])
data_df["irrep"] = data_df["irrep"].apply(lambda x: "[" + x[:3] + "] $" + x[4] + "_" + x[5] + ("^-" if x[-1]=="m" else "") + "$")
data_df["lvl"] = lvl_num


ecm_me = mJK.get_mean_error_array(ecm)
data_me = mJK.get_mean_error_array(data)

data_df["ecm"] = ecm_me[:,0]
data_df["data"] = data_me[:,0]
data_df["data_err"] = data_me[:,1]
data_df["data_perr"] = mJK.np.abs(data_df["data_err"]/data_df["data"])*100
data_df["x_labels"] = data_df["irrep"] + "\#" + data_df["lvl"].apply(str)

label_df = data_df.sort_values('ecm')
dataord = list(label_df.index)
dataordirrep = list(label_df["irrep"] + "\#" + label_df["lvl"].apply(str))


temp = []
for val in dataord:
    temp.append(data[val])

data_energy_ord = mJK.combine_ensemble_list(temp)

corry_m = mJK.cormatense(data_energy_ord, mJK.meanense(data_energy_ord))

#plot correlation
mJK.plt.rc('font', size=6)

fig,ax = mJK.plt.subplots(1,1,num=1)

mJK.correlation_plot(ax, dataordirrep, corry_m, label_all=True)

ax.set_xticklabels([str(dataordirrep[i]) for i in ax.get_xticks()], rotation=45, ha='left')

mJK.add_labels_matrix(ax, corry_m, hide_diag=True, size=6)

mJK.plt.savefig(f"{datafolder}/correl.pdf", transparent=True, bbox_inches='tight')

mJK.plt.rc('font', size=10)

#plot svds

fig,ax = mJK.plt.subplots(1,1,num=2,figsize=(4,3))
eig_val_data = mJK.np.sort(mJK.np.linalg.eigvals(corry_m))
eig_val_data /= mJK.np.max(eig_val_data)

mJK.plt.plot(eig_val_data[::-1],ls='',marker='s',mfc='w')

ax.semilogy()
ax.grid(ls=':',lw=0.6)
ax.set_title(f'{datafolder} correlation matrix eigenvalues')

mJK.plt.savefig(f"{datafolder}/corr_svd.pdf", transparent=True, bbox_inches='tight')

#plot percentage error
fig,ax = mJK.plt.subplots(1,1,num=3,figsize=(8,3))
label_df.plot(x="x_labels", y="data_perr", ax=ax, ls='', marker='o',color='k', mfc='w', label="\% error")


# ax.set_ylim([0,1.4])
ax.set_xticks(mJK.np.arange(len(dataordirrep)))
ax.set_xticklabels(label_df["x_labels"], rotation=90, ha='center')
ax.set_xlabel("")

ax.grid(ls=":", lw=0.7)

mJK.plt.savefig(f"{datafolder}/perc_err.pdf", transparent=True, bbox_inches='tight')




