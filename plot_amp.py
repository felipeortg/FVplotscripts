#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2021-09-08 13:49:21
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 1
# Plot the comparison between an amplitude (via Luescher) and the lattice levels

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 3 :
    print('Run in folder with this_lattice.py')
    print('Usage is: {0} amplitude_plot_label scat_devel_inter_lvl_file(s)'.format(sys.argv[0]))
    raise Exception('Too few arguments.')

# Import macros and the file with all lattice definitions
# These are in the same git repo, and therefore no need to append PATH
try:
    import spectrum as spec
    import lattice
except:
    print('Libraries spectrum and/or lattice not found')
    raise Exception('This script is part of a repo and needs the other libraries in the repo to work...')

# Import the name of this lattice_channel
sys.path.append('./')
import this_lattice as tl



# In[2]:


plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
plt.rc('font', family='serif')
plt.rc('font', size=14)
plt.rc('axes.formatter', useoffset = False)


amplitude = sys.argv[1]

# In[3]:
# Get the output files from scattering_devel with the interacting spectrum from an amplitude fit
scat_irreps_data = sys.argv[2:]

spectrum_pred = dict()

# scat_irreps_data = spec.reversemom_inirreps(irreps_data)

for scat_devel_file in scat_irreps_data:
    
#     scat_devel_file = header + "output_d" + irs + ".spectrum"

    irrep_int = spec.get_irrep_scat_devel_output(scat_devel_file)
    
    spectrum_pred[irrep_int] = np.array(spec.read_interacting_spectrum_scatdevel(scat_devel_file))


# In[4]:
# Get the values of the masses/names/Lattice properties/thresholds/plot options

lattice.init_lattice_props(tl.lattice)

masses = lattice.masses
names = lattice.names

LatticeLength = lattice.LatticeLength
if type(LatticeLength) != int:
    raise Exception('Not implemented yet for more than one volume')
chi = lattice.chi

Lrange = lattice.Lrange
Lspace = np.linspace(Lrange[0], Lrange[1], num=100)

#Styles in the plot
colors = lattice.colors
errorbarwidth = lattice.errorbarwidth


# In[5]:
# Get the values of the free levels

try:
    mode = dataset.mode
except:
    print("Default mode: scat_devel")
    mode = 'scat_devel'
    

if mode == 'scat_devel':
    free_file = lattice.free_file_sca

    two_meson_dict = spec.read_free_spectrum_scatdevel(free_file)

elif mode == 'redstar':
    free_file = lattice.free_file_red

    two_meson_dict = spec.read_redstar_file(free_file)

else:
    raise Exception("Unknown mode: available modes are scat_devel and redstar")




# In[7]:


def get_plot_thr_range(irreirrep):

    if irreirrep in ['T1m', 'A1']:
        tt = lattice.thresholdsm
        et = lattice.extrathresholdsm
        er = lattice.energyrangem
    elif irreirrep in ['T1p', 'A2']:
        tt = lattice.thresholdsp
        et = lattice.extrathresholdsp
        er = lattice.energyrangep
    else:
        tt = lattice.thresholdsmp
        et = lattice.extrathresholdsmp
        er = lattice.energyrangemp
    
    return tt, et, er


# In[9]:

irreps_data = list(spectrum_pred)

spectrum_int = dict()



# In[11]:

print('Begin plots: ', list(spectrum_pred))
for nn, irre in enumerate(spectrum_pred):
    print(irre)
    print(spectrum_pred[irre])
    irreP = spec.label2vec(irre[:3])
    irreirrep = irre[4:]
    
    thresholds,extrathresholds,energyrange  = get_plot_thr_range(irreirrep)

    plt.figure(nn,figsize=(4,6))
    axis = plt.gca()

    orderedmesonlist = spec.order_two_meson_list_at_xiL(two_meson_dict[irre], irreP, LatticeLength*chi, masses)

    # Plot the non-int levels  
    for jj, kdict in enumerate(orderedmesonlist):

        massp1 = masses[kdict['part1']]
        massp2 = masses[kdict['part2']]

        kp1 = spec.label2vec(kdict['k1'])
        kp2 = spec.label2vec(kdict['k2'])

        namep1 = names[kdict['part1']]
        namep2 = names[kdict['part2']]
        namepair = kdict['part1'] + 'xx' + kdict['part2']

        energyatupperL = spec.free_part_en(massp1, massp2, Lspace[-1]*chi,
                                       kp1, kp2, irreP)

        if energyatupperL > energyrange[-1]:
            continue                                         

        axis.plot(Lspace, spec.free_part_en(massp1, massp2, Lspace*chi,
                                       kp1, kp2, irreP),
        color=colors[namepair], label = r'$' + namep1 + '_{' + kdict['k1'] + '}'
                  + namep2 + '_{' + kdict['k2'] +'}$', zorder=4)

    # Plot thresholds
    for mm, particles in enumerate(thresholds):
        m1, m2 = particles.split('xx')
        plt.axhline(y=masses[m1] + masses[m2], linestyle = '--', color=colors[particles], label= r'$' + names[m1] + names[m2] + '$')

    for eth in extrathresholds:
        axis.axhline(y= eth[0], linestyle = eth[1], color=eth[2], label=eth[3])
        
    # Plot amplitude result
    plt.errorbar(spectrum_pred[irre][:,0], spectrum_pred[irre][:,1], yerr=spectrum_pred[irre][:,2], label = amplitude, ls='',
                     marker='o', mec='orangered', fillstyle='none', ecolor='orangered',elinewidth=1,capsize=5*errorbarwidth/0.8,zorder=3)  
        
        
    axis.legend(loc="center right", bbox_to_anchor=(2, 0.5))

    axis.set_title(r'$P = ' + str(irreP) + r'$, ' + irreirrep)
    axis.set_ylim(energyrange)
    axis.set_xlim(Lrange)
    plt.xticks(np.linspace(Lrange[0], Lrange[1], 3),
                    [str(Lrange[0]),str(int((Lrange[0]+Lrange[1])/2)),str(Lrange[1])])
    axis.set_xlabel(r'$L/a_s$')
    axis.set_ylabel(r'$a_t\,E_\text{cm}$')

    plt.subplots_adjust(left=0.2)
    plt.subplots_adjust(right=0.6)
    # plt.subplots_adjust(bottom=0.13)

    plt.savefig("plots/" + irreirrep +'_'+ spec.vec2label(irreP) + '.pdf', transparent=True)
        
