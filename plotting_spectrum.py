#!/usr/bin/env python3
# coding: utf-8
# @Date    : 2021-06-21 12:04:18
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)


import numpy as np
import matplotlib.pyplot as plt
import sys

# Import macros and the file with all lattice definitions
sys.path.append('/Users/Felipe/Google Drive/bin/')
import spectrum as spec
import lattice

# Import the name of this lattice_channel
sys.path.append('../../../')
import this_lattice as tl

# Get the values of the masses/names/Lattice properties/thresholds/plot options
lattice.init_lattice_props(tl.lattice)

masses = lattice.masses
names = lattice.names

LatticeLength = lattice.LatticeLength
chi = lattice.chi

Lrange = lattice.Lrange
Lspace = np.linspace(Lrange[0], Lrange[1], num=100)

colors = lattice.colors
errorbarwidth = lattice.errorbarwidth

# In[2]:


plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
plt.rc('font', family='serif')
plt.rc('font', size=14)
plt.rc('axes.formatter', useoffset = False)



# In[4]:
# Get the values of the free levels

if len(sys.argv) > 2:
    mode = sys.argv[2]
else:
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



# In[5]:

inter_file = sys.argv[1]

irrep, spectrum1, t0str = spec.read_reconfit_file(inter_file)

irreP = spec.label2vec(irrep[:3])
irreirrep = irrep[4:]

# In[6]:

if irreirrep in ['T1m', 'A1']:
    thresholds = lattice.thresholdsm
    extrathresholds = lattice.extrathresholdsm
    energyrange = lattice.energyrangem
elif irreirrep in ['T1p', 'A2']:
    thresholds = lattice.thresholdsp
    extrathresholds = lattice.extrathresholdsp
    energyrange = lattice.energyrangep
else:
    thresholds = lattice.thresholdsmp
    extrathresholds = lattice.extrathresholdsmp
    energyrange = lattice.energyrangemp


# In[7]:

spectrum = spec.clean_calc_spectrum(spectrum1, irrep, energyrange, chi, LatticeLength, errorbarwidth)
    


plt.figure(1,figsize=(4,6))
axis = plt.gca()

orderedmesonlist = spec.order_two_meson_list_at_xiL(two_meson_dict[irrep], irreP, LatticeLength*chi, masses)


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


# axis.axvline(x= LatticeLength, linestyle = '--', color='gray')

for mm, particles in enumerate(thresholds):
    m1, m2 = particles.split('xx')
    plt.axhline(y=masses[m1] + masses[m2], linestyle = '--', color=colors[particles], label= r'$' + names[m1] + names[m2] + '$')

for eth in extrathresholds:
    axis.axhline(y= eth[0], linestyle = eth[1], color=eth[2], label=eth[3])


for level in spectrum:
    plt.errorbar(level[0], level[1], yerr=level[2], 
                 marker='_', mec='k', fillstyle='none', ecolor='k',elinewidth=1,capsize=10,zorder=4)
    plt.plot()


axis.legend(loc="center right", bbox_to_anchor=(2, 0.5))

axis.set_title(r'$P = ' + str(irreP) + r'$, ' + irreirrep + '\n' + t0str)

axis.set_ylim(energyrange)

axis.set_xlim(Lrange)
plt.xticks(np.linspace(Lrange[0], Lrange[1], 3),
                [str(Lrange[0]),str(int((Lrange[0]+Lrange[1])/2)),str(Lrange[1])])
axis.set_xlabel(r'$L/a_s$')
axis.set_ylabel(r'$a_t\,E_\text{cm}$')

plt.subplots_adjust(left=0.2)
plt.subplots_adjust(right=0.6)
plt.subplots_adjust(bottom=0.13)

plt.savefig(irreirrep +'_'+ spec.vec2label(irreP) + '.pdf', transparent=True)



