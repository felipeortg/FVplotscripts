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
    print('Run with parent folder holding this_lattice.py and same folder as dataset with location of all reconfit results')
    print(f'Usage is: {sys.argv[0]} Ecm.xml scat_devel_inter_lvl_file(s)')
    raise Exception('Too few arguments.')

# Import macros and the file with all lattice definitions
# These files are in the same repo as these other libraries (python adds the script folder to the path)
# Also these are in the FVplotscripts git repo, already in the python path
# Therefore no need to append PATH
try:
    import spectrum as spec
    import lattice
except:
    print('Libraries spectrum and/or lattice not found')
    raise Exception('This script is part of a repo and needs the other libraries in the repo to work...')

# Import the name of this lattice_channel
sys.path.append('../')
import this_lattice as tl

# Location of the reconfit dataset to have all levels (used or not) in the plots
sys.path.append('./')
import dataset


# In[2]:
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
plt.rc('font', family='serif')
plt.rc('font', size=14)
plt.rc('axes.formatter', useoffset = False)


# Get the levels (Vol, Ecm, Ecm tot_err, number) in each irrep
# Get the fitted levels in a dictionary
# Key: V_irrep
# Value: dictionary: key lvl_num, value [Ecm, Ecm_err]
levels_in_fit = spec.read_Ecm_xml(sys.argv[1])

# levels_in_fit = spec.read_Ecm_ini(dataset.Ecm_ini)


# In[3]:
# Get the output files from scattering_devel with the interacting spectrum from an amplitude fit
scat_irreps_data = sys.argv[2:]

spectrum_pred = dict()

# scat_irreps_data = spec.reversemom_inirreps(irreps_data)

for scat_devel_file in scat_irreps_data:
    
#     scat_devel_file = header + "output_d" + irs + ".spectrum"

    irrep_int = spec.get_irrep_scat_devel_output(scat_devel_file)
    
    spectrum_pred[irrep_int] = spec.read_interacting_spectrum_scatdevel(scat_devel_file)


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

try:
    header = dataset.dataset
except:
    header ='../../spectrum/data_useme/'

irreps_data = list(spectrum_pred)

full_spectrum_reconfit = dict()

# irreps_int = []

print("Full reconfit out data from: ", header)
print("Reading files:")

for irs in irreps_data:    
    vols = set([lvl[0] for lvl in spectrum_pred[irs]])

    for nn, vol in enumerate(vols):
        inter_file = f"{header}V{vol}/{irs}/reconfit_out.txt"

        print(inter_file)

        irrept, spectrumt, t0strt = spec.read_reconfit_file(inter_file)
        
        energyrange = get_plot_thr_range(irrept[4:])[2]

        spectrum = spec.clean_calc_spectrum(spectrumt, irrept, energyrange, chi, LatticeLength, errorbarwidth)

        full_spectrum_reconfit[f"V{vol}_{irrept}"] = spectrum

    # print(inter_file, spectrumt, full_spectrum_reconfit)
# In[11]:

print('-----\nBegin plots: ', list(spectrum_pred))
for nn, irre in enumerate(spectrum_pred):
    print("Plot:", irre)
    # print(spectrum_pred[irre])
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
        plt.axhline(y=masses[m1] + masses[m2], linestyle = '--', color=colors[particles], lw=2, label= r'$' + names[m1] + names[m2] + '$')

    for eth in extrathresholds:
        axis.axhline(y= eth[0], linestyle = eth[1], color=eth[2], label=eth[3])
        
        
    # Plot lattice data
    # print(levels_in_fit[irre])

    vols = set([lvl[0] for lvl in spectrum_pred[irre]])

    for vol in vols:
        vol_irre = f"V{vol}_{irre}"
        labeled = False

        try:
            for nn,level in enumerate(full_spectrum_reconfit[vol_irre]):
                try:
                    # levels_in_fit contains dictionary with lvl nums, so check if its in it
                    if nn in levels_in_fit[vol_irre]:
                        # use energy an error from the Ecm.xml if in the fit
                        Ecm, Ecm_err = levels_in_fit[vol_irre][nn]
                        levcolor = 'k'
                    else:
                        # use the simple Ecm and error calculation when lvl not used in fit
                        Ecm, Ecm_err = level[1], level[2]
                        levcolor = 'grey'
                # if an irrep is not in the levels_in_fit it is because it was not used in the fit
                except KeyError:
                    Ecm, Ecm_err = level[1], level[2]
                    levcolor = 'grey'


                # add 'data' label when first black data point is added, or if none present make label grey
                if (not labeled and levcolor == 'k') or (not labeled and nn-1 == len(full_spectrum_reconfit[vol_irre])):        
                    plt.errorbar(level[0], Ecm, yerr=Ecm_err, 
                             marker='_', mec=levcolor, fillstyle='none', ecolor=levcolor,
                             elinewidth=1,capsize=10*errorbarwidth/0.8, label = 'data', zorder=4)
                    labeled = True

                # do not repeat label for rest of points
                else:
                    plt.errorbar(level[0], Ecm, yerr=Ecm_err, 
                             marker='_', mec=levcolor, fillstyle='none', ecolor=levcolor,
                             elinewidth=1,capsize=10*errorbarwidth/0.8,zorder=4)

        # no reconfit for this vol_irrep prediction
        except KeyError:
            print("No reconfit data for", vol_irre)
            print("!!! Assuming no fit in", vol_irre)

    # Plot fit result
    fit_result_irre = np.array(spectrum_pred[irre])
    plt.errorbar(fit_result_irre[:,0], fit_result_irre[:,1], yerr=fit_result_irre[:,2], label = 'fit', ls='',
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
        
