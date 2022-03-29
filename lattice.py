#!/usr/bin/env python3
# coding: utf-8
# @Date    : 2021-09-2 12:04:18
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)

drive_loc = '/Users/felipeortg/Google Drive/My Drive/'

LatticeLength = None
chi = None

masses = None
names = None


colors = None


thresholdsm = None
thresholdsp = None
thresholdsmp = None

ls_types = None

extrathresholdsm = None
extrathresholdsp = None
extrathresholdsmp = None

linestyles_types = ['solid', 'dotted', 'dashdot', (0, (3, 1, 1, 1, 1, 1))]

errorbarwidth = None

Lrange = None

energyrangem = None
energyrangep = None
energyrangemp = None

free_file_red = None

free_file_sca = None


def init_lattice_props(lattname):
    # particles
    global masses,names,thresholdsm,thresholdsp,thresholdsmp,extrathresholdsm,extrathresholdsp,extrathresholdsmp
    # Lattice
    global LatticeLength,chi,Lrange,colors,free_file_sca,free_file_red
    # Plots
    global energyrangem, energyrangep, energyrangemp, errorbarwidth

    # Parameters of the L=24 a_856 lattice
    if lattname == '24_a856_IG_1p':


        # TABLE I of 1904.03188
        mpion =  0.04735
        mkaon = 0.08659
        meta = 0.09602
        # Fits in the cluster: omega/fits_rge/000_T1mM.fewer2 t011
        momega = 0.1422
        mphi = 0.1709

        # eta/fits_rge/000_A1mP.fewer t010
        metap = 0.1373


        LatticeLength = 24
        chi = 3.455


        mesons = {
            'pion_proj0' : {'mass' : mpion, 'name' : r'\pi'},
            'Kneg_proj0' : {'mass' : mkaon, 'name' : r'K'},
            'Kbarneg_proj0' : {'mass' : mkaon, 'name' : r'\overline{K}'},
            'omega_proj0' : {'mass' : momega, 'name' : r'\omega'},
            'omega_proj1' : {'mass' : mphi, 'name' : r'\phi'},
            'pi' : {'mass' : mpion, 'name' : r'\pi'},
            'kaon' : {'mass' : mkaon, 'name' : r'K'},
            'kbar' : {'mass' : mkaon, 'name' : r'\overline{K}'},
            'omega' : {'mass' : momega, 'name' : r'\omega'},
            'phi' : {'mass' : mphi, 'name' : r'\phi'}
        }

        masses = dict()
        for key in mesons: masses[key] = mesons[key]['mass']

        names = dict()
        for key in mesons: names[key] = mesons[key]['name']


        colors = {
            'pion_proj0xxpion_proj0' : 'blue',
            'Kneg_proj0xxKbarneg_proj0' : 'red',
            'Kbarneg_proj0xxKneg_proj0' : 'red',
            'pion_proj0xxomega_proj0' : 'green',
            'pion_proj0xxomega_proj1' : 'brown',
            'pixxpi' : 'blue',
            'kaonxxkbar' : 'red',
            'kbarxxkaon' : 'red',
            'pixxomega' : 'green',
            'pixxphi' : 'brown'
        }


        thresholdsm = ['pixxpi','kaonxxkbar','pixxomega','pixxphi']
        thresholdsp = thresholdsm
        thresholdsmp = thresholdsm

        ls_types = ['solid', 'dotted', 'dashdot', (0, (3, 1, 1, 1, 1, 1))]

        extrathresholdsm = [
            [2*mpion + meta, ls_types[1], 'pink', r'$\pi\pi \eta$'],
            [mpion + 2*mkaon, ls_types[2], 'pink', r'$\pi K K$'],
            [4*mpion, ls_types[2], 'green', r'$\pi\pi \pi\pi$' ]
        ]
        extrathresholdsp = extrathresholdsm
        extrathresholdsmp = extrathresholdsm


        linestyles_types = ['solid', 'dotted', 'dashdot', (0, (3, 1, 1, 1, 1, 1))]

        Lrange = [20,28]

        energyrangem = [0.09,0.23]
        energyrangep = [0.17,0.25]
        energyrangemp = [0.09,0.25]


        errorbarwidth = 0.8

        free_file_red = drive_loc + 'Documents/JLab/time_like_form_factor/self notes/free_levels/a856/free_levels_redstar.txt'

        free_file_sca = drive_loc + 'Documents/JLab/time_like_form_factor/scattering_24_a856/free_levels/free_levels_scattering_devel_mod.txt'

    else:
        print(lattname)
        raise Exception('This lattice properties are not in here')
