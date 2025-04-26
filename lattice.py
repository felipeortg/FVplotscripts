#!/usr/bin/env python3
# coding: utf-8
# @Date    : 2021-09-2 12:04:18
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)

drive_loc = '/Users/felortga/Google Drive/My Drive/'

LatticeLength = None
chi = None
xi = None
nus = None
at = None

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
    errorflag = False
    # particles
    global masses,names,thresholdsm,thresholdsp,thresholdsmp,extrathresholdsm,extrathresholdsp,extrathresholdsmp
    # Lattice
    global LatticeLength,chi,xi,nus,at,Lrange,colors,free_file_sca,free_file_red
    # Plots
    global energyrangem, energyrangep, energyrangemp, errorbarwidth

    # Parameters of the L=24 a_856 lattice
    if lattname[:7] == '24_a856':


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
        nus = 4.3 / 3.4
        at = 0.2793 / 1.67245

        mesons = {
            'pion_proj0' : {'mass' : mpion, 'name' : r'\pi{}'},
            'Kneg_proj0' : {'mass' : mkaon, 'name' : r'K'},
            'Kbarneg_proj0' : {'mass' : mkaon, 'name' : r'\overline{K}'},
            'omega_proj0' : {'mass' : momega, 'name' : r'\omega{}'},
            'omega_proj1' : {'mass' : mphi, 'name' : r'\phi{}'},
            'pi' : {'mass' : mpion, 'name' : r'\pi{}'},
            'kaon' : {'mass' : mkaon, 'name' : r'K'},
            'kbar' : {'mass' : mkaon, 'name' : r'\overline{K}'},
            'omega' : {'mass' : momega, 'name' : r'\omega{}'},
            'phi' : {'mass' : mphi, 'name' : r'\phi{}'},
            'eta' : {'mass' : meta, 'name' : r'\eta{}'}
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
            'pixxphi' : 'brown',
            'pixxkaon' : 'red',
            'etaxxkaon' : 'green',
        }



        ls_types = ['solid', 'dotted', 'dashdot', (0, (3, 1, 1, 1, 1, 1))]




        if lattname[7:] == '_IG_1p':
            thresholdsm = ['pixxpi','kaonxxkbar','pixxomega','pixxphi']
            thresholdsp = ['pixxomega','pixxphi']
            thresholdsmp = thresholdsm

            extrathresholdsm = [
            [2*mpion + meta, ls_types[1], 'pink', r'$\pi\pi \eta$'],
            [mpion + 2*mkaon, ls_types[2], 'pink', r'$\pi K K$'],
            [4*mpion, ls_types[2], 'green', r'$\pi\pi \pi\pi$' ]
            ]
            extrathresholdsp = extrathresholdsm
            extrathresholdsmp = extrathresholdsm

            energyrangem = [0.09,0.23]
            energyrangep = [0.17,0.25]
            energyrangemp = [0.09,0.25]

            free_file_red = drive_loc + 'Documents/JLab/time_like_form_factor/self notes/free_levels/a856/free_levels_redstar.txt'

            free_file_sca = drive_loc + 'Documents/JLab/time_like_form_factor/scattering_24_a856/free_levels/free_levels_scattering_devel_mod.txt'

        elif lattname[7:] == '_IG_1o2':
            thresholdsm = ['pixxkaon','etaxxkaon']
            thresholdsp = []
            thresholdsmp = thresholdsm

            extrathresholdsm = [
            [2*mpion + mkaon, ls_types[2], 'pink', r'$\pi\pi K$'],
            ]
            extrathresholdsp = []
            extrathresholdsmp = extrathresholdsm

            energyrangem = [0.1,0.22]
            energyrangep = [0.17,0.25]
            energyrangemp = [0.1,0.22]

            free_file_sca = drive_loc + 'Documents/JLab/rhoformfactor/kinematic region per lattice/xcheck/856/free_levels_scattering_devel_mod.txt'
        elif lattname[7:] != "":
            errorflag = True


 
        linestyles_types = ['solid', 'dotted', 'dashdot', (0, (3, 1, 1, 1, 1, 1))]

        Lrange = [20,28]




        errorbarwidth = 0.8


    # Parameters of the a_840 lattice
    elif lattname[3:] == 'a840_IG_1p':


        # TABLE of 2008.06432 and references therein
        mpion =  0.06906 # 1203.6041
        mkaon = 0.09698
        meta = 0.10364
        
        momega = 0.15541
        mphi = 0.17949

        
        LatticeLength = int(lattname[:2])
        if not (LatticeLength in [16,20,24]):
            raise Exception('This volume: {0}, was not added'.format(LatticeLength))
        chi = 3.444 # 1203.6041


        mesons = {
            'pion_proj0' : {'mass' : mpion, 'name' : r'\pi{}'},
            'Kneg_proj0' : {'mass' : mkaon, 'name' : r'K'},
            'Kbarneg_proj0' : {'mass' : mkaon, 'name' : r'\overline{K}'},
            'omega_proj0' : {'mass' : momega, 'name' : r'\omega{}'},
            'omega_proj1' : {'mass' : mphi, 'name' : r'\phi{}'},
            'pi' : {'mass' : mpion, 'name' : r'\pi{}'},
            'kaon' : {'mass' : mkaon, 'name' : r'K'},
            'kbar' : {'mass' : mkaon, 'name' : r'\overline{K}'},
            'omega' : {'mass' : momega, 'name' : r'\omega{}'},
            'phi' : {'mass' : mphi, 'name' : r'\phi{}'},
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


        thresholdsm = ['pixxpi','kaonxxkbar','pixxomega']
        thresholdsp = ['pixxomega']
        thresholdsmp = thresholdsm

        ls_types = ['solid', 'dotted', 'dashdot', (0, (3, 1, 1, 1, 1, 1))]

        extrathresholdsm = []
        extrathresholdsp = extrathresholdsm
        extrathresholdsmp = extrathresholdsm


        linestyles_types = ['solid', 'dotted', 'dashdot', (0, (3, 1, 1, 1, 1, 1))]

        Lrange = [LatticeLength-4,LatticeLength+4]

        energyrangem = [0.125,0.2]
        energyrangep = [0.2,0.3]
        energyrangemp = [0.125,0.2]


        errorbarwidth = 0.8

        free_file_red = None

        free_file_sca = drive_loc + 'Documents/JLab/time_like_form_factor/self notes/free_levels/a840/free_levels_scattering_devel_mod_' + str(LatticeLength) + '.txt'

    elif lattname == '16_20_24_a840_IG_1p':


        # TABLE of 2008.06432 and references therein
        mpion =  0.06906
        mkaon = 0.09698
        meta = 0.10364
        
        momega = 0.15541
        mphi = 0.17949

        
        LatticeLength = 20

        chi = 3.444


        mesons = {
            'pion_proj0' : {'mass' : mpion, 'name' : r'\pi{}'},
            'Kneg_proj0' : {'mass' : mkaon, 'name' : r'K'},
            'Kbarneg_proj0' : {'mass' : mkaon, 'name' : r'\overline{K}'},
            'omega_proj0' : {'mass' : momega, 'name' : r'\omega{}'},
            'omega_proj1' : {'mass' : mphi, 'name' : r'\phi{}'},
            'pi' : {'mass' : mpion, 'name' : r'\pi{}'},
            'kaon' : {'mass' : mkaon, 'name' : r'K'},
            'kbar' : {'mass' : mkaon, 'name' : r'\overline{K}'},
            'omega' : {'mass' : momega, 'name' : r'\omega{}'},
            'phi' : {'mass' : mphi, 'name' : r'\phi{}'},
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


        thresholdsm = ['pixxpi','kaonxxkbar','pixxomega']
        thresholdsp = ['pixxomega']
        thresholdsmp = thresholdsm

        ls_types = ['solid', 'dotted', 'dashdot', (0, (3, 1, 1, 1, 1, 1))]

        extrathresholdsm = []
        extrathresholdsp = extrathresholdsm
        extrathresholdsmp = extrathresholdsm


        linestyles_types = ['solid', 'dotted', 'dashdot', (0, (3, 1, 1, 1, 1, 1))]

        Lrange = [LatticeLength-6,LatticeLength+6]

        energyrangem = [0.125,0.2]
        energyrangep = [0.2,0.3]
        energyrangemp = [0.125,0.2]


        errorbarwidth = 0.8

        free_file_red = None

        free_file_sca = drive_loc + 'Documents/JLab/time_like_form_factor/self notes/free_levels/a840/free_levels_scattering_devel_mod_24.txt'


    elif lattname == '24_a840_IG_2p':

        # TABLE of 2008.06432 and references therein
        mpion =  0.06906 # 1203.6041
        
        LatticeLength = int(lattname[:2])

        chi = 3.444 # 1203.6041
        nus = 4.3 / 3.4

        mesons = {
            'pion_proj0' : {'mass' : mpion, 'name' : r'\pi{}'},
            'pi' : {'mass' : mpion, 'name' : r'\pi{}'}}

        masses = dict()
        for key in mesons: masses[key] = mesons[key]['mass']

        names = dict()
        for key in mesons: names[key] = mesons[key]['name']

        colors = {
            'pion_proj0xxpion_proj0' : 'blue',
            'pixxpi' : 'blue',}

        thresholdsp = ['pixxpi']
        thresholdsm = thresholdsp
        thresholdsmp = thresholdsp

        ls_types = ['solid', 'dotted', 'dashdot', (0, (3, 1, 1, 1, 1, 1))]
        
        extrathresholdsm = [[4*mpion, ls_types[2], 'green', r'$\pi\pi \pi\pi$' ]]
        extrathresholdsp = extrathresholdsm
        extrathresholdsmp = extrathresholdsm
        linestyles_types = ['solid', 'dotted', 'dashdot', (0, (3, 1, 1, 1, 1, 1))]

        Lrange = [LatticeLength-4,LatticeLength+4]

        energyrangem = [0.125,0.35]
        energyrangep = [0.125,0.35]
        energyrangemp = [0.125,0.35]


        errorbarwidth = 0.8

        free_file_red = None
        free_file_sca = drive_loc + 'Documents/tmp_JLab/pipi_I2_scat/free_levels/free_levels_scattering_devel.txt'

    # Parameters of the L=24 a_850 lattice
    elif lattname == '24_a850_IG_1p':


        # TABLE I of 1904.03188
        mpion = 0.05593
        mkaon = 0.09027
        meta = 0.09790
        # Fits in the cluster: omega/fits_rge/000_T1mM.no_2/ t012
        momega = 0.14828
        mphi = 0.17459


        LatticeLength = 24
        chi = 3.456


        mesons = {
            'pion_proj0' : {'mass' : mpion, 'name' : r'\pi{}'},
            'Kneg_proj0' : {'mass' : mkaon, 'name' : r'K'},
            'Kbarneg_proj0' : {'mass' : mkaon, 'name' : r'\overline{K}'},
            'omega_proj0' : {'mass' : momega, 'name' : r'\omega{}'},
            'omega_proj1' : {'mass' : mphi, 'name' : r'\phi{}'},
            'pi' : {'mass' : mpion, 'name' : r'\pi{}'},
            'kaon' : {'mass' : mkaon, 'name' : r'K'},
            'kbar' : {'mass' : mkaon, 'name' : r'\overline{K}'},
            'omega' : {'mass' : momega, 'name' : r'\omega{}'},
            'phi' : {'mass' : mphi, 'name' : r'\phi{}'},
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
        thresholdsp = ['pixxomega','pixxphi']
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

        energyrangem = [0.09,0.2]
        energyrangep = [0.17,0.25]
        energyrangemp = [0.09,0.2]


        errorbarwidth = 0.8

        free_file_sca = drive_loc + 'Documents/JLab/time_like_form_factor/self notes/free_levels/a850/free_levels_scattering_devel_mod_24.txt'

    # Parameters of the a_860 lattice
    elif lattname[-10:] == 'a860_IG_1p':


        # from of https://arxiv.org/pdf/1507.02599.pdf
        mpion =  0.03928
        mkaon = 0.08344
        meta = 0.09299


        try:
            LatticeLength = int(lattname[:-11])
            if not (LatticeLength in [24, 32]):
                raise Exception('This volume: {0}, was not added'.format(LatticeLength))

            Lrange = [LatticeLength-4,LatticeLength+4]
        except:
            if lattname[:-11] == '24_32':
                LatticeLength = 28
            else:
                raise Exception('This volume: {0}, was not added'.format(lattname[:-9]))

            Lrange = [LatticeLength-8,LatticeLength+8]

        chi = 3.453


        mesons = {
            'pion_proj0' : {'mass' : mpion, 'name' : r'\pi{}'},
            'Kneg_proj0' : {'mass' : mkaon, 'name' : r'K'},
            'Kbarneg_proj0' : {'mass' : mkaon, 'name' : r'\overline{K}'},
            'pi' : {'mass' : mpion, 'name' : r'\pi{}'},
            'kaon' : {'mass' : mkaon, 'name' : r'K'},
            'kbar' : {'mass' : mkaon, 'name' : r'\overline{K}'},
        }

        masses = dict()
        for key in mesons: masses[key] = mesons[key]['mass']

        names = dict()
        for key in mesons: names[key] = mesons[key]['name']


        colors = {
            'pion_proj0xxpion_proj0' : 'blue',
            'Kneg_proj0xxKbarneg_proj0' : 'red',
            'Kbarneg_proj0xxKneg_proj0' : 'red',
            'pixxpi' : 'blue',
            'kaonxxkbar' : 'red',
            'kbarxxkaon' : 'red',
        }


        thresholdsm = ['pixxpi','kaonxxkbar']
        thresholdsp = []
        thresholdsmp = thresholdsm

        ls_types = ['solid', 'dotted', 'dashdot', (0, (3, 1, 1, 1, 1, 1))]

        extrathresholdsm = [[4*mpion, ls_types[2], 'green', r'$\pi\pi \pi\pi$' ]]
        extrathresholdsp = extrathresholdsm
        extrathresholdsmp = extrathresholdsm


        linestyles_types = ['solid', 'dotted', 'dashdot', (0, (3, 1, 1, 1, 1, 1))]

        

        energyrangem = [0.07,0.17]
        energyrangep = [0.2,0.3]
        energyrangemp = [0.07,0.17]


        errorbarwidth = 0.8

        free_file_red = None

        free_file_sca = drive_loc + 'Documents/JLab/rhoformfactor/kinematic region per lattice/interacting_levels/860/free_levels_scattering_devel_mod.txt'
    else:
        errorflag = True


    if errorflag:
        print(lattname)
        print("Saved lattice info are:")
        print("""
    24_a856
    24_a856_IG_1p
    24_a856_IG_1o2
    LL_a840_IG_1p
    16_20_24_a840_IG_1p
    24_a840_IG_2p
    24_a850_IG_1p
    LL_a860_IG_1p
    24_32_a860_IG_1p
    """)
        raise Exception('From errorflag: this lattice properties are not in here')

    # shorthand for back compatibility
    xi=chi
