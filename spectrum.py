#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2021-01-21 17:03:24
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 1.0
# Some macros to build spectrum plots

import numpy as np
import os
import xml.etree.ElementTree as ET

def label2vec(ks):
    vec = []
    for comp in range(len(ks)):
        vec.append(int(ks[comp]))
    
    return np.array(vec)

def vec2label(ks):
    string = ''
    for comp in ks:
        string += str(np.abs(comp))
    
    return str(string)

def reversemom_inirreps(irreps):
    new_irreps = []
    for irre in irreps:
        d = list(label2vec(irre[:3]))
        d.reverse()

        name = irre[3:]

        new_irreps.append(vec2label(d)  + name)

    return new_irreps

def free_part_en(m1,m2,L,k1,k2,P):
    """ Get the energy for m1, m2, L, momenta 1, 2, target"""
    k1arr = np.array(k1)
    k2arr = np.array(k2)
    Parr = np.array(P)

    k1sq = np.sum(k1arr**2)
    k2sq = np.sum(k2arr**2)
    P2 = np.sum(Parr**2)


    p1en = np.sqrt(m1**2 + k1sq * (2*np.pi/L)**2)
    p2en = np.sqrt(m2**2 + k2sq * (2*np.pi/L)**2)

    return np.sqrt( (p1en + p2en)**2 - P2 * (2*np.pi/L)**2 )

def ecm_prop_unc(energyunc, P, Lchi):
    """ Propagate the uncertainty in lattice energy to the CM frame"""
    energy, unc = energyunc
    
    ecm = np.sqrt(energy**2 - np.dot(P,P)*(2*np.pi/Lchi)**2)
    
    unccm = energy * unc / ecm
    
    return [ecm, unccm]

def latt_prop_unc(energyunc, P, Lchi):
    """ Propagate the uncertainty in cm energy to the lattice frame"""
    energy, unc = energyunc
    
    elatt = np.sqrt(energy**2 + np.dot(P,P)*(2*np.pi/Lchi)**2)
    
    unclatt = energy * unc / elatt
    
    return [elatt, unclatt]

def order_two_meson_list_at_xiL(twomesonlist, irreP, xiL, mesonmasses):
    """Order the list in the dictionary """

    energy_list = []
    for twomesons in twomesonlist:
        energy = free_part_en(mesonmasses[twomesons['part1']], mesonmasses[twomesons['part2']], xiL,
                                       label2vec(twomesons['k1']), label2vec(twomesons['k2']), irreP)

        energy_list.append([energy, twomesons])

    energy_list.sort(key=lambda en_mes: en_mes[0])

    sortedlist = [en_mes[1] for en_mes in energy_list]

    return sortedlist

# The following macros take the key from the free energy spectrum and return a dictionary like
# redstar:  {'part1': 'Kneg_proj0', 'k1': '210', 'part2': 'Kbarneg_proj0', 'k2': '210'}
# scat_devel: {'part1': 'kaon', 'k1': '100', 'part2': 'kbar', 'k2': '200'}

def split_redstar_key(string):
    """ Go from a key to particle names, momenta, and target irrep/momentum"""
    # string format is: flav1_proj#_pnml_Hlirrep__mom1xxflav2_proj#_pnml_Hlirrep__mom2__tF,embF_tirrep,embi__momt
    string_pieces = string.split('xx')

    # string_pieces: [flav1_proj#_pnml_Hlirrep__mom1, flav2_proj#_pnml_Hlirrep__mom2__tF,embF_tirrep,embi__momt]
    if len(string_pieces) > 2:
        raise Exception("This key has multiple xxs: " + string)

    # We know the first particle is given before the xx's
    part1 = string_pieces[0].split('__')
    # part1: [flav1_proj#_pnml_Hlirrep, mom1]
    
    remaining_pieces = string_pieces[1].split('__')
    # remaining_pieces: [flav2_proj#_pnml_Hlirrep, mom2, tF,embF_tirrep,embi, momt]
    
    # The second particle and its momentum begin immediately
    part2 = remaining_pieces[0:2]
    # part2: [flav2_proj#_pnml_Hlirrep, mom2]
    
    # We do not care about the flavor, but only about the target irrep
    targetirrep = remaining_pieces[2].split('_')[1].split(',')[0]
    # irrep: (tF,embF_tirrep,embi) -> [tF,embF, tirrep,embi] -> tirrep
    
    targetmom = remaining_pieces[-1]
    
    # Strip irrep and momentum from particle key
    
    parts = [part1, part2]
    newparts = dict()
    for n, part in enumerate(parts):
        # [flav_proj#_pnml_Hlirrep, mom]
        if len(part) < 2:
            print(part)
            Exception('Why dont I have momentum?')

        if part[1] == '000': #we dont use helicity for states at rest
            fullname = part[0][:-3]
        else:
            # Use the helicity to get rid of irrep
            fullname = part[0].split('_H')[0]

        # fullname: flav_proj#_pnml
        
        partname = fullname[:-5]
        mom_in_name = fullname[-3:]

        # partname: flav_proj#
        # mom_in_name: nml
        
        if mom_in_name == part[1]:
            # No surprise here
            partmom = mom_in_name
        else:
            raise Exception("Why is the " + mom_in_name + " not equal to " + part[1] + " ?" )
            
        lab1 = 'part' + str(n+1)
        momlab = 'k' + str(n+1)
        newparts[lab1] = partname
        newparts[momlab] = partmom
    
    return newparts

def split_scat_devel_key(string ):
    """ Go from a key to particle names, momenta, and target irrep/momentum"""
    # string format is: flav1_irrep_mom,emb__flav2_irrep_mom,emb___embtarget
    string_pieces = string.split('__')

     # string_pieces: [flav1_irrep_mom,emb, flav2_irrep_mom,emb, _embtarget]
    if len(string_pieces) > 3:
        raise Exception("This key has multiple __s: " + string)

    # We know the first particle is given before the __'s
    part1 = string_pieces[0].split('_')
    # part1: [flav1, irrep, mom,emb]
    
    # The second particle and its momentum begin immediately
    part2 = string_pieces[1].split('_')
    # part2: [flav2, irrep, mom,em]
    
    # Strip irrep and momentum from particle key
    
    parts = [part1, part2]
    newparts = dict()

    for n, part in enumerate(parts):
        # [flav, irrep, mom,em]
        if len(part) < 3:
            print(part)
            Exception('Why dont I have momentum?')

        partname = part[0]
        partmom = part[2].split(',')[0]

        # partname: flav
        # partmom: nml
            
        lab1 = 'part' + str(n+1)
        momlab = 'k' + str(n+1)
        newparts[lab1] = partname
        newparts[momlab] = partmom

    return newparts

def print_twomeson_dict(particles):
    return particles['part1'] + ' ' + particles['k1'] + ', ' + particles['part2'] + ' ' + particles['k2']


# These macros read spectrum files produced by various programs

def read_reconfit_file(filename, noise = .1):
    """ 
    Read a spectrum file calculated via reconfit
    Return:
    irrep: string with irrep name
    spectrum: list with [Ecm, Ecm_err]
    t0_str: string of the form t0 = #
    """

    spectrum = []

    with open(filename, 'r') as f:
        # Fill this list with all our levels

        states = False
        for nn, line in enumerate(f):

            # Irrep given in the first line
            if nn == 0:
                irrep = line[:-2] # drop the new line and the G-parity
            if nn == 1:
                t0str = line[5:13]

            # begin reading after state
            if line[:5] == 'state':
                states = True
                continue

            # end reading at empty line
            elif line == '\n':
                states = False
                continue

            if states:

                elems = line.split('| ')

                mass, unc = elems[2].split('+/- ')

                if float(mass) == 0:
                    print("###########\n WARNING: Level has zero mass")
                    print(elems)
                    continue

                if float(unc)/float(mass) > noise: #ignore noise results
                    print("###########\n WARNING: Noisy level ignored")
                    print(line)
                else:
                    spectrum.append([float(mass),float(unc)])


    return irrep, spectrum, t0str

def read_redstar_file(filename):
    """ Read a spectrum file generated by redstar """

    two_mesons = dict()

    with open(filename, 'r') as f:
        # Fill this list with all our operators


        for line in f:

            # Skip lines with not even a target
            if len(line) < 8 : continue 

            target, strparticles = line.split('  ')

            target = target[:-1] #drop G-parity

            particles = split_redstar_key(strparticles.rstrip())

            if target in two_mesons:
                two_mesons[target].append(particles)
            else:
                two_mesons[target] = [particles]

    return two_mesons

def read_free_spectrum_scatdevel(filename):
    """ Read a free spectrum file generated by scatdevel """


    # Fill this dict with all our operators
    two_mesons = dict()

    with open(filename, 'r') as f:

        tmp_ene = 0

        inirrep = False
        for line in f:
            if line[0] == '#' or line[:5] == '---->': # ignore comments and printed announcements
                continue

            if line[0] == '*':
                if inirrep:
                    inirrep = False
                    continue
                else:
                    inirrep = True
                    continue


            if inirrep:     
                irrepparams = line.split()
                
                irrep_name = irrepparams[3].split('=')[1]
                
                irrep_P = irrepparams[2].split('=')[1][1:-2].split(',')


                target = ''.join(irrep_P) + '_' + irrep_name

                # scatdevel prints the momentum reversed
                target= reversemom_inirreps([target])[0]

                continue


            splitted = line.split()
            
            if len(splitted) == 0: #ignore empty lines
                continue
            
            if tmp_ene == splitted[0]:
                if len(splitted) < 3 or splitted[-1] != 'embedding':
                    print(splitted)
                    print('WARNING: There is a degeneracy in Irrep: '+ target 
                    + '\n with particles:' + print_twomeson_dict(particles)
                     + ' and particles ' + splitted[1])

            tmp_ene = splitted[0]

            particles = split_scat_devel_key(splitted[1].rstrip())
            

            if target in two_mesons:
                two_mesons[target].append(particles)
            else:
                two_mesons[target] = [particles]

    return two_mesons

def get_irrep_scat_devel_output(file):
    """
    Get the irrep name from get_finite_volume_spectrum_with_errors output
    Reverse the momentum to canonical
    """

    file = os.path.basename(file)

    name, ext = file.split('.')

    elems = name.split('_')

    dd = elems[1][1:]

    irr = elems[2]

    irrep = dd + '_' + irr

    return reversemom_inirreps([irrep])[0]

def read_interacting_spectrum_scatdevel(filename, error = True):
    """ 
    Read an interacting spectrum file calculated via scatt_devel
    That means from get_finite_volume_spectrum[_with_errors] output
    Return a list with [Vol, Ecm, Ecm_err]
    """

    # Fill this list with all our levels
    spectrum = []


    with open(filename, 'r') as f:
        # Fill this list with all our levels
        for nn, line in enumerate(f):
            elems = line.split(' ')

            # the output has 3 spaces between volume and energy

            if error:
                Lmunc = [int(elems[0]),float(elems[3]),float(elems[4])]
            else:
                Lmunc = [int(elems[0]),float(elems[3])]

            # Lmunc = [float(vals) for vals in strings]
            
            spectrum.append(Lmunc)


    return spectrum

def clean_calc_spectrum(dirty_spec, irrep, erange, chi, Lsize, Lwidth):
    """ 
    Take a spectrum from reconfit and dismiss dirty levels and those out of the range of interest
    Also move around the L-value when overlapping for clarity
    Return list with elements [newL, Ecm, Ecm_err]
    The list is made so that levels remain in same order as input
    """
    temp_list = []
    place = 0


    irreP = label2vec(irrep[:3])

    for level in dirty_spec:
        En = level[0]

        mass, unc = ecm_prop_unc(level, irreP, Lsize*chi)

        if mass + unc > erange[0] and mass - unc < erange[1]: #leave only levels in the energy ranges
        
            temp_list.append([Lsize, mass, unc, place])
            place += 1

            # print(mass + unc, erange[0], mass - unc ,erange[1], mass)

    clean_spec = [None]*len(temp_list)
    # move around levels when overlapping
    while len(temp_list) > 0:

        overlaps = [0]

        for jj in range(1, len(temp_list)): # check overlap in remaining levels
            dist = np.abs(temp_list[0][1]-temp_list[jj][1])
            sumunc = temp_list[0][2]+temp_list[jj][2]

            if dist / sumunc < 1: # add them whenever they are too close
                overlaps.append(jj)

        overlapping_lvls = len(overlaps)
        if overlapping_lvls == 1:
            clean_spec[temp_list[0][3]] = (temp_list[0][:3])

        else:
            plotwidth = overlapping_lvls*Lwidth/2
            Lpositions = np.linspace(Lsize - plotwidth, Lsize + plotwidth, num = overlapping_lvls)
            for nn, L in enumerate(Lpositions):
                pl = overlaps[nn]
                clean_spec[temp_list[pl][3]] = [L, temp_list[pl][1], temp_list[pl][2]]

        overlaps.reverse()

        for ov in overlaps:
            del temp_list[ov]


    return clean_spec


def read_Ecm_ini(filename):
    """ Read an Ecm_data.ini file and return the level numbers in each irrep """

    # Fill this dict number of levels
    levelsinirrep = dict()


    with open(filename, 'r') as f:
        for line in f:
            elems = line.split(' ')

            irrep = elems[1] + '_' + elems[2]

            if '(*)' in line: #starred levels in ecm ini are to be ignored
                continue

            if irrep in levelsinirrep:
                levelsinirrep[irrep].append(int(elems[3]))
            else:
                levelsinirrep[irrep] = [int(elems[3])]


    return levelsinirrep


def read_Ecm_xml(filename):
    """
    Read an xml with the Ecm data,
    Return a dictionary with
    Key: V_irrep
    Value: dictionary: key lvl_num, value [Ecm, Ecm_err]
    """
    # read the xml, binary since ET knows how to take care of that
    with open(filename, 'rb') as xml_file:
        tree = ET.parse(xml_file)

    # this is how the xml in python is read
    root = tree.getroot()

    # will save all read levels in here [Key, Val]
    levels = []

    # specific tree structure of Ecm_xmls
    for elem in root.find('Energy_Levels/elem/Val/DATA'):
        key = dict()
        vals = dict()

        for k in elem.find('Key'):
            key[k.tag] = k.text

        for v in elem.find('Val'):
            vals[v.tag] = v.text

        levels.append([key, vals])


    # create a dictionary with the fitted levels
    # Key is V_irrep, Val is list, elements [V, Ecm, Ecm_err, lvl_num]
    fitted_levels = dict()

    for level in levels:
        irrepP = level[0]['d'].replace(" ","")
        irrep = irrepP +"_"+level[0]['irrep']

        # scattering_devel uses reference momenta...
        vol_irrep = "V" + level[0]['V'] + "_" + reversemom_inirreps([irrep])[0]

        lvl_plot_info = [int(level[0]['level_num']), float(level[1]['value']), float(level[1]['error'])]

        if vol_irrep not in fitted_levels:
            fitted_levels[vol_irrep] = []

        fitted_levels[vol_irrep].append( lvl_plot_info )

    fitted_levels_lvl_dict = dict()
    for vol_irrep in fitted_levels:
        lvl_dict = dict()

        for lvl in fitted_levels[vol_irrep]:
            lvl_dict[lvl[0]] = [lvl[1], lvl[2]]

        fitted_levels_lvl_dict[vol_irrep] = lvl_dict 

    return fitted_levels_lvl_dict


def read_ini_file(filename):
    """ Read an Ecm_data.ini or mel_list file and return a dictionary per irrep """

    # Fill this dict number of levels
    levelsinirrep = dict()


    with open(filename, 'r') as f:
        for line in f:
            elems = line.split(' ')

            irrep = elems[1] + '_' + elems[2]

            if '(*)' in line: #starred levels in ini are to be ignored
                continue

            data = dict()
            data['vol'] = elems[0]
            data['mom'] = label2vec(elems[1])
            data['irre'] = elems[2]
            data['lvl_num'] = elems[3]
            data['file_name'] = elems[4].rstrip()


            if irrep in levelsinirrep:
                levelsinirrep[irrep].append(data)
            else:
                levelsinirrep[irrep] = [data]


    return levelsinirrep

def label_state(state_dict):
    """ return a string labeling a level dictionary """

    # make sure it is canonical momentum
    momentum = state_dict['mom']
    momentum.sort()

    return "V" + state_dict['vol'] + "_" + state_dict['irre'] + "_" + vec2label(momentum[::-1]) + "_#" + state_dict['lvl_num']

