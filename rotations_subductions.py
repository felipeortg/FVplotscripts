#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2023-05-04 20:47:15
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 1.0
# Rotation conventions, useful for kinematic factor calculations

import spectrum as spec
import numpy as np

def Ry(theta):
    return np.array([[ np.cos(theta), 0, np.sin(theta)],
                   [ 0           , 1, 0           ],
                   [-np.sin(theta), 0, np.cos(theta)]])
  
def Rz(theta):
    return np.array([[ np.cos(theta), -np.sin(theta), 0 ],
                   [ np.sin(theta), np.cos(theta) , 0 ],
                   [ 0           , 0            , 1 ]])

def Euler_rot(phi, theta, psi):
    return Rz(phi) @ Ry(theta) @ Rz(psi)


# This can be used as checks
# Euler_rot(np.pi/2, np.pi/4, -np.pi/2) @ [0,0,np.sqrt(2)]
# Euler_rot(np.pi/4, np.arccos(1/np.sqrt(3)), 0) @ [0,0,np.sqrt(3)]
# Euler_rot(np.pi/2, np.arccos(2/np.sqrt(5)), 0) @ [0,0,np.sqrt(5)]
# Euler_rot(-3*np.pi/4, -np.arccos(np.sqrt(2/3)), 0) @ [0,0,np.sqrt(6)]

sq2 = np.sqrt(2)
pol_vec = {
    "p1":np.array([-1./sq2, -1.j/sq2,0]),
    "0":np.array([0,0,1]),
    "m1":np.array([1./sq2, -1.j/sq2,0]),
    }

# from tabII of 1107.1930
subduction_J1m = {
    #000
    "T1_0" : {"r1": pol_vec['p1'],
              "r2": pol_vec['0'],
              "r3": pol_vec['m1']
             },
    #100
    "A1_1" : {"r1": pol_vec['0']},
    "E2_1" : {"r1": (pol_vec['p1'] + pol_vec['m1'])/sq2,
              "r2": (pol_vec['p1'] - pol_vec['m1'])/sq2
             },
    # 110
    "A1_2" : {"r1": pol_vec['0']},
    "B1_2" : {"r1": (pol_vec['p1'] + pol_vec['m1'])/sq2
             },
    "B2_2" : {"r1": (pol_vec['p1'] - pol_vec['m1'])/sq2
             },
    #111
    "A1_3" : {"r1": pol_vec['0']},
    "E2_3" : {"r1": (pol_vec['p1'] + pol_vec['m1'])/sq2,
              "r2": (pol_vec['p1'] - pol_vec['m1'])/sq2
             },
    #200
    "A1_4" : {"r1": pol_vec['0']},
    "E2_4" : {"r1": (pol_vec['p1'] + pol_vec['m1'])/sq2,
              "r2": (pol_vec['p1'] - pol_vec['m1'])/sq2
             }
}

#from 1203.6041v2 supplemental material
# also save to GDrive/Documents/JLab/software/rotations and CG phases/summary
# to add more use that
rotation_lat = {
    "000": np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    "100": np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]]),
    "110": np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]]),
    "111": np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    "200": np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
}

#from 1203.6041v2 and tab X
rotation_ref = {
    "0": Euler_rot(0,0,0),
    "00n": Euler_rot(0,0,0),
    "0nn": Euler_rot(np.pi/2, np.pi/4, -np.pi/2),
    "nnn": Euler_rot(np.pi/4, np.arccos(1/np.sqrt(3)), 0),
    "0n2n": Euler_rot(np.pi/2, np.arccos(2/np.sqrt(5)), 0),
    "nn2n": Euler_rot(-3*np.pi/4, -np.arccos(np.sqrt(2/3)), 0)
}

def rotation(mom_lab):
    mom_vec = spec.label2vec(mom_lab)

    mom_type = vec2mom_type( mom_vec )

    rot_ref = rotation_ref[mom_type]

    return rotation_lat[mom_lab] @ rot_ref


def print_mom_type_comp(comp):
    if comp == 1:
        return ""
    else:
        return int(comp)

def vec2mom_type(vec):
    pos_vec = np.sort(np.abs(vec))

    # case 000
    if pos_vec[-1] == 0:
        return "0"

    # case 00n
    if pos_vec[-2] == 0:
        return "00n"

    # case 0mn
    if pos_vec[-3] == 0:
        # case 0nn
        if pos_vec[1] == pos_vec[2]:
            return "0nn"
        else:
            div = np.gcd(pos_vec[1], pos_vec[2])
            return f"0{print_mom_type_comp(pos_vec[1]/div)}n{print_mom_type_comp(pos_vec[2]/div)}n"

    # case nnm
    if pos_vec[0] == pos_vec[1]:
        if pos_vec[1] == pos_vec[2]:
            return "nnn"

    # case  nmp
    div = np.gcd.reduce(pos_vec)
    return f"{print_mom_type_comp(pos_vec[0]/div)}n{print_mom_type_comp(pos_vec[1]/div)}n{print_mom_type_comp(pos_vec[2]/div)}n"






