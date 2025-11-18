#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2025-10-18
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)

import sys
import jack_utility as JK

if len(sys.argv) < 2:
    raise Exception(f"Usage is: {sys.argv[0]} jack_file")

# -----------------
jack_file = str(sys.argv[1])

ensem = JK.get_ensem(jack_file)

print(JK.calc(ensem))

