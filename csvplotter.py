#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2019-02-21 09:46:45
# @Author  : Felipe Ortega (felipeortegagama@gmail.com)
# @Version : 1.0

# Give the name of the file
# Give the numbers of columns of the data to plot

import numpy as np
import matplotlib.pyplot as plt
import csv
import sys

if len(sys.argv) < 2:
    print("Usage is {0} filename xcol ycol [yerr]".format(sys.argv[0]))
    exit()

filename = sys.argv[1]

xcol = int(sys.argv[2])

ycol = int(sys.argv[3])

if len(sys.argv) > 4:
    print("printing error bars \n")
    yerr = int(sys.argv[4])

with open(filename, 'r') as f:
    data = csv.reader(f, delimiter=' ') #change delimiter, default is a comma

    next(data, None)  # skip the headers

    try:
        if len(sys.argv) > 4:
            plotdata = np.array([[float(row[xcol]),float(row[ycol]),float(row[yerr])] for row in data])
        else:
            plotdata = np.array([[float(row[xcol]),float(row[ycol])] for row in data])

    except Exception as e:
        plotdata = []
        with open(filename, 'r') as f2:
            datanew = csv.reader(f2, delimiter=' ')
            next(datanew, None)  # skip the headers
            for nn, row in enumerate(datanew):
                try:
                    if len(sys.argv) > 4:
                        plotdata.append([float(row[xcol]),float(row[ycol]),float(row[yerr]),float(row[yerr])])
                    else:
                        plotdata.append([float(row[xcol]),float(row[ycol]),float(row[yerr])])

                except Exception as e:  
                    print('##### WARNING:')
                    if len(sys.argv) > 4:
                        print('Ignoring row {2}: x {0} and y {1} +/- {4}. It reads {3}'.format(xcol, ycol, nn, row, yerr))
                        print(e,row[xcol],row[ycol], row[yerr],'.')
                    else:
                        print('Ignoring row {2}: x {0} and y {1}. It reads {3}'.format(xcol, ycol, nn, row))
                        print(e,row[xcol],row[ycol],'.')

        plotdata = np.array(plotdata)

if len(sys.argv) > 4:
    plt.plot(plotdata[:,0], plotdata[:,1])
    plt.fill_between(plotdata[:,0], plotdata[:,1]-plotdata[:,2],plotdata[:,1]+plotdata[:,2],alpha=0.5)
else:
    plt.plot(plotdata[:,0], plotdata[:,1],'.')


plt.show()

