#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2022-02-01 16:13:13
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 2.0
# Macros to run read/write jackknife files

import csv
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import re

def type_Jack_file(filename):
    with open(filename, 'r') as f:
        # Fill this list with all our levels
        data = csv.reader(f, delimiter=' ') #change from default comma

        row = data.__next__()


        cfgs, tl, comp = [int(row[0]), int(row[1]), int(row[2])] 

        if comp != 0 and comp != 1:
            print("Not prepared for {0}".format(row))
            print("Format should be: Ncfg Nt (0:r 1:c) 0 1")
            raise ValueError

    return cfgs, tl, comp


def read_Jack_file(filename, r_comp=0):
    """
    Read a Jack file from Roberts format: Ncfg Nt (0:r 1:c) 0 1
    return cfgs, npoints, [xdata, ydata]
    r_comp {
        0 : real data
        1 : comp data real part
        2 : comp data imag part
    }
    """

    imag=0
    if r_comp== 2:
        r_comp=1
        imag=1

    with open(filename, 'r') as f:
        # Fill this list with all our levels
        data = csv.reader(f, delimiter=' ') #change from default comma

        dataarray = []
        count = 0
        # thiscfgdata = []

        for nn, row in enumerate(data):
            # print(row)
            if nn == 0:
                cfgs, tl, comp = [int(row[0]), int(row[1]), int(row[2])] 
                print("There are {0} configs, with {1} timeslices each".format(cfgs,tl))
                if comp != 0 and comp != 1:
                    print("Not prepared for {0}".format(row))
                    print("Format should be: Ncfg Nt (0:r 1:c) 0 1")
                    raise ValueError
                try:
                    if comp != r_comp:
                        err = []
                        err.append("The type of file does not match the requested read")
                        err.append("File: {0}, requested {1}".format(comp,r_comp))
                        err.append("Use r_comp={0:real,1:real part of complex,2:imag part of complex")
                        # print("\n".join(err))
                        raise ValueError("\n".join(err))
                except ValueError as e:
                    if r_comp != 0:
                        raise e
                    print(e)
                    print("!!!!! WARNING: reading real part by default")


            else:
                if count == 0:
                    thiscfgdata = []
                thiscfgdata.append([int(row[-2-comp]),float(row[-1-comp+imag])])
                count += 1

                if count == tl:
                    dataarray.append(thiscfgdata)
                    count -= tl


    return cfgs, tl, dataarray


def write_Jack_file(newfile, dataarray, w_comp=0):
    """
    w_comp {
        0 : real data
        1 : comp data
    }
    """

    cfgs, xlen = np.shape(dataarray)

    xlist = np.arange(xlen)

    with open(newfile, 'w') as f:
        if w_comp:
            print("Writing a complex valued ensemble in: " + newfile)
        else:
            print("Writing a real valued ensemble in: " + newfile)
            
        data = csv.writer(f, delimiter=' ') #change from default comma
        data.writerow([cfgs, xlen, w_comp, 0, 1])

        for datacfg in dataarray:
            if w_comp:
                data_w_x = [[x,data.real, data.imag] for x,data in zip(xlist,datacfg)]
            else:
                data_w_x = [[x,data] for x,data in zip(xlist,datacfg)]
            
            for row in data_w_x:
                data.writerow(row)

    return newfile

def prune_Jack_file(filename, pruned_cfgs, r_comp=0):
    """
    Read a Jack file from Roberts format: Ncfg Nt (0:r 1:c) 0 1
    Write it again without the pruned cfgs
    """
    if filename[-5:] == '.jack':
        newfile = filename[:-5] + '_pruned.jack'
    else:
        newfile + '.pruned'

    cfgs, tl, dataarray = read_Jack_file(filename, r_comp)

    newcfgs = cfgs - len(pruned_cfgs)

    with open(newfile, 'w') as f:
        # Fill this list with all our levels
        data = csv.writer(f, delimiter=' ') #change from default comma


        # thiscfgdata = []
        print("Writing a real valued ensemble")

        data.writerow([newcfgs, tl, 0, 0, 1])

        for nn, datacfg in enumerate(dataarray):
            np1 = nn + 1
            if np1 in pruned_cfgs:
                continue

            for row in datacfg:
                data.writerow(row)

    print("Changed to {0} configs, with {1} timeslices each.".format(newcfgs, tl))
    return newfile


def get_mask_from_noise(filename,nrcutoff,r_comp=0):
    """
    provide a noise ratio cutoff to create a mask for timeslices surpassing that
    """

    cfgs, npoints, rawdata = read_Jack_file(filename, r_comp)

    # Assume x data is fixed
    print("Assuming the xdata does not have configuration variation")
    xdata = [rawdata[0][tt][0] for tt in range(npoints)]

    ydata = np.array(rawdata)[:,:,1]

    m_ydata = meanense(ydata)

    e_ydata = errormean(ydata, m_ydata)

    noiseratio = e_ydata/m_ydata

    mask = []
    for x, nr in zip(xdata,noiseratio):

        if nr > nrcutoff:
            mask.append(x)

    return mask

def sort_ensem_with_xd(xdata, ensem):
    """
    Sort ensemble by increasing values of xdata

    Parameters
    ----------
    xdata : array / list

    ensem : array (cfgs, len(xdata))
        Ensemble of data

    Returns 
    -------
    xdata, ensem :
        sorted arrays
    """

    sorted_indices = np.argsort(xdata)

    sorted_ensem = []
    sorted_xdata = []

    for ind in sorted_indices:
        sorted_ensem.append(ensem[:, ind])
        sorted_xdata.append(xdata[ind])


    return sorted_xdata, np.transpose(np.array(sorted_ensem))



def mask_xd_ensem(xdata, ensem, mask=[]):
    """
    Mask an xd and an ensemble according to a boolean mask array
    True elements stay, False elements are masked away
    Parameters
    ----------
    xdata : array

    ensem : array (cfgs, len(xdata))
        Ensemble of data

    mask : list
        boolean values to be masked

    Returns
    -------
    xdata, ydata, xmasked, ymasked :
        array and ensemble
    """
    if len(mask) == 0: # masking of the data
            return xdata, ydata, [], []


    xmasked = []
    ymasked = []

    xnew = xdata.copy()
    ensemnew = np.copy(ensem)

    lenxdata = len(xdata)
    npoints = lenxdata

    for nn, keep in enumerate(mask[-1::-1]): # need to go in reverse as we will remove vals
        if keep: # no need to remove
            continue

        torem = lenxdata - nn - 1

        # print(nn, torem, len(xdata), xdata[torem])
        xmasked.append(xdata[torem])
        ymasked.append(ensem[:,torem])

        xnew.pop(torem)
        newydata = np.delete(ensemnew,(torem), axis=1)
        ensemnew = newydata

        npoints-=1

    print(f"There are {npoints} data values left after masking.")

    return xnew, ensemnew, xmasked[-1::-1], np.transpose(np.array(ymasked[-1::-1]))

def maskdata(filename, mask=[], r_comp=0):
    """ take File and return cfgs, npoints, xdata, ydata(ensemble), xmasked, ymasked(ensemble)
    r_comp {
        0 : real data
        1 : comp data real part
        2 : comp data imag part
    }
    """
    cfgs, npoints, rawdata = read_Jack_file(filename, r_comp)

    # Assume x data is fixed
    print("Assuming the xdata does not have configuration variation")
    xdata = [rawdata[0][tt][0] for tt in range(npoints)]

    ydata = np.array(rawdata)[:,:,1]

    xmasked = []
    ymasked = []

    if len(mask) != 0: # masking of the data
        for elem in mask:
            try:
                torem = xdata.index(elem)
                
                xmasked.append(xdata[torem])
                ymasked.append(ydata[:,torem])

                xdata.pop(torem)
                newydata = np.delete(ydata,(torem), axis=1)
                ydata = newydata
                npoints-=1
            except Exception as e:
                print("Error in mask: ", mask, ", if an element is twice it was removed only once")
                print(e)

        print(f"{npoints} timeslices left after masking.")
    
    return cfgs, npoints, xdata, ydata, xmasked, np.transpose(np.array(ymasked))

def get_data(filename, r_comp=0):
    """ 
    take File and return cfgs, npoints, xdata, ydata(ensemble)
    r_comp {
        0 : real data
        1 : comp data real part
        2 : comp data imag part
    }
    """
    cfgs, npoints, rawdata = read_Jack_file(filename, r_comp)

    # Assume x data is fixed
    print("Assuming the xdata does not have configuration variation")
    xdata = [rawdata[0][tt][0] for tt in range(npoints)]

    ydata = np.array(rawdata)[:,:,1]
    
    return cfgs, npoints, xdata, ydata

def get_ensem(filename, r_comp=0):
    """ 
    take File and return ydata(ensemble)
    r_comp {
        0 : real data
        1 : comp data real part
        2 : comp data imag part
    }
    """
    cfgs, npoints, rawdata = read_Jack_file(filename, r_comp)

    return np.array(rawdata)[:,:,1]

###################
###################
# Jackknife functions that assume that the axis = 0 is over the sampling variation
###################
###################

def jackdown(ensemble):
    ensemble = np.array(ensemble)
    avg = np.mean(ensemble, axis=0)
    
    f = - (ensemble.shape[0] - 1)
    
    return avg + (ensemble - avg)/f

def jackup(ensemble):
    ensemble = np.array(ensemble)
    avg = np.mean(ensemble, axis=0)
    
    f = - (ensemble.shape[0] - 1)
    
    return avg + (ensemble - avg) * f

def meanense(ense):
    return  np.mean(ense, axis = 0)

###################
# Assume the ense is two dimensional for the covariance matrix, implementation for more dimensions not immediate
# This is the covariance of the mean, not the intrinsic covariance of the data
###################
def covmatense(ense, meanense, jackknifed = False): 
    """This is the covariance of the mean, not the intrinsic covariance of the data"""
    ense = np.array(ense)
    indepindices = ense.shape[1:]
    
    if len(indepindices) > 1:
        raise Exception("This only calculates covariance matrix of ensemble with one dimension")
        
    
    sizeense = ense.shape[0]
    sizedim = indepindices[0]
    sigma = np.empty([sizedim,sizedim])
    
    for ii in range(sizedim):

        for jj in range(ii, sizedim):

            sigma[ii,jj] = np.sum( (ense[:,ii] - meanense[ii])*(ense[:,jj] - meanense[jj]), axis=0)
            
            if jackknifed:
                sigma[ii,jj] *= (sizeense-1)/(sizeense)
                
            else:            
                sigma[ii,jj] /= (sizeense*(sizeense-1))
            
            if ii != jj:
                sigma[jj,ii] = sigma[ii,jj]
                                
    return sigma


def cormatense(ense, meanense, jackknifed = False):
    """Calculate the correlation of the data"""

    covmat = covmatense(ense, meanense, jackknifed)

    return cov2cor(covmat)

def cov2cor(cov):
    sig = np.sqrt(np.diag(cov))

    def check_zeros(s1_s2):
        if s1_s2 == 0:
            return 1
        else:
            return 1/s1_s2

    sig2weight = np.array([[check_zeros(s1*s2) for s1 in sig] for s2 in sig])

    return cov * sig2weight # this does element-wise multiplication of the arrays   


def errormean(ense, meanense, jackknifed = False):
    """This is the error of the mean, not the intrinsic covariance of the data"""                
    #return  np.diag(covmatense(ense, meanense, jackknifed))**(1/2) # This is a waste of resources no?

    ense = np.array(ense)
    indepindices = ense.shape[1:]
    
    if len(indepindices) > 1:
        raise Exception("This only calculates covariance matrix of ensemble with one dimension")
        
    sizeense = ense.shape[0]
    sizedim = indepindices[0]
    sigma = np.empty([sizedim])
    
    for ii in range(sizedim):

        sigma[ii] = np.sum( (ense[:,ii] - meanense[ii])*(ense[:,ii] - meanense[ii]), axis=0)
            
        if jackknifed:
            sigma[ii] *= (sizeense-1)/(sizeense)
            
        else:            
            sigma[ii] /= (sizeense*(sizeense-1))
                                
    return (sigma)**(1/2)

def calc(ense, jackknifed = False):
    """ Mock the calc command line utility
    Return an array with shape (indep_index, 2)
    """

    mean = meanense(ense)
    err = errormean(ense, mean, jackknifed)

    return np.transpose([mean, err])

def get_mean_error_array(list_of_JK_ensems):

    me_arr = []

    for JK in list_of_JK_ensems:

        mean = meanense(JK)

        error = errormean(JK,mean)

        me_arr.append([mean[0], error[0]])

    return np.array(me_arr)

def combine_ensemble_list(list_of_ensems):
    """ Take a list of 1-dim arrays, eg from get_data,
    and put them into an ensemble array shape (confs, indep_index)
    """
    cfgs = np.shape(np.array(list_of_ensems))[1]
    print("There are {1} data values, with {0} configs each.".format(cfgs,len(list_of_ensems)))
    return np.transpose(np.array(list_of_ensems)[:,:,0])

def change_mean_error(ensem, mean=None, error=None, jackknifed = False):
    """ Change the mean and/ or error of an ensemble
    the procedure is identical whether the ensemble is scaled up or down
    but add the jackknifed option for consistency with other macros
     """


    if mean == None and error == None:
        return ensem

    orig_mean = meanense(ensem)
    orig_err = errormean(ensem, orig_mean)

    if mean == None:
        new_mean = orig_mean
    else:
        new_mean = mean

    if error == None:
        ratio_err = 1
    else:
        ratio_err = error/orig_err


    return ratio_err*(ensem - orig_mean) + new_mean




def ensemble_op(fun, ensemble, jackknifed = False):
    """Do an operation over an ensemble with jackknife statistics
    always return the upscale ensemble, but can receive a down scaled ensemble
    """

    if not jackknifed:
        ensedown = jackdown(ensemble)
    else:
        ensedown = ensemble

    enseres = fun(ensedown)

    return jackup(enseres)

def two_ensemble_op(operation, ensemble1, ensemble2, jackknifed = False):
    """Do an operation(ens1, ens2) with jackknife statistics
    always return the upscale ensemble, but can receive a downscaled ensembles """

    if not jackknifed:
        ensedown1 = jackdown(ensemble1)
        ensedown2 = jackdown(ensemble2)
    else:
        ensedown1 = ensemble1
        ensedown2 = ensemble2

    op_ense = operation(ensedown1, ensedown2)

    return jackup(op_ense)


def td_ensemble_op(td_fun, operation, xd, ensemble, jackknifed = False):
    """Do a time-dependant operation(fun(t), ens(t)) over an ensemble with jackknife statistics
    always return the upscale ensemble, but can receive a downscaled ensemble 
    
    Example of a td_fun
    def td_fun(t):
        mass = .2612
        t0 = 10
        
        return np.exp(mass*(t-t0))
    
    and to multiply it by the jackknife you would do
    td_ensemble_op(td_fun, lambda x,y: x*y, xd, ensemble)

    """

    if not jackknifed:
        ensedown = jackdown(ensemble)
    else:
        ensedown = ensemble

    td_factor = [td_fun(t) for t in xd]

    cfgs = np.shape(ensedown)[0]
    
    factorjk = np.array([td_factor for i in range(cfgs)])

    op_ense = operation(factorjk, ensedown)

    return jackup(op_ense)
        

def eff_mass(correl, dt = 3, jackknifed = False):
    """
    Calculate the jackknife ensemble of the effective mass of a correlator

    """  
    cc = correl[:,:-dt]
    ccdt = correl[:,dt:]
    effmass = two_ensemble_op(lambda x,y: np.log(np.abs(x/y))/dt, cc, ccdt, jackknifed)
    return effmass


def mprint(mat,round=100):
    """Print a matrix rounding to 1/round parts"""
    print(np.round(round*mat)/round)


###################
###################
# Plot functions
###################
###################
def dummy_factor(x):
    return 1

def plot_line_model(ax, lab, nn, xd, model, params_ense, factor=dummy_factor):
    """
    Plot the value of a model
    input: axes, label, nn: color number, xd: min and max to get range, model function, params ensemble, [factor in xspace]

    """

    xarr = np.linspace(min(xd), max(xd), num=100)

    paramsjk = jackdown(params_ense)

    pred = []
    for paramset in paramsjk:
        pred.append([model(xval, *paramset) * factor(xval) for xval in xarr])

    predup = jackup(pred)

    mean_pred = meanense(predup)
    error_pred = errormean(predup, mean_pred)


    out2 = ax.fill_between(xarr, mean_pred-error_pred, mean_pred+error_pred, alpha=0.5, facecolor = 'C'+ str(nn))
    out1 = ax.plot(xarr, mean_pred, c = 'C'+str(nn),label=lab, lw=1, zorder=1)

    return out1, out2

def plot_line_model_boot(ax, lab, nn, xd, model, params_ense, factor=dummy_factor):
    """
    Plot the value of a model with bootstrap resampling
    input: axes, label, nn: color number, xd: min and max to get range, model function, params ensemble, [factor in xspace]

    """

    xarr = np.linspace(min(xd), max(xd), num=100)

    pred = []
    for paramset in params_ense:
        pred.append([model(xval, *paramset) * factor(xval) for xval in xarr])


    mean_pred = meanense(pred)
    error_pred = np.std(pred, axis=0, ddof=1)


    out2 = ax.fill_between(xarr, mean_pred-error_pred, mean_pred+error_pred, alpha=0.5, facecolor = 'C'+ str(nn))
    out1 = ax.plot(xarr, mean_pred, c = 'C'+str(nn),label=lab, lw=1, zorder=1)

    return out1, out2

def plot_line_model_color(ax, lab, xd, model, params_ense, color, factor=dummy_factor):
    """
    Plot the value of a model
    input: axes, label, xd: min and max to get range, model function, params ensemble, color_kw, [factor in xspace]

    """

    xarr = np.linspace(min(xd), max(xd), num=100)

    paramsjk = jackdown(params_ense)

    pred = []
    for paramset in paramsjk:
        pred.append([model(xval, *paramset) * factor(xval) for xval in xarr])

    predup = jackup(pred)

    mean_pred = meanense(predup)
    error_pred = errormean(predup, mean_pred)


    out = ax.fill_between(xarr, mean_pred-error_pred, mean_pred+error_pred, alpha=0.5, facecolor = color)
    out = ax.plot(xarr, mean_pred, c = color,label=lab, lw=1, zorder=1)

    return out

def plot_error_model(ax, lab, nn, xd, model, params_ense, factor=dummy_factor):
    """
    Plot the value of a model
    input: axes, label, nn: color number, xd: min and max to get range, model function, params ensemble, [factor in xspace]

    """

    xarr = np.linspace(min(xd), max(xd), num=100)

    paramsjk = jackdown(params_ense)

    pred = []
    for paramset in paramsjk:
        pred.append([model(xval, *paramset) * factor(xval) for xval in xarr])

    predup = jackup(pred)

    mean_pred = meanense(predup)
    error_pred = errormean(predup, mean_pred)


    out = ax.fill_between(xarr, -error_pred, error_pred, alpha=0.5, facecolor = 'C'+str(nn))
    # out = ax.plot(xarr, mean_pred, c = 'C'+str(nn),label=lab, lw=1, zorder=1)

    return out

def plot_line_model_mean(ax, lab, nn, xd, model, params_ense, factor=dummy_factor):
    """
    Plot the value of a model
    input: axes, label, nn: color number, xd: min and max to get range, model function, params ensemble, [factor in xspace]

    """

    xarr = np.linspace(min(xd), max(xd), num=100)

    paramsjk = jackdown(params_ense)

    pred = []
    for paramset in paramsjk:
        pred.append([model(xval, *paramset) * factor(xval) for xval in xarr])

    predup = jackup(pred)

    mean_pred = meanense(predup)

    out = ax.plot(xarr, mean_pred, c = 'C'+str(nn),label=lab, lw=0.6, zorder=1)

    return out

def print_line_model(xd, model, params_ense, factor=dummy_factor):
    """
    Print the value of a model
    input: xarray, model function, params ensemble, [factor in xspace]

    """

    xarr = np.linspace(min(xd), max(xd), num=105)

    paramsjk = jackdown(params_ense)

    pred = []
    for paramset in paramsjk:
        pred.append([model(xval, *paramset) * factor(xval) for xval in xarr])

    predup = jackup(pred)

    mean_pred = meanense(predup)
    error_pred = errormean(predup, mean_pred)

    print("Model")
    print(np.transpose(np.array([xarr,mean_pred])))
    print("Model + S")
    print(np.transpose(np.array([xarr,mean_pred+error_pred])))
    print("Model - S")
    print(np.transpose(np.array([xarr,mean_pred-error_pred])))

def get_line_model(xd, model, params_ense, factor=dummy_factor):
    """
    Return the value of a model at xd from the ensem of parameters
    input: xarray, model function, params ensemble, [factor in xspace]

    """

    xarr = np.asarray(xd)

    paramsjk = jackdown(params_ense)

    pred = []
    for paramset in paramsjk:
        pred.append(model(xarr, *paramset))

    # this always returns an upscaled ensemble
    predup = td_ensemble_op(factor, lambda x,y: x*y, xarr, pred, jackknifed = True)

    mean_pred = meanense(predup)
    error_pred = errormean(predup, mean_pred)

    return mean_pred, error_pred

def plot_fromensem(axs, xdata, ensem, nn=0, mask = []):
    """ Get the data from a ensemble and plot:
    The mean with errors
    The error to value ratio
    The correlation and covariance
    option to mask the input data

    Parameters
    ----------
    axs : list
        List of Axes where the plots are drawn

    xdata : array
       The x data

    ensem : array (cfgs, len(xdata))
       The ensemble to plot

    nn : integer
       The color in the cycler for ax1

    mask : list
        x-values to be masked

    Returns
    -------
    out : list
        list of artists added to ax1
    """

    if len(axs) == 0:
        print("Axis list is empty, no plots")
        return 1

    # Check if the ensem and xdata array are compatible
    if len(xdata) != np.shape(ensem)[1]:
        raise Exception(f"The sizes of xdata and ensem are not compatible: {len(xdata)}, {np.shape(ensem)[1]}")

    m_data = meanense(ensem)

    e_data = errormean(ensem, m_data)

    axs[0].set_title('Mean Data')

    out = plot_data(axs[0], xdata, m_data, e_data, nn)

    if len(axs) > 1:

        axs[1].plot(xdata, np.abs(e_data/m_data), 'x', c = 'C'+str(nn))
        axs[1].axhline(0,lw=1,color='k')
        axs[1].axhline(0.12,lw=1,color='k')

        xlims = axs[1].get_xlim()
        ylims = axs[1].get_ylim()

        axs[1].set_title('|err/val| ratio')

    if len(axs) > 2:

        corr = cormatense(ensem, m_data)

        correlation_plot(axs[2], xdata, corr)
        
        axs[2].set_title('Correlation matrix')


        if len(axs) > 3:
            cov = covmatense(ensem, m_data)

            matrix_plot(axs[3], xdata, np.log(np.abs(cov)))

            axs[3].set_title('Log |Covariance matrix|')

    return out


def plot_fromfile(filename, axs, nn=0, mask = []):
    """ Get the data from a file and plot:
    The mean with errors
    The error to value ratio
    The correlation and covariance
    option to mask the input data

    Parameters
    ----------
    filename : string
        The name of the file with the data

    axs : list
        List of Axes where the plots are drawn

    nn : integer
       The color in the cycler for ax1

    mask : list
        x-values to be masked

    Returns
    -------
    out : list
        list of artists added to ax1
    """

    if len(axs) == 0:
        print("Axis list is empty, no plots")
        return 1

    cfgs, npoints, xdata, ydata, xmasked, ymasked = maskdata(filename,mask )

    print("---")
    print("xdata, shape(ydata):")
    print(len(xdata),np.shape(ydata))
    print("---")

    m_ydata = meanense(ydata)

    e_ydata = errormean(ydata, m_ydata)

    axs[0].set_title('Mean Data')

    out = plot_data(axs[0], xdata, m_ydata, e_ydata, nn)

    xlims = axs[0].get_xlim()
    ylims = axs[0].get_ylim()

    if len(xmasked) > 0:

        m_ymasked = meanense(ymasked)

        e_ymasked = errormean(ymasked, m_ymasked)

        plot_data(axs[0], xmasked, m_ymasked, e_ymasked, nn, alpha=0.5)

        axs[0].set_xlim(xlims)
        axs[0].set_ylim(ylims)

    if len(axs) > 1:

        axs[1].plot(xdata, np.abs(e_ydata/m_ydata), 'x', c = 'C'+str(nn))
        axs[1].axhline(0,lw=1,color='k')
        axs[1].axhline(0.12,lw=1,color='k')

        xlims = axs[1].get_xlim()
        ylims = axs[1].get_ylim()


        if len(xmasked) > 0:

            axs[1].plot(xmasked, np.abs(e_ymasked/m_ymasked), 'x', c = 'C'+str(nn), alpha=0.5)

            axs[1].set_xlim(xlims)
            axs[1].set_ylim(ylims)

        axs[1].set_title('|err/val| ratio')

    if len(axs) > 2:

        corr = cormatense(ydata, m_ydata)

        correlation_plot(axs[2], xdata, corr)
        
        axs[2].set_title('Correlation matrix')


        if len(axs) > 3:
            cov = covmatense(ydata, m_ydata)

            matrix_plot(axs[3], xdata, np.log(np.abs(cov)))

            axs[3].set_title('Log |Covariance matrix|')



    return out

def corr_mat_from_file(filename, axs, mask = [], label_all=False):
    """ Get the data from a file and plot
    The correlation
    option to mask the input data

    Parameters
    ----------
    filename : string
        The name of the file with the data

    axs : plt axis
        Axis where the plot is drawn

    mask : list
        x-values to be masked

    Returns
    -------
    out : correlation matrix
    """

    try:
        cfgs, npoints, xdata, ydata, xmasked, ymasked = maskdata(filename,mask )
    except ValueError as e:
        if str(e).split("\n")[1] == 'File: 1, requested 0':
            print("Data is complex, real part will be read")
            cfgs, npoints, xdata, ydata, xmasked, ymasked = maskdata(filename,mask, r_comp=1 )
        else:
            raise e

    print("---")
    print("xdata, shape(ydata):")
    print(str(len(xdata))+",",np.shape(ydata))
    print("---")

    m_ydata = meanense(ydata)

    corr = cormatense(ydata, m_ydata)

    cmap = mpl.cm.RdBu_r
    norm = mpl.colors.Normalize(vmin=-1, vmax=1)

    matrix_plot(axs, xdata, corr, cmap=cmap, norm=norm, label_all=label_all)
    
    axs.set_title('Correlation matrix')

    return corr


def matrix_plot(ax, xd, matrix, cmap=None, norm=None, label_all=False):
    """
    A helper function to make a matrix plot

    Parameters
    ----------
    ax : Axes
        The axes to draw to

    xd : array
       The x data

    matrix : 2d array
       The matrix to plot

    cmap : colormap
        defaults to None

    norm : Normalize
        defaults to None

    If you need to rotate labels try:
    ax.set_xticklabels([str(xd[t]) for t in ax.get_xticks()] rotation=45, ha='left')

    Returns
    -------
    out : list
        list of artists added
    """

    mat = ax.matshow(matrix, cmap = cmap, norm = norm)
    
    plt.colorbar(mat, ax=ax)

    # Assume xdata is nicely ordered
    if len(xd) < 6 or label_all:
        dist = 1
    else:
        dist = np.floor(len(xd)/6) 
    ticks = range(0, len(xd), int(dist))
    tl = [str(xd[t]) for t in ticks]

    ax.set_xticks(ticks)
    ax.set_xticklabels(tl)
    ax.set_yticks(ticks)
    ax.set_yticklabels(tl)

    return mat


def correlation_plot(ax, xd, corr, label_all=False):
    """
    A helper function to make a correlation plot

    Parameters
    ----------
    ax : Axes
        The axes to draw to

    xd : array
       The x data

    corr : 2d array
       The correlation matrix to plot

    Returns
    -------
    out : list
        list of artists added
    """

    cmap = mpl.cm.RdBu_r
    norm = mpl.colors.Normalize(vmin=-1, vmax=1)

    out = matrix_plot(ax, xd, corr, cmap=cmap, norm=norm, label_all=label_all)

    return out

def plot_data(ax, xd, yd, yerr, nn, **kwargs):
    """
    A helper function to make a graph

    Parameters
    ----------
    ax : Axes
        The axes to draw to

    xd : array
       The x data

    yd : array
       The y data

    ye : array
       The y data error

    nn : integer
       The color in the cycler

    Returns
    -------
    out : list
        list of artists added
    """


    # Assume the data and errors were calculated 

    out = ax.errorbar(xd, yd, yerr,
                 ls='', c = 'C'+str(nn),
                 marker='s',mfc='w',#fillstyle='none',
                 ms=4, 
                 elinewidth=1, capsize=3.5, **kwargs)
    return out

def plot_data_xerr(ax, xd, xerr, yd, yerr, nn, **kwargs):
    """
    A helper function to make a graph

    Parameters
    ----------
    ax : Axes
        The axes to draw to

    xd : array
       The x data

    yd : array
       The y data

    ye : array
       The y data error

    nn : integer
       The color in the cycler

    Returns
    -------
    out : list
        list of artists added
    """


    # Assume the data and errors were calculated 

    out = ax.errorbar(xd, yd, yerr, xerr,
                 ls='', c = 'C'+str(nn),
                 marker='.',mfc='w',ms=4, 
                 elinewidth=1, capsize=3.5, **kwargs)
    return out


def plot_ensemble_mean_err(ax, xd, ensem, nn=0, **kwargs):
    """
    A helper function to plot the mean and error of an ensemble

    Parameters
    ----------
    ax : Axes
        The axes to draw to

    xd : array
       The x data

    ensem : array (cfgs, len(xd))
       The ensemble to plot as y data and error

    nn : integer
       The color in the cycler

    Returns
    -------
    out : list
        list of artists added
    """   

    # Check if the ensem and xd array are compatible
    if len(xd) != np.shape(ensem)[1]:
        raise Exception(f"The sizes of xd and ensem are not compatible: {len(xd)}, {np.shape(ensem)[1]}")

    mean_err = calc(ensem)
    yd = mean_err[:,0]
    yerr = mean_err[:,1]


    return plot_data(ax, xd, yd, yerr, nn, **kwargs)

def plot_eff_mass(ax, correl, **kwargs):
    """
    A helper function to plot the eff mass of correlator

    Parameters
    ----------
    ax : Axes
        The axes to draw to

    correl : correlator
       The y data

    Returns
    -------
    out : artist
        artist added to ax
    """

    effmass = eff_mass(correl)

    mean_error = calc(effmass)

    t = np.arange(np.shape(mean_error)[0])

    out = ax.errorbar(t, mean_error[:,0], mean_error[:,1],
                 ls='', elinewidth=1, capsize=3.5, **kwargs)
    return out


def plot_cfg_fromfile(filename, axs, mask = [], scaled = True):
    """
    Input: jackfile, list of axs, only the first one will be used.
    Plot the individual configurations, and the contours over several standard deviations
    Outliers should be easy to spot with this visualization.
    the scaled option is to subtract the mean and scale by the error
    """

    cfgs, npoints, xdata, ydata, xmasked, ymasked = maskdata(filename,mask )

    print("---")
    print("xdata, shape(ydata):")
    print(len(xdata),np.shape(ydata))
    print("---")

    if len(axs) == 0:
        print("Axis list is empty, no plots")
        return 1

    m_ydata = meanense(ydata)

    e_ydata = errormean(ydata, m_ydata)

    if len(xdata) == 1:
        if scaled:
            axs[0].hist((ydata - m_ydata)/(e_ydata*np.sqrt(cfgs)), bins=int(cfgs/10))
        else:
            axs[0].hist(ydata, bins=int(cfgs/10))
        return axs[0]

    cmap = mpl.cm.gist_rainbow
    cmap = cmap.reversed()
    norm = mpl.colors.Normalize(vmin=1, vmax=cfgs)

    gcmap=plt.cm.get_cmap('Greys', 10)
    gind = [i/10 for i in range(1,6)]
    gind.extend(np.arange(5,0,-1)/10)
    cmapc = mpl.colors.ListedColormap([gcmap(xx/.63) for xx in gind])

    for nn, ycfg in enumerate(ydata):
        if scaled:
            ycfg = (ycfg - m_ydata)/(e_ydata*np.sqrt(cfgs))

        axs[0].plot(xdata, ycfg, '_', ls='', color = cmap( norm( nn + 1 )))

    xlims = axs[0].get_xlim()
    ylims = axs[0].get_ylim()


    Ny = 100

    y = np.linspace(ylims[0], ylims[1], num= Ny)
    X, Y = np.meshgrid(xdata, y)
    if scaled:
        Z = Y
    else:
        Z = (Y - m_ydata )/(e_ydata*np.sqrt(cfgs))

    CS = axs[0].contourf(X, Y, Z,levels=np.arange(-9,11,2), cmap=cmapc)
    CS1 = axs[0].contour(CS,colors='k',linestyles= 'solid', linewidths=.8)
    labs = axs[0].clabel(CS1, inline=True, fontsize=10)

    axs[0].set_xlim(xlims)

    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("bottom", size="3%", pad=0.3)

    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                 cax=cax, label='Configuration', orientation = 'horizontal')

    return axs[0]

def plot_cfg_fromens(ensem, ax, scaled = True):
    """
    Input: ensemble, an ax.
    Plot the individual configurations, and the contours over several standard deviations
    Outliers should be easy to spot with this visualization.
    the scaled option is to subtract the mean and scale by the error
    """

    cfgs = len(ensem[:,0])
    npoints = len(ensem[0,:])
    xdata = np.arange(0,npoints)
    ydata = ensem

    print("---")
    print("xdata, shape(ydata):")
    print(len(xdata),np.shape(ydata))
    print("---")


    m_ydata = meanense(ydata)

    e_ydata = errormean(ydata, m_ydata)

    cmap = mpl.cm.gist_rainbow
    cmap = cmap.reversed()
    norm = mpl.colors.Normalize(vmin=1, vmax=cfgs)

    gcmap=plt.cm.get_cmap('Greys', 10)
    gind = [i/10 for i in range(1,6)]
    gind.extend(np.arange(5,0,-1)/10)
    cmapc = mpl.colors.ListedColormap([gcmap(xx/.63) for xx in gind])

    for nn, ycfg in enumerate(ydata):
        if scaled:
            ycfg = (ycfg - m_ydata)/(e_ydata*np.sqrt(cfgs))

        ax.plot(xdata, ycfg, '_', ls='', color = cmap( norm( nn + 1 )))

    xlims = ax.get_xlim()
    ylims = ax.get_ylim()


    Ny = 100

    y = np.linspace(ylims[0], ylims[1], num= Ny)
    X, Y = np.meshgrid(xdata, y)
    if scaled:
        Z = Y
    else:
        Z = (Y - m_ydata )/(e_ydata*np.sqrt(cfgs))

    CS = ax.contourf(X, Y, Z,levels=np.arange(-9,11,2), cmap=cmapc)
    CS1 = ax.contour(CS,colors='k',linestyles= 'solid', linewidths=.8)
    labs = ax.clabel(CS1, inline=True, fontsize=10)

    ax.set_xlim(xlims)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="3%", pad=0.3)

    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                 cax=cax, label='Configuration', orientation = 'horizontal')

    return ax

def add_labels_correlation(covax, filename, mask=[], **kwargs):
    """
    Write values of matrix in a matshow axis when plotting correlation directly from a file
    ax: Axis where matshow was used
    mat: the values of the matrix in question
    kwargs for the text
    """

    cfgs, npoints, xdata, ydata, xmasked, ymasked = maskdata(filename,mask )
    m_ydata = meanense(ydata)

    e_ydata = errormean(ydata, m_ydata)


    corr = cormatense(ydata, m_ydata)

    for (j,i),label in np.ndenumerate(corr):
        if i < j:
            continue
        col = np.abs(label)
        
        # Some aesthetic, we dont need 1.00, 1 is enough, we dont need 0.7 or -0.7, .7 and -.7
        text = f"${label:.2f}$".replace("0.", ".").replace("1.00","1").replace("-.00","0").replace("-0","0")        
        if col < 0.7:
            covax.text(i,j,text,ha='center',va='center', c='k', **kwargs)

        else:
            covax.text(i,j,text,ha='center',va='center', c='w', **kwargs)


def add_labels_matrix(ax, mat, sym=True, hide_diag=False, **kwargs):
    """
    Write values of matrix in a matshow axis
    ax: Axis where matshow was used
    mat: the values of the matrix in question
    sym: boolean to skip the lower triangular possibly redundant labels
    hide_diag: boolean to skip the label over the diagonal of the matrix
    kwargs for the text
    """

    col_thres = np.max(mat)

    for (j,i),label in np.ndenumerate(mat):
        if i < j and sym:
            continue
        if i == j and hide_diag:
            continue
        col = np.abs(label)
        # Some aesthetic, we dont need 1.00, 1 is enough, we dont need 0.7 or -0.7, .7 and -.7
        text = re.sub("^-?0.", lambda x:'.' if x.group(0)=='0.' else '-.', f"{label:.2f}")
        text = text.replace("1.00","1").replace("-.00","0").replace("-0","0")
        if col < 0.7*col_thres:
            ax.text(i,j,text,ha='center',va='center', c='k', **kwargs)

        else:
            ax.text(i,j,text,ha='center',va='center', c='w', **kwargs)
