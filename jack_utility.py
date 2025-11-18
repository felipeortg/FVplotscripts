#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2022-02-01 16:13:13
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 2.0
# Macros to run read/write jackknife files

import csv
import numpy as np

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
