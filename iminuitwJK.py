#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2022-02-01 16:13:13
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 2.0
# Macros to run read/write jackknife files, perform operations and fits, make visualizations.

# qcdi version: no iminuit (and fitting) support

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import csv
# import subprocess
# from iminuit import Minuit
# from iminuit.util import describe, make_func_code

from packaging import version
# from iminuit import __version__ as imversion
# if version.parse(imversion) < version.parse("2.6"):
#    raise Exception(f"iminuit version 2.6 or newer is required, version {imversion} found.")


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

    print("Changed to {0} configs, with {1} timeslices each".format(newcfgs, tl))
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
# Fit functions
###################
###################

# class LeastSquares:
#     """
#     Generic least-squares cost function with error.
#     """

#     errordef = Minuit.LEAST_SQUARES # for Minuit to compute errors correctly

#     def __init__(self, model, x, y, invcov):
#         self.model = model  # model predicts y for given x
#         self.func_code = make_func_code(describe(model)[1:])
#         self.x = np.asarray(x)
#         self.y = np.asarray(y)
#         self.invcov = np.asarray(invcov)

#     def __call__(self, *par):  # we accept a variable number of model parameters
#         ym = self.model(self.x, *par)
#         return (self.y - ym) @ self.invcov @ (self.y - ym)

#     @property
#     def ndata(self):
#         return len(self.x)

# def do_fit(model, xd, yd, invcov, **kwargs):
#     """
#     Do a fit over a set of data
#     """
#     lsq = LeastSquares(model, xd, yd, invcov)
#     m = Minuit(lsq, **kwargs)
#     m.migrad()
#     m.hesse()
#     # m.minos()

#     return m

# def do_fit_limits(model, xd, yd, invcov, limits, **kwargs):
#     """
#     Do a fit over a set of data, place limits on the parameters
#     """
#     lsq = LeastSquares(model, xd, yd, invcov)
#     m = Minuit(lsq, **kwargs)
#     m.limits = limits
#     m.migrad()
#     m.hesse()
#     # m.minos()

#     return m

# def do_fit_limits_fixed(model, xd, yd, invcov, limits, fixed, **kwargs):
#     """
#     Do a fit over a set of data, place limits on the parameters or fix some of them
#     """
#     lsq = LeastSquares(model, xd, yd, invcov)
#     m = Minuit(lsq, **kwargs)
#     m.limits = limits
#     m.fixed = fixed
#     m.migrad()
#     m.hesse()
#     # m.minos()

#     return m

# def do_fit_fixed(model, xd, yd, invcov, fixed, **kwargs):
#     """
#     Do a fit over a set of data, place limits on the parameters or fix some of them
#     """
#     lsq = LeastSquares(model, xd, yd, invcov)
#     m = Minuit(lsq, **kwargs)
#     m.fixed = fixed
#     m.migrad()
#     m.hesse()
#     # m.minos()

#     return m

# def fit_data_input_cov(xdata, ydata_and_m, fun, inv_cov_ydata, fitfun=dict(fitfun="do_fit"), verb = 0, **kwargs):
#     """
#     Perform a JK fit to ydata=fun(xdata)
#     covariance matrix is provided as input
#     provide the initial value of the fit parameters as kwargs
#     This function should be used internally
#     """

#     ydata, m_ydata = ydata_and_m

#     if verb > 1:
#         print("Mean of data, and size of cov matrix")
#         print(m_ydata)
#         print(np.shape(inv_cov_ydata))
#         print("-------")

#     if verb > 2:
#         print("Cov matrix")
#         mprint(np.linalg.inv(inv_cov_ydata),1e9)
#         print("-------")
#         print(np.linalg.svd(inv_cov_ydata,compute_uv=False))
#         print("-------")
#         print(1/np.linalg.svd((covmatense(ydata, m_ydata)),compute_uv=False))

    
#     if verb > 0:
#         print("Correlation matrix of the data")
#         mprint(cormatense(ydata, m_ydata))

#         print("-------")     
#         print("Correlation matrix used for the fit")
#         mprint(cov2cor(np.linalg.inv(inv_cov_ydata)))
#         print("-------")


#     # Do fit to mean to get priors
#     # Two types of fit so far, can extend it arbitrarily
#     if fitfun["fitfun"] == "do_fit" :
#         mean_fit = do_fit(fun, xdata, m_ydata, inv_cov_ydata, **kwargs)
#     elif fitfun["fitfun"] == "do_fit_limits":
#         mean_fit = do_fit_limits(fun, xdata, m_ydata, inv_cov_ydata, fitfun["limits"], **kwargs)
#     elif fitfun["fitfun"] == "do_fit_fixed":
#         mean_fit = do_fit_fixed(fun, xdata, m_ydata, inv_cov_ydata, fitfun["fixed"], **kwargs)
#     elif fitfun["fitfun"] == "do_fit_limits_fixed":
#         mean_fit = do_fit_limits_fixed(fun, xdata, m_ydata, inv_cov_ydata, fitfun["limits"], fitfun["fixed"], **kwargs)
#     else:
#         raise ValueError("Fit unsupported, options are do_fit, do_fit_limits and do_fit_limits_fixed")


#     if not mean_fit.valid:
#         print(mean_fit)
#         raise Exception("The fit to the mean data was not valid")

#     if verb > 0:
#         print("Fit to the mean result:")
#         print(mean_fit)
#         print("-------")

#     priors = mean_fit.values
#     priordict = {}
#     for nn, ele in enumerate(mean_fit.parameters):
#         priordict[ele] = priors[nn]

#     # Do fit to each Jackknife sample
#     y_down = jackdown(ydata)

#     chi2 = []
#     paramsjk = []
#     failed = 0

#     for nn, jk_y in enumerate(y_down):
#         # already checked the types in the dictionary, no need to add "try:" here
#         if fitfun["fitfun"] == "do_fit" :
#             m = do_fit(fun, xdata, jk_y, inv_cov_ydata, **priordict)
#         elif fitfun["fitfun"] == "do_fit_limits":
#             m = do_fit_limits(fun, xdata, jk_y, inv_cov_ydata, fitfun["limits"], **priordict)
#         elif fitfun["fitfun"] == "do_fit_fixed":
#             m = do_fit_fixed(fun, xdata, jk_y, inv_cov_ydata, fitfun["fixed"], **priordict)
#         elif fitfun["fitfun"] == "do_fit_limits_fixed":
#             m = do_fit_limits_fixed(fun, xdata, jk_y, inv_cov_ydata, fitfun["limits"], fitfun["fixed"], **priordict)

#         if m.valid:
#             chi2.append(m.fval)
#             paramsjk.append(m.values)
#         else:
#             failed += 1
#             print("The fit {0} did not converge".format(nn))
#             print(jk_y)
#             print(m)
#             print("-------")

#     siz = ydata.shape[0]
#     print("{0} out of {1} JK fits successful".format(siz-failed,siz))

#     params_up = jackup(paramsjk)

#     # the average is independent of whether or not the ensemble is scaled up or down
#     # chi2_up = jackup(chi2)
#     mean_chi2 = meanense(chi2)


#     return mean_fit, params_up, mean_chi2

# def invert_cov(cov_ydata, svd_reset=None):
#     """
#     svd_reset can be, either use stores both to use later:
#     * num_reset: keep the n largest singular values
#     * rat_reset: keep (correlation) singular values above a ratio wrt the largest value
#     """

#     if svd_reset == None:
    
#         # Invert the covariance matrix
#         try:
#             inv_cov_ydata=np.linalg.inv(cov_ydata)
#         except:
#             print("Eigenvalues of the covariance matrix should be non-zero, and by def positive")
#             print(np.linalg.eigvalsh(cov_ydata))
#             raise Exception('Inversion of covariance data failed, maybe try an svd reset with inversion keyword')

#     else:
#         u, s, vh = np.linalg.svd(cov_ydata, hermitian = True)

#         cor_ydata = cov2cor(cov_ydata)
#         scor = np.linalg.svd(cor_ydata, compute_uv=False, hermitian = True)

#         s_reset = np.zeros(len(s))

#         if "num_reset" in svd_reset:

#             kept_svds = len(s) - svd_reset["num_reset"]
 
#             # get the ratio reset in case we need later
#             if kept_svds == len(s):
#                 svd_reset["rat_reset"] = 0
#             else:
#                 svd_reset["rat_reset"] = scor[kept_svds]/scor[0]

#             if kept_svds < 1:
#                 raise Exception("Cannot remove {0} singular values, the dof are only {1}".format(svd_reset["num_reset"], len(s)))

#             s_reset[:kept_svds] = 1/s[:kept_svds]

#         elif "rat_reset" in svd_reset:
#             s_reset[0] = 1./s[0]
#             kept_svds = 1
#             for nn, sing_val in enumerate(scor[1:]):
#                 if sing_val/scor[0] > svd_reset["rat_reset"]:
#                     kept_svds += 1
#                     s_reset[nn+1] = 1/s[1+nn]
#                 else:
#                     break

#             svd_reset["num_reset"] = len(s) - kept_svds

#         else:
#             raise ValueError(
#             """
#             svd_reset only has:
#             * num_reset: keep the n largest singular values
#             * rat_reset: keep (correlation) singular values above a ratio wrt the largest value
#             """)

#         inv_cov_ydata = u @ np.diag(s_reset) @ vh 
#         print(svd_reset)

#     return inv_cov_ydata

# def fit_data(xdata, ydata, fun, verb = 0, limits=[None], inversion=None, **kwargs):
#     """
#     Perform a JK fit to ydata=fun(xdata) using the ydata ensemble to calculate covariance
#     provide the initial value of the fit parameters as kwargs
    
#     * Limits can be provided in a list, the list needs to be the same length as the variables
#         - None for no limit
#         - use 1 value to fix variable
#         - [low_lim, up_lim], one of them can be None if no up/low lim
#             if low_lim=up_lim the variable will also be fixed

#     *Inversion for covariance:
#         - None for normal inversion
#         - "diag" to ignore off diagonal covariance elements in inversion
#         - dict(num_reset) to reset a certain number of smallest *covariance* eigvals 
#         - dict(rat_reset) to reset *correlation* eigvals smaller than rat_reset * (bigger eigval)
#     """
#     # print(xdata,ydata)

#     # check for limits and fixed
#     lim_fix = dict(lims = [], fixes = [])
#     num_lim = 0
#     num_fixed = 0
#     for limit in limits:
#         if limit == None:
#             lim_fix['lims'].append(limit)
#             lim_fix['fixes'].append(False)

#         elif len(limit) == 1 or limit[0] == limit[1]:
#             lim_fix['lims'].append(None)
#             lim_fix['fixes'].append(True)
#             num_fixed += 1
#             num_lim += 0

#         else:
#             lim_fix['lims'].append(limit)
#             lim_fix['fixes'].append(False)
#             num_lim +=1

#     if num_lim:
#         if num_fixed:
#             fitfun = dict(fitfun="do_fit_limits_fixed",limits=lim_fix['lims'],fixed=lim_fix['fixes'])
#         else:
#             fitfun = dict(fitfun="do_fit_limits",limits=lim_fix['lims'])

#     if num_fixed:
#         fitfun = dict(fitfun="do_fit_fixed",fixed=lim_fix['fixes'])        
#     else:
#         fitfun = dict(fitfun="do_fit")


#     # the average is used for the initial fit, and to compute the covariance
#     m_ydata = meanense(ydata)

#     # check for how to invert covariance
#     if inversion == "diag":
#         cov_ydata = np.diag(np.diag(covmatense(ydata, m_ydata)))
#         inv_cov = invert_cov(cov_ydata)

#     else:
#         cov_ydata = covmatense(ydata, m_ydata)
#         inv_cov = invert_cov(cov_ydata, inversion)


#     return fit_data_input_cov(xdata, (ydata, m_ydata), fun, inv_cov, fitfun=fitfun, verb=verb, **kwargs)


def order_number(number, p = 1):
    """ Get the order of a number given p significant digits 
    If p = 1 
    Assign order 0 to numbers [0.95, 9.5)
    Assign order 1 to numbers [9.5, 95) ...
    If p = 2
    Assign order 0 to numbers [0.995, 9.95)
    Assign order 1 to numbers [9.95, 99.5) ...
    ...
    """
    number = np.abs(number)
    numord = int(np.floor(np.log10(number)))

    if np.round( number *10**(-numord + p - 1)) == 10**p:
        numord += 1

    return numord

def p_significant_digits(num, p=1):
    """
    For 1 sd [0.95, 1.5) -> 1, [1.5, 2.5) -> 2, ..., [8.5, 9.5) -> 9
    For 2 sd [0.995, 1.05) -> 1, [1.05, 1.15) -> 1.1, ..., [9.85, 9.95) -> 9.9
    For 3 sd [0.9995, 1.005) -> 1, [1.005, 1.015) -> 1.01, ..., [9.985, 9.995) -> 9.99
    Returns numbers in [1,10) range
    """
    if num == 0:
        return 0

    if p <= 0:
        raise ValueError(f"p should be positive, not {p}")

    numord = order_number(num, p)

    f = numord - (p - 1)

    # Note that round chooses the closest even number 3.5 and 4.5 -> 4
    return  np.round(num/10**f) * 10 **(1 - p)
        

def percent_significant_digits(num, percent=9):
    """ Write a number up to the significant digits
    that do not change its value more than the given percent
    return number [1,10), exponent, and number of sd """
    if np.abs(percent - 50) >=  50:
        raise ValueError(f"percent should be in (0,100), not {percent}")

    if num == 0:
        return num, 0, 1

    p = 1
    val1 = p_significant_digits(num, p)
    order1 = order_number(num, p)
    
    p += 1
    val2 = p_significant_digits(num, p)
    order2 = order_number(num, p)
    # print(val1, val2)
    change = np.abs( (val2*10**order2 - val1*10**order1)/(val1*10**order1) ) * 100

    while change > percent:
        val1 = val2
        order1 = order2

        p += 1
        val2 = p_significant_digits(num, p)
        order2 = order_number(num, p)
        # print(val1, val2)
        change = np.abs( (val2*10**order2 - val1*10**order1)/(val1*10**order1) ) * 100

    digits = p - 1

    return val1, order1, digits


def simple_order_number(number):
    """ 
    Get the order of a number: floor of the log10 of the number
    Assign order 0 to numbers [1, 10)
    Assign order 1 to numbers [10, 100)
    ...
    """
    number = np.abs(number)
    numord = int(np.floor(np.log10(number)))

    return numord

def PDG_rounding_rules(num):
    """
    Implement the rounding rules from the PDG,
    add our implementation of p_significant_digits.
    Take the first three significant digits
    * If <355, round to two significant digits.
    * If >=355 and <949, round to one significant digit.
    * If >=950, round to 1000 and keep _two_ significant digits.
    -------
    Returns
    -------
    out: [rounded number in (1,10], number of significant digits, order]
    if input is 0: [0, 1, 0]
    """
    # special case of 0
    if num == 0:
        return [0, 1, 0]

    # order = order_number(num, 3)
    # threedig = 100*p_significant_digits(num, 3)
    order = simple_order_number(num)
    threedig = num/10**(order - 2)

    if threedig < 355:
        # return [np.round(threedig/10) * 10**(order - 1), 2, order]
        return [np.round(threedig/10)/10, 2, order]
    elif threedig < 950:
        # return [np.round(threedig/100) * 10**(order), 1, order]
        return [np.round(threedig/100), 1, order]
    else:
        # return [np.round(threedig/100) * 10**(order), 2, order]
        return [np.round(threedig/100)/10, 2, order+1]


def value_error_rounding(value, error):
    """ 
    Get a value and error and return them in the rounding convention 

    Returns
    -------
    out : dict with the following elements
        vsd : value with decimal point at printing convention
        vpr : precision to write value
        esd : error with decimal point at printing convention
        epr : precision to write error
        ord : order of scientific notation ( if any )
        ntn : string pm or sh for notation to use
        sci : boolean to use or not scientific notation
    """
    
    esd, errd, eord = PDG_rounding_rules(error)
    
    # special case zero error (eg fixed params)
    if esd==0:
        out = dict(ntn = 'sh')
        vord = simple_order_number(value)
        vsd = value/10**vord
        
        # scientific notation, keep at most keep at most 3 decimals
        if vord < -3 or 2 < vord:
            vdecimals = 3
            while vdecimals:
                if np.round(vsd * 10**vdecimals) % 10 == 0:
                    vdecimals -= 1
                else:
                    break
                
            out.update(vsd = vsd, vpr = vdecimals, esd = 0, epr = 0, ord = vord, sci= True)
            return out
        
        # keep at most 3 decimals, or 4sd
        else:
            vdecimals = max(3, - vord  + 3)
            while vdecimals:
                if np.round(value * 10**vdecimals) % 10 == 0:
                    vdecimals -= 1
                else:
                    break
                    
            out.update(vsd = value, vpr = vdecimals, esd = 0, epr = 0, ord = vord, sci= False)
            return out
            
        
    # greater eord than vord
    if np.abs(value) < 9.5*10**(eord - 1):
        # in this case we do not need pm, only short hand notation
        out = dict(ntn = 'sh')
        
        # keep only 1 digit
        vsd = p_significant_digits(value, 1)
        vord = order_number(value, 1)
        
        # scientific notation
        if vord < -3 or 2 < vord:
            # adjust error to be in written wrt to value order
            epadded = esd * 10**(eord - vord)
            out.update(vsd = vsd, vpr = 0, esd = epadded, epr = 0, ord = vord, sci= True)
            return out
        
        # numbers with no decimals, have both written in their order
        elif vord >= 0:
            vpadded = vsd * 10**(vord)
            epadded = esd * 10**(eord)
            out.update(vsd = vpadded, vpr = 0, esd = epadded, epr = 0, ord = vord, sci= False)
            return out
        
        # numbers with decimals, value with decimals, error wrt to value order
        else:
            
            vpadded = vsd * 10**(vord)
            vdecimals = -vord
            epadded = esd * 10**(eord - vord)
            out.update(vsd = vpadded, vpr = vdecimals, esd = epadded, epr = 0, ord = vord, sci= False)
            return out
            
        
    # values same order as error use (in general) pm notation
    elif np.abs(value) < 9.95 * 10**(eord):
        vsd = np.round( value / 10**(eord - 1) ) / 10
        vord = eord
        vald = 2

        # get rid of unnecessary zeros
        if np.round(vsd*10) % 10 == 0:
            vald = 1
        if np.round(esd*10) % 10 == 0:
            errd = 1

        # scientific notation
        if vord < -3 or 2 < vord:
            out = dict(sci = True)

            # sh notation for *one* digit in value and error
            if vald*errd == 1:
                out.update(vsd = vsd, vpr = 0, esd = esd, epr = 0, ntn = 'sh', ord = vord)
                return out
            # pm notation otherwise
            else:
                out.update(vsd = vsd, vpr = vald-1, esd = esd, epr = errd - 1, ntn = 'pm', ord = vord)
                return out

        # NO scientific notation
        out = dict(sci = False)
        # pm notation (in general) for value order zero
        if vord == 0:
            
            # sh notation for *one* digit in value and error
            if vald*errd == 1:
                out.update(vsd = vsd, vpr = 0, esd = esd, epr = 0, ntn = 'sh', ord = vord)
                return out
            # pm notation otherwise
            else:
                out.update(vsd = vsd, vpr = vald - 1, esd = esd, epr = errd - 1, ntn = 'pm', ord = vord)
                return out

        # sh notation for all others
        # write in their order for numbers without decimals
        elif vord > 0:
            vpadded = vsd * 10**vord
            epadded = esd * 10**eord
            out.update(vsd = vpadded, vpr = 0, esd = epadded, epr = 0, ntn = 'sh', ord = vord)
            return out

        # value with decimals
        else:
            vpadded = vsd * 10**vord
            
            # can get rid of zeros for *one* digit in value and error
            if vald*errd == 1:
                vdecimals = - vord
                epadded = esd
            
            # otherwise get all decimals and pad error
            else:
                vdecimals = - vord  + 1
                epadded = esd * 10

            out.update(vsd = vpadded, vpr = vdecimals, esd = epadded, epr = 0, ntn = 'sh', ord = vord)
            return out
    
    # no need for else statement
    # values order of magnitude greater than error will (in general) use sh notation
    valround = np.round( value / 10**(eord - (errd - 1)) ) * 10**(eord - (errd - 1) )
    vord = simple_order_number(valround)
    vsd = valround/10**vord
    vald = vord - eord + errd
    

    # errors with 2 sd digits can discriminate with more precision
    if errd == 2:
        # get rid of unnecessary zeros
        # check if last digits of value and error are zeros
        unneeded_zeros = False
        if np.round(vsd * 10**(vald - 1)) % 10== 0 and np.round(esd*10) % 10 == 0:
            unneeded_zeros = True

        # scientific notation, sh notatioon
        if vord < -3 or 2 < vord:
            out = dict(sci = True, ntn = 'sh')
            vdecimals = vald - 1

            # unneeded zeros we can decrease the decimals
            if unneeded_zeros:
                vdecimals -= 1
            # otherwise we need to multiply error by 10
            else:
                esd *= 10
                
            out.update(vsd = vsd, vpr = vdecimals, esd = esd, epr = 0, ord = vord)
            return out

        # NO scientific notation
        out = dict(sci = False)
        
        # pm notation for eord == 0 and no unneeded_zeros
        if eord == 0:
            # sh for uneeded zeros case
            if unneeded_zeros:
                out.update(vsd = valround, vpr = 0, esd = esd, epr = 0, ntn ='sh', ord = vord)
                return out
            # pm when needed decimals
            else:
                out.update(vsd = valround, vpr = vald - 1, esd = esd, epr = errd - 1, ntn ='pm', ord = vord)
                return out
        
        # sh notation, and value and error in full for eord > 0 (vord = 2)
        # eord cannot be 2 or more, since that would mean vord > 2, and that case has been considered
        elif eord > 0:
            # consistency check
            if vord != 2:
                print(value, error)
                raise Exception
                
            out.update(vsd = valround, vpr = 0, esd = esd*10**eord, epr = 0, ntn ='sh', ord = vord)
            return out

        # sh notation for all other cases (eord < 0)
        else:
            vdecimals = - eord + errd - 1

            # if last of both digits are zeros we can ignore them
            if unneeded_zeros:
                vdecimals -= 1
            # otherwise we should multiply esd * 10
            else:
                esd = esd * 10

            out.update(vsd = valround, vpr = vdecimals, esd = esd, epr = 0, ntn ='sh', ord = vord)
            return out

    # errors with 1 sd digits are more loose
    elif errd == 1:
        # 1 sd digit always uses sh notation
        out = dict(ntn = 'sh')
        # scientific notation
        if vord < -3 or 2 < vord:
            out.update()
            vdecimals = vald - 1

            out.update(vsd = vsd, vpr = vdecimals, esd = esd, epr = 0, ord = vord, sci = True)
            return out
        
        # NO scientific notation
        else:
            # if eord is negative we need decimals
            vdecimals = max(0, -eord)
            
            # if eord is positive we need to write error in its order
            if eord > 0:
                esd *= 10**eord

            out.update(vsd = valround, vpr = vdecimals, esd = esd, epr = 0, ord = vord, sci = False)
            return out          



def ve_dict2string(ve_dict):
    """ 
    Change from a dictionary with
    vsd, vpr, esd, epr, ntn, sci, vord, name
    To latex printable, eg
    $num(unc)x10^n$ 
    ...
    """

    if ve_dict['sci']:
        if ve_dict['ntn'] == 'pm':
            return f"({ve_dict['vsd']:.{ve_dict['vpr']}f} \\pm {ve_dict['esd']:.{ve_dict['epr']}f})\\times 10^{{{ve_dict['ord']}}}"

        elif ve_dict['ntn'] == 'sh':
            return f"{ve_dict['vsd']:.{ve_dict['vpr']}f}({ve_dict['esd']:.{ve_dict['epr']}f})\\times 10^{{{ve_dict['ord']}}}"

        else:
            ntn = ve_dict['ntn']
            st= f"Unsupported {ntn} option, only pm and sh notation supported"
            raise ValueError(st)

    else:
        if ve_dict['ntn'] == 'pm':
            return f"{ve_dict['vsd']:.{ve_dict['vpr']}f} \\pm {ve_dict['esd']:.{ve_dict['epr']}f}"

        elif ve_dict['ntn'] == 'sh':
            return f"{ve_dict['vsd']:.{ve_dict['vpr']}f}({ve_dict['esd']:.{ve_dict['epr']}f})"

        else:
            ntn = ve_dict['ntn']
            st= f"Unsupported {ntn} option, only pm and sh notation supported"
            raise ValueError(st)



def add_fit_info_ve(p, value, error):
    """
    Get a formatted string of the form p = value(error) or p = value +/- error
    where p is the name of the variable
    """    
    ve_dict = value_error_rounding(value, error)
    ve_dict['name'] = p
    return f"${p} = " + ve_dict2string(ve_dict) + "$"

def add_fit_info_ve_old(p, value, error):
    """
    Get a formatted string of the form p = value(error) or p = value +/- error
    where p is the name of the variable
    """
    
    if np.abs(value) >= error:
        
        # round error to 9% significance
        esd, eord, errd = percent_significant_digits(error, 9)
        
        # special case of 0 error, e.g. for fixed params
        if esd == 0:
            vsd, vord, vald = percent_significant_digits(value, 9)
            valround = vsd * 10**vord
            eord = vord
        else:
            # round value to the nearest integer in (order of the error minus 1)
            valround = np.round(value / 10**(eord - 1)) * 10**(eord - 1)
            vord = int(np.floor(np.log10(valround)))
            vsd = valround/10**vord
        
        # the number of signigicant digits should be the order difference + 2
        vald = vord - eord + 2
        # except if there are zeros to the right, those should not count as vald
        for order in range(eord-1,vord):
            if np.floor(valround/10**order) % 10 == 0: 
                vald -= 1
            else:
                break    
        
        # print(vsd, vord, vald)
        # print(esd, eord, errd)
        
        # scientific notation
        if vord > 2 or -3 > vord:
             
            # p/m notation when same order and either has more than 1 digit
            if vord == eord and vald*errd != 1:
                fit_info = f"${p} = ({vsd:.{vald-1}f} \\pm {esd:.{errd-1}f})\\times 10^{{{vord}}}$"

            # shorthand notation for all other cases
            else:
                # have value with as many decimals as the maximum of val or error (wrt vord)
                decimals = max(vord + vald - vord - 1, vord + errd - eord - 1)

                # if error has more than 1 digit
                # or the value decimals exceed the error last digit position, by def at most by 1
                # multiply the esd by 10
                if errd != 1 or decimals != vord + errd - eord - 1:
                    esd = esd * 10
                fit_info = f"${p} = {vsd:.{decimals}f}({esd:.0f})\\times 10^{{{vord}}}$"

        # non-scientific notation        
        else:
            if eord == 0:
                # p/m notation when error has two digits or value has non-zero decimal
                if errd == 2 or vald > vord + 1:
                    # value is at least order 0, decimal places determined by val > vord + 1
                    vald = max(vald - vord - 1, 0)
                    fit_info = f"${p} = {valround:.{vald}f} \\pm {esd:.{errd-1}f}$"
                
                # shorthand notation for no decimals in either
                else:
                    fit_info = f"${p} = {valround:.0f}({esd:.0f})$"
                
            # error positive order, i.e. neither value or error will have decimals
            elif eord > 0:
                esd = esd * 10**(eord)
                fit_info = f"${p} = {valround:.0f}({esd:.0f})$"
                
            # error below decimal
            else:
                # have value with as many decimals as the maximum of val or error
                decimals = max(vald - vord - 1, errd - eord - 1)
                
                # if error has more than 1 digit
                # or the value decimals exceed the error last digit position, by def at most by 1
                # multiply the esd by 10
                if errd != 1 or decimals != errd - eord - 1:
                    esd = esd * 10

                fit_info = f"${p} = {valround:.{decimals}f}({esd:.0f})$"
                                
    # cases when error dominates
    else:
        # round value to 9% significance
        vsd, vord, vald =  percent_significant_digits(value, 9)
        
        # special case for value 0, eg when fixing val
        if vsd == 0:
            esd, eord, errd = percent_significant_digits(error, 9)
            errround = esd * 10**eord
            vord = eord
        else:
            # round error to nearest integer in (order of the value minus 1)    
            errround = np.round(error / 10**(vord - 1)) * 10**(vord - 1)
            eord = int(np.floor(np.log10(errround)))
            esd = errround/10**eord
                
        # if error is order of magnitude bigger round to 9% significance 
        if eord > vord:
            esd, eord, errd = percent_significant_digits(error, 9)
            
            # for much larger errors keep only 1 sd in value
            if eord - 1 > vord:
                vsd = p_significant_digits(value, 1)
                vord = order_number(value, 1)
                vald = 1           
        else:
            # for errors with same order as value remove possible trailing zero for errd
            if (errround/ 10**(vord - 1)) % 10 == 0:
                errd = 1
            else:
                errd = 2   
                       
        # print(vsd, vord, vald)
        # print(esd, eord, errd)    
        
        # numbers that need scientific notation
        if vord > 2 or vord < -3:
            # p/m notation when same order and either has more than 1 digit
            if eord == vord and vald*errd != 1:
                fit_info = f"${p} = ({vsd:.{vald-1}f} \\pm {esd:.{errd-1}f})\\times 10^{{{vord}}}$"
            
            # when error is greater, which means vald = 1
            # or (same order both val and err have 1 digit)
            else:
                esd = esd * 10 ** (eord - vord)
                fit_info = f"${p} = {vsd:.0f}({esd:.0f})\\times 10^{{{vord}}}$"
        
        # no scientific notation
        else:
            vsd = vsd * 10**vord
            
            # when value order 0
            if vord == 0:
                # p/m notation if value has two digits or error has non-zero decimal
                if vald == 2 or errd - eord > 1:
                    errd = max(errd - eord - 1, 0)
                    fit_info = f"${p} = {vsd:.{vald-1}f} \\pm {esd:.{errd}f}$"
                    
                else:
                    fit_info = f"${p} = {vsd:.0f}({esd:.0f})$"
                        
            # value above decimal, and hence error too
            if vord > 0:
                esd = esd * 10**eord
                fit_info = f"${p} = {vsd:.0f}({esd:.0f})$"
                
            # value below decimal
            if vord < 0:
                # cover either error or value total decimals
                decimals = max(vald - vord - 1, errd - eord - 1)
                
                # when error has more decimlas we need times 10
                if eord == vord and errd > vald:
                    esd *= 10
                else:
                    # otherwise we need times order difference and the value digits - 1
                    esd = esd * 10 ** (eord - vord + vald - 1)
                
                fit_info = f"${p} = {vsd:.{decimals}f}({esd:.0f})$"
                
    return fit_info

def summarize_fit_result(fit_data_result, pretty_vars = [], svd_reset=None):
    return summarize_result(fit_data_result[0], fit_data_result[1], fit_data_result[2], range(int(fit_data_result[0].ndof + fit_data_result[0].nfit)), pretty_vars, svd_reset)


def summarize_result(mean_fit, params_ense, mean_chi2, xdata, pretty_vars = [], svd_reset=None):
    mean_pars = meanense(params_ense)
    cov_pars = covmatense(params_ense,mean_pars)
    cor_pars = cormatense(params_ense,mean_pars)

    if svd_reset == None or svd_reset["num_reset"]==0:
        fit_info = [f"$\\chi^2 / n_\\mathrm{{dof}} = {mean_chi2:.1f} / ({len(xdata)} - {mean_fit.nfit}) = {mean_chi2/(len(xdata) - mean_fit.nfit):.2f}$"]
    else:
        rmv_dof = svd_reset["num_reset"]
        fit_info = [f"$\\chi^2 / n_\\mathrm{{dof}} = {mean_chi2:.1f} / ([{len(xdata)} - {rmv_dof}] - {mean_fit.nfit}) = {mean_chi2/(len(xdata) - rmv_dof - mean_fit.nfit):.2f}$"]
        

    names = mean_fit.parameters

    if len(mean_pars) == len(pretty_vars):
        names = pretty_vars

    nn = 0
    for p, v, e in zip(names, mean_pars, np.diag(cov_pars)**(1/2)):
        if mean_fit.fixed[nn]:
            e = 0
            p = p + "\\mathrm{\\ [fixed]}"
        
        
        fit_info.append(add_fit_info_ve(p, v, e))
        
        nn +=1

    return [mean_pars, np.diag(cov_pars)**(1/2)], cor_pars, fit_info


def summarize_result_corr(pars, mean_chi2, lenxdata, pretty_vars):
    mean_pars = pars[0]
    err_pars = pars[1]
    cor_pars = pars[2]

    cov_pars = np.copy(cor_pars)

    for ii, err1 in enumerate(err_pars):
        for jj, err2 in enumerate(err_pars):
            cov_pars[ii,jj] *= err1*err2

    ffit = len(mean_pars)

    fit_info = [f"$\\chi^2$ / $n_\\mathrm{{dof}}$ = {mean_chi2:.1f} / ({lenxdata} - {ffit}) = {mean_chi2/(lenxdata - ffit):.2f}"]


    names = pretty_vars

    for p, v, e in zip(names, mean_pars, np.diag(cov_pars)**(1/2)):

        fit_info.append(add_fit_info_ve(p, v, e))

    return [mean_pars, np.diag(cov_pars)**(1/2)], cor_pars, fit_info    




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


    out = ax.fill_between(xarr, mean_pred-error_pred, mean_pred+error_pred, alpha=0.5, color = 'C'+str(nn), ec =None)
    out = ax.plot(xarr, mean_pred, c = 'C'+str(nn),label=lab, lw=1, zorder=1)

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

    xarr = xd

    paramsjk = jackdown(params_ense)

    pred = []
    for paramset in paramsjk:
        pred.append([model(xval, *paramset) * factor(xval) for xval in xarr])

    predup = jackup(pred)

    mean_pred = meanense(predup)
    error_pred = errormean(predup, mean_pred)

    return mean_pred, error_pred

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

        axs[1].set_title('err/val ratio')

    if len(axs) > 2:

        corr = cormatense(ydata, m_ydata)

        cmap = mpl.cm.RdBu_r
        norm = mpl.colors.Normalize(vmin=-1, vmax=1)

        matrix_plot(axs[2], xdata, corr, cmap=cmap, norm=norm)
        
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
    A helper function to make a matrix plot

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
                 marker='s',fillstyle='none',ms=4, 
                 elinewidth=1, capsize=3.5, **kwargs)
    return out

def plot_eff_mass(ax, correl, **kwargs):

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

    for (j,i),label in np.ndenumerate(mat):
        if i < j and sym:
            continue
        if i == j and hide_diag:
            continue
        col = np.abs(label)
        # Some aesthetic, we dont need 1.00, 1 is enough, we dont need 0.7 or -0.7, .7 and -.7
        text = f"{label:.2f}".replace("0.", ".").replace("1.00","1").replace("-.00","0").replace("-0","0")
        if col < 0.7:
            ax.text(i,j,text,ha='center',va='center', c='k', **kwargs)

        else:
            ax.text(i,j,text,ha='center',va='center', c='w', **kwargs)

"""
An idea for an image of the fit summary:
fit_info = []
fit_info.extend(summ[-1])
axs[0].text(0,.5,"\n".join(fit_info), transform=ax.transAxes, bbox=dict(ec='None',fc="None"), ha='left', va='center')


An idea for the correlation plots

axs[0].axis('off')

mJK.correlation_plot(axs[1], ["$"+var+"$" for var in ppvars], summ[1])
axs[1].set_yticklabels("")

mJK.add_labels_matrix(axs[1], summ[1], size=12)
"""


