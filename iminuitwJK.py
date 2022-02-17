#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2022-02-01 16:13:13
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 1.0
# Template to run a jackknife fit


import numpy as np
import matplotlib.pyplot as plt
import csv
import subprocess
from iminuit import Minuit
from iminuit.util import describe, make_func_code


def read_Jack_file(filename):
    """
    Read a Jack file from Roberts format: Ncfg Nt (0:r 1:c) 0 1
    """

    with open(filename, 'r') as f:
        # Fill this list with all our levels
        data = csv.reader(f, delimiter=' ') #change from default comma

        dataarray = []
        count = 0
        # thiscfgdata = []

        for nn, row in enumerate(data):
            # print(row)
            if nn == 0:
                cfgs, tl, comp = [int(row[0]), int(row[1]), int(row[4])] 
                print("There are {0} configs, with {1} timeslices each".format(cfgs,tl))
                if comp != 1:
                    print("Not prepared for {0}".format(comp))
                    raise ValueError
            else:
                if count == 0:
                    thiscfgdata = []
                thiscfgdata.append([int(row[-2]),float(row[-1])])
                count += 1

                if count == tl:
                    dataarray.append(thiscfgdata)
                    count -= tl


    return cfgs, tl, dataarray


def get_mask_from_noise(filename,nrcutoff):

    cfgs, npoints, rawdata = read_Jack_file(filename)

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


def maskdata(filename,mask=[]):
    cfgs, npoints, rawdata = read_Jack_file(filename)

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

    covmat = covmatense(ense, meanense, jackknifed)

    sig = np.sqrt(np.diag(covmat))

    sig2weight = np.array([[1/(s1*s2) for s1 in sig] for s2 in sig])

    return covmat * sig2weight


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

def ensemble_op(fun, ensemble, jackknifed = False):
    """Do an operation over an ensemble with jackknife statistics
    always return the upscale ensemble, but can receive a down scaled ensemble """

    if not jackknifed:
        ensedown = jackdown(ensemble)
    else:
        ensedown = ensemble

    enseres = fun(ensedown)

    return jackup(enseres)

# Example of an ensemble operation
# def expmult(jkarray):
#     npoints = np.shape(jkarray)[1]
    
#     factor = [np.exp(.2612*((t+tmin)-10)) for t in range(npoints)]
    
#     cfgs = np.shape(jkarray)[0]
    
#     factorjk = np.array([factor for i in range(cfgs)])
    
#     return factorjk*jkarray


def mprint(mat,round=100):
    """Print a matrix rounding to 1/round parts"""
    print(np.round(round*mat)/round)

###################
###################
# Fit functions
###################
###################

class LeastSquares:
    """
    Generic least-squares cost function with error.
    """

    errordef = Minuit.LEAST_SQUARES # for Minuit to compute errors correctly

    def __init__(self, model, x, y, invcov):
        self.model = model  # model predicts y for given x
        self.func_code = make_func_code(describe(model)[1:])
        self.x = np.asarray(x)
        self.y = np.asarray(y)
        self.invcov = np.asarray(invcov)

    def __call__(self, *par):  # we accept a variable number of model parameters
        ym = self.model(self.x, *par)
        return (self.y - ym) @ self.invcov @ (self.y - ym)

def do_fit(model, xd, yd, invcov, **kwargs):
    lsq = LeastSquares(model, xd, yd, invcov)
    m = Minuit(lsq, **kwargs)
    m.migrad()
    m.hesse()
    # m.minos()

    return m

def fit_data(xdata, ydata, fun, verb = 0, **kwargs):

    # print(xdata,ydata)

    m_ydata = meanense(ydata)

    cov_ydata = covmatense(ydata, m_ydata)

    if verb > 1:
        print("Mean of data, and size of cov matrix")
        print(m_ydata)
        print(np.shape(cov_ydata))
        print("-------")

    if verb > 2:
        print("Cov matrix")
        mprint(cov_ydata,1e9)
        print("-------")

    
    if verb > 0:
        print("Correlation matrix of the data")
        mprint(cormatense(ydata, m_ydata))
        print("-------")

    # Invert the covariance matrix
    try:
        inv_cov_ydata=np.linalg.inv(cov_ydata)
    except:
        print("Eigenvalues of the covariance matrix should be non-zero, and by def positive")
        print(np.linalg.eigvalsh(cov_ydata))
        raise Exception('Inversion of covariance data failed')

    # Do fit to mean to get priors

    mean_fit = do_fit(fun, xdata, m_ydata, inv_cov_ydata, **kwargs)

    if not mean_fit.valid:
        print(mean_fit)
        raise Exception("The fit to the mean data was not valid")

    if verb > 0:
        print("Fit to the mean result:")
        print(mean_fit)
        print("-------")

    priors = mean_fit.values
    priordict = {}
    for nn, ele in enumerate(mean_fit.parameters):
        priordict[ele] = priors[nn]

    # Do fit to each Jackknife sample
    y_down = jackdown(ydata)

    chi2 = []
    paramsjk = []
    failed = 0

    for nn, jk_y in enumerate(y_down):
        m = do_fit(fun, xdata, jk_y, inv_cov_ydata, **priordict)

        if m.valid:
            chi2.append(m.fval)
            paramsjk.append(m.values)
        else:
            failed += 1
            print("The fit {0} did not converge".format(nn))
            print(jk_y)
            print(m)
            print("-------")

    siz = ydata.shape[0]
    print("{0} out of {1} JK fits successful".format(siz-failed,siz))

    params_up = jackup(paramsjk)


    chi2_up = jackup(chi2)
    mean_chi2 = meanense(chi2_up)


    return mean_fit, params_up, mean_chi2


def summarize_result(mean_fit, params_ense, mean_chi2, xdata):
    mean_pars = meanense(params_ense)
    cov_pars = covmatense(params_ense,mean_pars)
    cor_pars = cormatense(params_ense,mean_pars)


    fit_info = [f"$\\chi^2$ / $n_\\mathrm{{dof}}$ = {mean_chi2:.1f} / ({len(xdata)} - {mean_fit.nfit}) = {mean_chi2/(len(xdata) - mean_fit.nfit):.2f}"]

    for p, v, e in zip(mean_fit.parameters, mean_pars, np.diag(cov_pars)**(1/2)):
        if np.round(e*1e3) == 0:
            fit_info.append(f"${p} = {v:.3f} \\pm {e:.0e}$")
        else:
            fit_info.append(f"${p} = {v:.3f} \\pm {e:.3f}$")

    return [mean_pars, np.diag(cov_pars)**(1/2)], cor_pars, fit_info



###################
###################
# Plot functions
###################
###################

def plot_line_model(ax, lab, nn, xd, model, params_ense, factor):

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

def print_line_model(xd, model, params_ense, factor):

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

def plot_fromfile(filename, ax1, ax2, nn=0, mask = []):

    cfgs, npoints, rawdata = read_Jack_file(filename)

    # Assume x data is fixed
    print("Assuming the xdata does not have configuration variation")
    xdata = [rawdata[0][tt][0] for tt in range(npoints)]

    ydata = np.array(rawdata)[:,:,1]

    if len(mask) != 0: # masking of the data
        for elem in mask:
            torem = xdata.index(elem)
            xdata.pop(torem)
            newydata = np.delete(ydata,(torem), axis=1)
            ydata = newydata

    print(len(xdata),np.shape(ydata))

    m_ydata = meanense(ydata)

    e_ydata = errormean(ydata, m_ydata)

    ax1.set_title('Mean Data')

    ax2.plot(xdata, e_ydata/m_ydata, 'x')
    ax2.axhline(0,lw=1,color='k')
    ax2.axhline(0.12,lw=1,color='k')

    ax2.set_title('err/val ratio')

    cov = covmatense(ydata, m_ydata)
    corr = cormatense(ydata, m_ydata)

    plt.matshow(np.log(np.abs(cov)))
    plt.title('Log |Covariance matrix|')

    plt.matshow(corr)
    plt.title('Correlation matrix')

    return plot_data(ax1, xdata, m_ydata, e_ydata, nn)


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

    out = ax.errorbar(xd, yd, yerr,
                 ls='', c = 'C'+str(nn),
                 marker='s',fillstyle='none',ms=4, 
                 elinewidth=1, capsize=3.5, **kwargs)
    return out
