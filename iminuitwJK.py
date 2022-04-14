#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2022-02-01 16:13:13
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 1.0
# Template to run a jackknife fit


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import csv
# import subprocess
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
                try:
                    cfgs, tl, comp = [int(row[0]), int(row[1]), int(row[2])] 
                    print("There are {0} configs, with {1} timeslices each".format(cfgs,tl))
                    if comp != 0:
                        print("Not prepared for {0}".format(row))
                        raise ValueError
                except Exception as e:
                    print("!!!  Error while reading header !!! \n{0}".format(row))
                    print("Format should be: Ncfg Nt (0:r 1:c) 0 1")
                    print("This reader only accepts real option")
                    raise Exception

            else:
                if count == 0:
                    thiscfgdata = []
                thiscfgdata.append([int(row[-2]),float(row[-1])])
                count += 1

                if count == tl:
                    dataarray.append(thiscfgdata)
                    count -= tl


    return cfgs, tl, dataarray

def prune_Jack_file(filename, pruned_cfgs):
    """
    Read a Jack file from Roberts format: Ncfg Nt (0:r 1:c) 0 1
    Write it again without the pruned cfgs
    """
    if filename[-5:] == '.jack':
        newfile = filename[:-5] + '_pruned.jack'
    else:
        newfile + '.pruned'

    cfgs, tl, dataarray = read_Jack_file(filename)

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
    always return the upscale ensemble, but can receive a downscaled ensemble """

    if not jackknifed:
        ensedown = jackdown(ensemble)
    else:
        ensedown = ensemble

    td_factor = [td_fun(t) for t in xd]

    cfgs = np.shape(ensedown)[0]
    
    factorjk = np.array([td_factor for i in range(cfgs)])

    op_ense = operation(factorjk, ensedown)

    return jackup(op_ense)
        

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


def order_number(number, p = 1):
    """ Get the order of a number given p significant digits """
    # If p = 1 
    # Assign order 0 to numbers [0.95, 9.5)
    # Assign order 1 to numbers [9.5, 95) ...
    # If p = 2
    # Assign order 0 to numbers [0.995, 9.95)
    # Assign order 1 to numbers [9.95, 99.5) ...
    # ...

    numord = int(np.floor(np.log10(number)))
    if np.round( number *10**(-numord + p - 1)) == 10**p:
        numord += 1

    return numord

def p_significant_digits(num, p=1):
    # For 1 sd [0.95, 1.05) -> 1, [1.05, 1.15) -> 1.1
    # For 2 sd [0.995, 1.005) -> 1, [1.005, 1.015) -> 1.01
    # Exits numbers in [1,10) range

    numord = order_number(num, p)

    f = numord - (p - 1)

    if p > 0:

        return np.round(num/10**f) * 10 **(1 - p)
        # Note that round chooses the closest even number 3.5 and 4.5 -> 4

    else:
        raise ValueError(f"p should be positive, not {p}")

def percent_significant_digits(num, percent=10):
    """ Write a number up to the significant digits
    that do not change its value more than the given percent
    return number [1,10) and exponent """
    if np.abs(percent - 50) >=  50:
        raise ValueError(f"percent should be in (0,100), not {percent}")

    if num == 0:
        return num, 0

    change = 100
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


def add_fit_info_ve(p, v, e):
    if e >= v:
        # Only use 1 significant digit if error is bigger
        v1sd = p_significant_digits(v, 1)
        # Use 2 significant digits
        e2sd = p_significant_digits(e, 2)

        vord = order_number(v, 1)
        eord = order_number(e, 2)

        if 3 > vord and vord > -1 :

            einvord = e2sd * 10 ** (eord)
            v1in0 = v1sd * 10 ** (vord)
            fit_info = f"${p} = {v1in0:.3g}({einvord:.3g})$"

        elif vord > -5 :

            if vord == eord:
                # Use 1 significant digits in this case
                e1sd = p_significant_digits(e, 1)
                eord = order_number(e, 1)

                einvord = e1sd * 10 ** (eord - vord)
                v1in0 = v1sd * 10 ** (vord)
                fit_info = f"${p} = {v1in0:.3g}({einvord:.3g})$"

            else:
                einvord = e2sd * 10 ** (eord - vord )
                v1in0 = v1sd * 10 ** (vord)
                fit_info = f"${p} = {v1in0:.3g}({einvord:.3g})$"

        else:
            einvord = e2sd * 10 ** (eord - vord)
            fit_info = f"${p} = {v1sd:.0f}({einvord:.1f})\\times 10^{vord}$"

    else:
        # Use as many significant digits for value as error 9p precision
        esd, eord, errd = percent_significant_digits(e, percent=9)

        valpro0 = eord - (errd - 1)

        vord = order_number(v, 1)

        dd = 0

        while vord - dd != valpro0:
            print(valpro0, vord, dd)
            dd += 1
            vord = order_number(v, dd + 1)

        valsd = p_significant_digits(v, dd + 1)

        if 3 > vord and vord > -1 :

            einvord = esd * 10 ** (eord)
            v1in0 = valsd * 10 ** (vord)
            fit_info = f"${p} = {v1in0:.3g}({einvord:.3g})$"

        elif vord > -5 :

            if vord == eord:
                # Use 1 significant digits in this case
                v1sd = p_significant_digits(v, 1)
                e1sd = p_significant_digits(e, 1)
                eord = order_number(e, 1)

                einvord = e1sd * 10 ** (eord - vord)
                v1in0 = v1sd * 10 ** (vord)
                fit_info = f"${p} = {v1in0:.3g}({einvord:.3g})$"

            else:
                einvord = esd * 10 ** (eord - vord )
                v1in0 = valsd * 10 ** (vord)
                fit_info = f"${p} = {v1in0:.3g}({einvord:.3g})$"

        else:
            einvord = esd * 10 ** (eord - vord)
            fit_info = f"${p} = {valsd:.{dd}f}({einvord:.1f})\\times 10^{vord}$"

    return fit_info


def summarize_result(mean_fit, params_ense, mean_chi2, xdata):
    mean_pars = meanense(params_ense)
    cov_pars = covmatense(params_ense,mean_pars)
    cor_pars = cormatense(params_ense,mean_pars)


    fit_info = [f"$\\chi^2$ / $n_\\mathrm{{dof}}$ = {mean_chi2:.1f} / ({len(xdata)} - {mean_fit.nfit}) = {mean_chi2/(len(xdata) - mean_fit.nfit):.2f}"]

    for p, v, e in zip(mean_fit.parameters, mean_pars, np.diag(cov_pars)**(1/2)):

        #fit_info.append(add_fit_info_ve(fit_info, p, v, e))

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

        axs[1].plot(xmasked, np.abs(e_ymasked/m_ymasked), 'x', c = 'C'+str(nn), alpha=0.5)

        axs[1].set_xlim(xlims)
        axs[1].set_ylim(ylims)

        axs[1].set_title('err/val ratio')

    if len(axs) > 2:

        cov = covmatense(ydata, m_ydata)
        corr = cormatense(ydata, m_ydata)

        cmap = mpl.cm.RdBu_r
        norm = mpl.colors.Normalize(vmin=-1, vmax=1)
        mat1 = axs[2].matshow(corr, cmap = cmap, norm = norm)
        
        axs[2].set_title('Correlation matrix')
        divider = make_axes_locatable(axs[2])
        cax = divider.append_axes("right", size="5%", pad=0.1)
        plt.colorbar(mat1, cax=cax)

        # Assume xdata is nicely ordered
        if len(xdata) < 6:
            dist = 1
        else:
            dist = np.floor(len(xdata)/6) 
        ticks = range(0, len(xdata), int(dist))
        tl = [str(xdata[t]) for t in ticks]

        axs[2].set_xticks(ticks)
        axs[2].set_xticklabels(tl)
        axs[2].set_yticks(ticks)
        axs[2].set_yticklabels(tl)


        if len(axs) > 3:

            mat2 = axs[3].matshow(np.log(np.abs(cov)))
            axs[3].set_title('Log |Covariance matrix|')
            divider = make_axes_locatable(axs[3])
            cax = divider.append_axes("right", size="5%", pad=0.1)
            plt.colorbar(mat2, cax=cax)


            axs[3].set_xticks(ticks)
            axs[3].set_xticklabels(tl)
            axs[3].set_yticks(ticks)
            axs[3].set_yticklabels(tl)



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

def plot_cfg_fromfile(filename, axs, nn=0, mask = [], plots = ['all']):


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


    for nn, ycfg in enumerate(ydata):

        axs[0].plot(xdata, ycfg, '_', ls='', color = cmap( norm( nn + 1 )))


    axs[0].fill_between(xdata, m_ydata + 1*np.sqrt(cfgs) * e_ydata,
             m_ydata - 1*np.sqrt(cfgs) * e_ydata,
              color = 'gray',
              alpha=1, zorder=1.7)        

    axs[0].fill_between(xdata, m_ydata + 3*np.sqrt(cfgs) * e_ydata,
             m_ydata - 3*np.sqrt(cfgs) * e_ydata,
              color = 'gray',
              alpha=0.3, zorder=1.6)        

    axs[0].fill_between(xdata, m_ydata + 5*np.sqrt(cfgs) * e_ydata,
             m_ydata - 5*np.sqrt(cfgs) * e_ydata,
              color = 'gray',
              alpha=0.1, zorder=1.5)

    axs[0].fill_between(xdata, m_ydata + 7*np.sqrt(cfgs) * e_ydata,
             m_ydata - 7*np.sqrt(cfgs) * e_ydata,
              color = 'gray',
              alpha=0.1, zorder=1.4)

    axs[0].fill_between(xdata, m_ydata + 9*np.sqrt(cfgs) * e_ydata,
             m_ydata - 9*np.sqrt(cfgs) * e_ydata,
              color = 'gray',
              alpha=0.1, zorder=1.3)

    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("bottom", size="3%", pad=0.3)

    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                 cax=cax, label='Configuration', orientation = 'horizontal')

    # cax = divider.append_axes("bottom", size="3%", pad=0.3)

    # plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
    #              cax=cax, label='Configuration', orientation = 'horizontal')

    return axs[0]
