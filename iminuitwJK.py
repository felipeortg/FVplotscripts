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


def read_Jack_file(filename, r_comp=0):
    """
    Read a Jack file from Roberts format: Ncfg Nt (0:r 1:c) 0 1
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

                if comp != r_comp:
                    err = []
                    err.append("The type of file does not match the requested read")
                    err.append("File: {0}, requested {1}".format(comp,r_comp))
                    err.append("Use r_comp={0:real,1:real part of complex,2:imag part of complex")
                    # print("\n".join(err))
                    raise ValueError("\n".join(err))

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

    sig2weight = np.array([[1/(s1*s2) for s1 in sig] for s2 in sig])

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

    mean = meanense(ense)
    err = errormean(ense, mean, jackknifed)

    return np.array([mean, err])

def get_mean_error_array(list_of_JK_ensems):

    me_arr = []

    for JK in list_of_JK_ensems:

        mean = meanense(JK)

        error = errormean(JK,mean)

        me_arr.append([mean[0], error[0]])

    return np.array(me_arr)

def combine_ensemble_list(list_of_ensems):
    """ take a list of 1-dim arrays, eg from get_data, and put them into a ensemble array shape (confs, indep_index) """

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
    
    Example of an td_fun operation
    def td_fun(t):
        tmin = 3
        mass = .2612
        t0 = 10
        
        return np.exp(mass*((t+tmin)-t0))
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

def do_fit_limits(model, xd, yd, invcov, limits, **kwargs):
    lsq = LeastSquares(model, xd, yd, invcov)
    m = Minuit(lsq, **kwargs)
    m.limits = limits
    m.migrad()
    m.hesse()
    # m.minos()

    return m

def fit_data_input_cov(xdata, ydata_and_m, fun, inv_cov_ydata, fitfun=dict(fitfun="do_fit"), verb = 0, **kwargs):
    """
    Perform a JK fit to ydata=fun(xdata)
    covariance matrix is provided as input
    provide the initial value of the fit parameters as kwargs
    """

    ydata, m_ydata = ydata_and_m

    if verb > 1:
        print("Mean of data, and size of cov matrix")
        print(m_ydata)
        print(np.shape(inv_cov_ydata))
        print("-------")

    if verb > 2:
        print("Cov matrix")
        mprint(np.linalg.inv(inv_cov_ydata),1e9)
        print("-------")
        print(np.linalg.svd(inv_cov_ydata,compute_uv=False))
        print("-------")
        print(1/np.linalg.svd((covmatense(ydata, m_ydata)),compute_uv=False))

    
    if verb > 0:
        print("Correlation matrix of the data")
        mprint(cormatense(ydata, m_ydata))

        print("-------")     
        print("Correlation matrix used for the fit")
        mprint(cov2cor(np.linalg.inv(inv_cov_ydata)))
        print("-------")


    # Do fit to mean to get priors
    try:
        # Two types of fit so far, can extend it arbitrarily
        if fitfun["fitfun"] == "do_fit" :
            mean_fit = do_fit(fun, xdata, m_ydata, inv_cov_ydata, **kwargs)
        elif fitfun["fitfun"] == "do_fit_limits":
            mean_fit = do_fit_limits(fun, xdata, m_ydata, inv_cov_ydata, fitfun["limits"] **kwargs)
        else:
            raise ValueError("Fit supported for generic do_fit, and do_fit_limits")

    except Exception as e:
        print(e)
        raise Exception("Fix the above exception in fitfun to proceed to the fit")


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
        # already checked the types in the dictionary, no need to add "try:" here
        if fitfun["fitfun"] == "do_fit" :
            m = do_fit(fun, xdata, jk_y, inv_cov_ydata, **priordict)
        elif fitfun["fitfun"] == "do_fit_limits":
            m = do_fit_limits(fun, xdata, jk_y, inv_cov_ydata, fitfun["limits"] **priordict)

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

    # the average is independent of whether or not the ensemble is scaled up or down
    # chi2_up = jackup(chi2)
    mean_chi2 = meanense(chi2)


    return mean_fit, params_up, mean_chi2

def invert_cov(cov_ydata, svd_reset=None):

    if svd_reset != None:
        
        u, s, vh = np.linalg.svd(cov_ydata, hermitian = True)

        cor_ydata = cov2cor(cov_ydata)
        scor = np.linalg.svd(cor_ydata, compute_uv=False, hermitian = True)

        s_reset = np.zeros(len(s))

        if "num_reset" in svd_reset:

            kept_svds = len(s) - svd_reset["num_reset"]
 
            # get the ratio reset to use later
            if kept_svds == len(s):
                svd_reset["rat_reset"] = 0
            else:
                svd_reset["rat_reset"] = scor[kept_svds]/scor[0]

            if kept_svds < 1:
                raise Exception("Cannot remove {0} singular values, the dof are only {1}".format(svd_reset["num_reset"], len(s)))

            s_reset[:kept_svds] = 1/s[:kept_svds]

        elif "rat_reset" in svd_reset:
            s_reset[0] = [1./s[0]]
            kept_svds = 1
            for nn, sing_val in enumerate(scor[1:]):
                if sing_val/scor[0] > svd_reset["rat_reset"]:
                    kept_svds += 1
                    s_reset[nn+1] = 1/s[1+nn]
                else:
                    break

            svd_reset["num_reset"] = len(s) - kept_svds

        else:
            raise ValueError(
            """
            svd_reset only has:
            * num_reset: keep the n largest singular values
            * rat_reset: keep singular values above a ratio wrt the largest value
            """)

        inv_cov_ydata = u @ np.diag(s_reset) @ vh  

    else:
        # Invert the covariance matrix
        try:
            inv_cov_ydata=np.linalg.inv(cov_ydata)
        except:
            print("Eigenvalues of the covariance matrix should be non-zero, and by def positive")
            print(np.linalg.eigvalsh(cov_ydata))
            raise Exception('Inversion of covariance data failed')

    return inv_cov_ydata

def fit_data(xdata, ydata, fun, verb = 0, **kwargs):
    """
    Perform a JK fit to ydata=fun(xdata) using the ydata ensemble to calculate covariance
    provide the initial value of the fit parameters as kwargs
    """
    # print(xdata,ydata)
    m_ydata = meanense(ydata)

    cov_ydata = covmatense(ydata, m_ydata)

    return fit_data_input_cov(xdata, (ydata, m_ydata), fun, invert_cov(cov_ydata), verb=verb, **kwargs)

def fit_data_limits(xdata, ydata, fun, limits, verb = 0, **kwargs):
    """
    Perform a JK fit to ydata=fun(xdata) using the ydata ensemble to calculate covariance
    provide the initial value of the fit parameters as kwargs
    provide limits as a list of [low,high] in parameter order, use None for no limit
    """

    m_ydata = meanense(ydata)

    cov_ydata = covmatense(ydata, m_ydata)

    fitfun = dict(fitfun="do_fit",limits=limits)

    return fit_data_input_cov(xdata, (ydata, m_ydata), fun, invert_cov(cov_ydata), fitfun=fitfun, verb=verb, **kwargs)

def fit_data_diag_cov(xdata, ydata, fun, verb = 0, **kwargs):
    """
    Perform a JK fit to ydata=fun(xdata) using the diagonal elements of the covariance matrix
    provide the initial value of the fit parameters as kwargs
    """
    # print(xdata,ydata)

    m_ydata = meanense(ydata)

    cov_ydata = np.diag(np.diag(covmatense(ydata, m_ydata)))

    fitfun = dict(fitfun="do_fit")

    return fit_data_input_cov(xdata, (ydata, m_ydata), fun, invert_cov(cov_ydata), verb=verb, **kwargs)

def fit_data_svd(xdata, ydata, fun, svd_reset, verb = 0, **kwargs):
     """
     Perform a JK fit to ydata=fun(xdata) using svd resetting of the covariance matrix
     provide the initial value of the fit parameters as kwargs
     svd
     """

     m_ydata = meanense(ydata)
     cov_ydata = covmatense(ydata, m_ydata)


     return fit_data_input_cov(xdata, (ydata, m_ydata), fun, invert_cov(cov_ydata, svd_reset), verb=verb, **kwargs)


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
    For 1 sd [0.95, 1.05) -> 1, [1.05, 1.15) -> 1.1
    For 2 sd [0.995, 1.005) -> 1, [1.005, 1.015) -> 1.01
    Returns numbers in [1,10) range
    """

    sign = np.sign(num)
    num = np.abs(num)

    numord = order_number(num, p)

    f = numord - (p - 1)

    if p > 0:

        return sign * np.round(num/10**f) * 10 **(1 - p)
        # Note that round chooses the closest even number 3.5 and 4.5 -> 4

    else:
        raise ValueError(f"p should be positive, not {p}")

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


def add_fit_info_ve(p, v, e):
    """
    Get a formatted string of the form p = v(e) or p = v +/- e
    where p is the name of the variable
    v is value, e is error
    """

    if e >= np.abs(v):
        # get value to 9% significance
        vsd, vord, vp = percent_significant_digits(v, percent=9)

        # Use 2 significant digits for error
        e2sd = p_significant_digits(e, 2)
        errd = 2
        # could have a trailing zero in the error
        if (e2sd * 10)%10 == 0:
            errd = 1

        eord = order_number(e, 2)

        # numbers that need scientific notation
        if vord > 2 or vord < -3:
            if vord == eord:
                # separate into +/- if same order, e.g. 2.1 +/- 1.1

                fit_info = f"${p} = ({vsd:.{vp-1}f} \\pm {e2sd:.{errd-1}f})\\times 10^{{{vord}}}$"

            else:
                # with error order bigger than value we can do
                e2sd = e2sd * 10 ** (eord - vord)
                fit_info = f"${p} = {vsd:.0f}({e2sd:.0f})\\times 10^{{{vord}}}$"

        # formats without scientific notation
        else:
            # when none need a decimal point
            if vord - ( vp - 1 ) > -1 and eord - 1 > -1:
                val = vsd * 10 ** vord
                err = e2sd * 10 ** eord

                fit_info = f"${p} = {val:.0f}({err:.0f})$"

            # other option is that both are below the decimal point
            elif vord < 0 and eord < 0:
                val = vsd * 10 ** vord
                vald = - vord + (vp - 1)

                # equal order, non-zero second digit error, and zero second digit val
                # we need to ensure val has enough digits 
                if vord == eord and (e2sd * 10)%10 != 0 and (val * 10)%10 == 0:
                    vald+=1
                
                esd = e2sd * 10 ** (eord + vald)

                fit_info = f"${p} = {val:.{vald}f}({esd:.0f})$"

            # all other cases there is a decimal point we have to worry about
            else:
                val = vsd * 10 ** vord
                vald = 0

                # check if val needs the decimal point
                if vord < 1:
                    vald = - vord + (vp - 1)

                err = e2sd * 10 ** eord
                errd = 0

                # check if err needs the decimal point
                if eord == 0 and (e2sd * 10)%10 != 0:
                    errd = 1

                fit_info = f"${p} = {val:.{vald}f} \\pm {err:.{errd}f}$"

    else:
        # Check order of value to 2 precision digits
        # print(v)
        vsd = p_significant_digits(v, p=2)
        vord = order_number(v, p=2)
        vald = 2

        # Check error to 9 percent precision
        esd, eord, errd = percent_significant_digits(e, percent=9)

        # Increase the precision of v until all important digits are in: ov - oe + errd
        vp = vald
        while vord - eord + errd > vp:
            vp += 1
            vord = order_number(v, p=vp)

        vsd = p_significant_digits(v, p=vp)

        if vord < eord:
            esd = p_significant_digits(e, 2)
            eord = order_number(e, 2)
            errd = 2
            if vord < eord:
                raise Exception("There is a weird case where v:{0} > e:{1} but the orders are reversed...".format(v,e))


        # numbers that need scientific notation
        if vord > 2 or vord < -3:
            if errd == 1:
                # when only one digit for error
                if vald == 2 and vord == eord:
                    # only separate into +/- if two digits in value, same order, e.g. 2.1 +/- 1 
                    use_format = 'pm scientific'

                else:
                    # with one digit in error having (esd) works most of the time
                    use_format = 'shorthand scientific'

            else:
                # when two digits for error
                if vord == eord:
                    # separate into +/- if same order, e.g. 2.1 +/- 1.1
                    use_format = 'pm scientific'

                else:
                    # with value bigger than error we can do v.(ov-oe+1)f(esd)
                    use_format = 'shorthand scientific'

            if use_format == 'pm scientific':
                vald = vp - 1
                fit_info = f"${p} = ({vsd:.{vp-1}f} \\pm {esd:.{errd-1}f})\\times 10^{{{vord}}}$"

            elif use_format == 'shorthand scientific':
                vald = vord - eord + (errd - 1)
                esd = esd * 10 ** (errd - 1)
                fit_info = f"${p} = {vsd:.{vald}f}({esd:.0f})\\times 10^{{{vord}}}$"

            else:
                raise Exception(f"The value v:{v} was expected scientific, but with e:{e} did not match any category")

        # formats without scientific notation
        else:
            # When the two numbers have their sd above the decimal place
            if vord - ( vald - 1 ) > -1 and eord - (errd - 1) > -1:
                val = vsd * 10 ** vord
                err = esd * 10 ** eord

                fit_info = f"${p} = {val:.0f}({err:.0f})$"

            # When both numbers are below the decimal point
            elif eord == 0:
                # print('eord0')
                # Since the option when writing as integers has been taken
                val = vsd * 10 ** vord
                vald = vald - 1

                err = esd * 10 ** eord
                errd = errd - 1

                fit_info = f"${p} = {val:.{vald}f} \\pm {err:.{errd}f}$"

            else:
                # print('sh')
                # All other ones can use short hand notation with decimal point
                val = vsd * 10 ** vord

                # fix the case where .95 +/- .46 -> 1.0(5) to do 0.95(50)
                # specially since .95 +/- .45 -> .95(45) 
                if vord == eord and errd == 1:
                    errd += 1

                vald = - eord + (errd - 1)
                
                esd = esd * 10 ** (errd - 1)

                fit_info = f"${p} = {val:.{vald}f}({esd:.0f})$"

    return fit_info


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

    for p, v, e in zip(names, mean_pars, np.diag(cov_pars)**(1/2)):

        fit_info.append(add_fit_info_ve(p, v, e))

        # if np.round(e*1e3) == 0:
        #     fit_info.append(f"${p} = {v:.3f} \\pm {e:.0e}$")
        # else:
        #     fit_info.append(f"${p} = {v:.3f} \\pm {e:.3f}$")

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

        # if np.round(e*1e3) == 0:
        #     fit_info.append(f"${p} = {v:.3f} \\pm {e:.0e}$")
        # else:
        #     fit_info.append(f"${p} = {v:.3f} \\pm {e:.3f}$")

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
    out : list
        list of artists added to ax1
    """

    cfgs, npoints, xdata, ydata, xmasked, ymasked = maskdata(filename,mask )

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

    gcmap=plt.cm.get_cmap('Greys', 10)
    gind = [i/10 for i in range(1,6)]
    gind.extend(np.arange(5,0,-1)/10)
    cmapc = mpl.colors.ListedColormap([gcmap(xx/.63) for xx in gind])

    for nn, ycfg in enumerate(ydata):

        axs[0].plot(xdata, ycfg, '_', ls='', color = cmap( norm( nn + 1 )))

    xlims = axs[0].get_xlim()
    ylims = axs[0].get_ylim()


    Ny = 100

    y = np.linspace(ylims[0], ylims[1], num= Ny)
    X, Y = np.meshgrid(xdata, y)
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
        
        
        if col < 0.7:
            covax.text(i,j,f"{label:.2f}",ha='center',va='center', c='k', **kwargs)

        else:
            covax.text(i,j,f"{label:.2f}",ha='center',va='center', c='w', **kwargs)


def add_labels_matrix(ax, mat, sym=True, hide_diag=False, **kwargs):
    """
    Write values of matrix in a matshow axis
    ax: Axis where matshow was used
    mat: the values of the matrix in question
    kwargs for the text
    """

    for (j,i),label in np.ndenumerate(mat):
        if i < j and sym:
            continue
        if i == j and hide_diag:
            continue
        col = np.abs(label)
        # Some aesthetic, we dont need 1.00, 1 is enough, we dont need 0.7 or -0.7, .7 and -.7
        text = f"${label:.2f}$".replace("0.", ".").replace("1.00","1").replace("-.00","0").replace("-0","0")
        if col < 0.7:
            ax.text(i,j,text,ha='center',va='center', c='k', **kwargs)

        else:
            ax.text(i,j,text,ha='center',va='center', c='w', **kwargs)

"""
An idea for a image of the fit summary:
fit_info = []
fit_info.extend(summ[-1])
axs[0].text(0,.5,"\n".join(fit_info), transform=ax.transAxes, bbox=dict(ec='None',fc="None"), ha='left', va='center')


An idea fro the correlation plots

axs[0].axis('off')

mJK.correlation_plot(axs[1], ["$"+var+"$" for var in ppvars], summ[1])
axs[1].set_yticklabels("")

mJK.add_labels_matrix(axs[1], summ[1], size=12)
"""


