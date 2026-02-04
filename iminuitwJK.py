#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2022-02-01 16:13:13
# @Author  : Felipe G. Ortega-Gama (felipeortegagama@gmail.com)
# @Version : 2.0
# Macros to run read/write jackknife files, perform operations and fits, make visualizations.


import numpy as np

import csv

# import subprocess
from iminuit import Minuit
from iminuit.util import describe, make_func_code

from packaging import version
from iminuit import __version__ as imversion
if version.parse(imversion) < version.parse("2.6"):
    raise Exception(f"iminuit version 2.6 or newer is required, version {imversion} found.")

from jack_utility import *

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
        # self.func_code = make_func_code(describe(model)[1:])
        self._parameters = {k: None for k in describe(model)[1:]} # new version, limits can be given here ...
        self.x = np.asarray(x)
        self.y = np.asarray(y)
        self.invcov = np.asarray(invcov)

    def __call__(self, *par):  # we accept a variable number of model parameters
        ym = self.model(self.x, *par)
        return (self.y - ym) @ self.invcov @ (self.y - ym)

    @property
    def ndata(self):
        return len(self.y)

class LeastSquares_bayesian_priors:
    """
    Generic least-squares cost function with error.
    """

    errordef = Minuit.LEAST_SQUARES # for Minuit to compute errors correctly

    def __init__(self, model, x, y, invcov, priors_mean, priors_var):
        self.model = model  # model predicts y for given x
        # self.func_code = make_func_code(describe(model)[1:])
        self._parameters = {k: None for k in describe(model)[1:]} # new version, limits can be given here ...
        self.x = np.asarray(x)
        self.y = np.asarray(y)
        self.invcov = np.asarray(invcov)

        if len(priors_mean) != len(priors_var):
            raise Exception("Mean and variance of priors should be same length")

        self.pars_bay_prior = np.asarray(priors_mean)
        self.pars_bay_prior_var = np.asarray(priors_var)

    def __call__(self, *par):  # we accept a variable number of model parameters
        ym = self.model(self.x, *par)

        priors = np.array(self.pars_bay_prior - [*par])**2/self.pars_bay_prior_var

        return (self.y - ym) @ self.invcov @ (self.y - ym) + np.sum(priors)

    @property
    def ndata(self):
        return len(self.y)

def do_fit_priors(model, xd, yd, invcov, bayesian_mean, bayesian_var, **kwargs):
    """
    Do a fit over a set of data
    """
    lsq = LeastSquares_bayesian_priors(model, xd, yd, invcov, bayesian_mean, bayesian_var)
    m = Minuit(lsq, **kwargs)
    m.migrad()
    m.hesse()
    # m.minos()

    return m

def do_fit(model, xd, yd, invcov, **kwargs):
    """
    Do a fit over a set of data
    """
    lsq = LeastSquares(model, xd, yd, invcov)
    m = Minuit(lsq, **kwargs)
    m.migrad()
    m.hesse()
    # m.minos()

    return m

def do_fit_limits(model, xd, yd, invcov, limits, **kwargs):
    """
    Do a fit over a set of data, place limits on the parameters
    """
    lsq = LeastSquares(model, xd, yd, invcov)
    m = Minuit(lsq, **kwargs)
    m.limits = limits
    m.migrad()
    m.hesse()
    # m.minos()

    return m

def do_fit_limits_fixed(model, xd, yd, invcov, limits, fixed, **kwargs):
    """
    Do a fit over a set of data, place limits on the parameters or fix some of them
    """
    lsq = LeastSquares(model, xd, yd, invcov)
    m = Minuit(lsq, **kwargs)
    m.limits = limits
    m.fixed = fixed
    m.migrad()
    m.hesse()
    # m.minos()

    return m

def do_fit_fixed(model, xd, yd, invcov, fixed, **kwargs):
    """
    Do a fit over a set of data, place limits on the parameters or fix some of them
    """
    lsq = LeastSquares(model, xd, yd, invcov)
    m = Minuit(lsq, **kwargs)
    m.fixed = fixed
    m.migrad()
    m.hesse()
    # m.minos()

    return m

def fit_data_input_cov(xdata, ydata_and_m, fun, inv_cov_ydata, fitfun=dict(fitfun="do_fit"), verb = 0, **kwargs):
    """
    Perform a JK fit to ydata=fun(xdata)
    covariance matrix is provided as input
    provide the initial value of the fit parameters as kwargs
    This function should be used internally
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
    # Two types of fit so far, can extend it arbitrarily
    if fitfun["fitfun"] == "do_fit" :
        mean_fit = do_fit(fun, xdata, m_ydata, inv_cov_ydata, **kwargs)
    elif fitfun["fitfun"] == "do_fit_limits":
        mean_fit = do_fit_limits(fun, xdata, m_ydata, inv_cov_ydata, fitfun["limits"], **kwargs)
    elif fitfun["fitfun"] == "do_fit_fixed":
        mean_fit = do_fit_fixed(fun, xdata, m_ydata, inv_cov_ydata, fitfun["fixed"], **kwargs)
    elif fitfun["fitfun"] == "do_fit_limits_fixed":
        mean_fit = do_fit_limits_fixed(fun, xdata, m_ydata, inv_cov_ydata, fitfun["limits"], fitfun["fixed"], **kwargs)
    else:
        raise ValueError("Fit unsupported, options are do_fit, do_fit_limits and do_fit_limits_fixed")


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
            m = do_fit_limits(fun, xdata, jk_y, inv_cov_ydata, fitfun["limits"], **priordict)
        elif fitfun["fitfun"] == "do_fit_fixed":
            m = do_fit_fixed(fun, xdata, jk_y, inv_cov_ydata, fitfun["fixed"], **priordict)
        elif fitfun["fitfun"] == "do_fit_limits_fixed":
            m = do_fit_limits_fixed(fun, xdata, jk_y, inv_cov_ydata, fitfun["limits"], fitfun["fixed"], **priordict)

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
    """
    svd_reset can be, either use stores both to use later:
    * num_reset: keep the n largest singular values
    * rat_reset: keep (correlation) singular values above a ratio wrt the largest value
    """

    if svd_reset == None:
    
        # Invert the covariance matrix
        try:
            inv_cov_ydata=np.linalg.inv(cov_ydata)
        except:
            print("Eigenvalues of the covariance matrix should be non-zero, and by def positive")
            print(np.linalg.eigvalsh(cov_ydata))
            raise Exception('Inversion of covariance data failed, maybe try an svd reset with inversion keyword')

    else:
        ydata_error_inv = np.diag(1/np.sqrt(np.diag(cov_ydata)))

        cor_ydata = cov2cor(cov_ydata)
        u, s, vh = np.linalg.svd(cor_ydata, hermitian = True)

        s_inv_reset = np.zeros(len(s))

        if "num_reset" in svd_reset:

            kept_svds = len(s) - svd_reset["num_reset"]
 
            # get the ratio reset in case we need later
            if kept_svds == len(s):
                svd_reset["rat_reset"] = 0
            else:
                svd_reset["rat_reset"] = s[kept_svds]/s[0]

            if kept_svds < 1:
                raise Exception("Cannot remove {0} singular values, the dof are only {1}".format(svd_reset["num_reset"], len(s)))

            s_inv_reset[:kept_svds] = 1/s[:kept_svds]

        elif "rat_reset" in svd_reset:
            s_inv_reset[0] = 1./s[0]
            kept_svds = 1
            for nn, sing_val in enumerate(s[1:]):
                if sing_val/s[0] > svd_reset["rat_reset"]:
                    kept_svds += 1
                    s_inv_reset[nn+1] = 1/s[1+nn]
                else:
                    break

            svd_reset["num_reset"] = len(s) - kept_svds

        else:
            raise ValueError(
            """
            svd_reset only has:
            * num_reset: keep the n largest singular values
            * rat_reset: keep (correlation) singular values above a ratio wrt the largest value
            """)

        inv_cov_ydata = ydata_error_inv @ u @ np.diag(s_inv_reset) @ vh @ ydata_error_inv
        print(svd_reset)

    return inv_cov_ydata

def fit_data(xdata, ydata, fun, verb = 0, limits=[None], inversion=None, **kwargs):
    """
    Perform a JK fit to ydata=fun(xdata) using the ydata ensemble to calculate covariance
    provide the initial value of the fit parameters as kwargs
    
    * Limits can be provided in a list, the list needs to be the same length as the variables
        - None for no limit
        - use 1 value to fix variable
        - [low_lim, up_lim], one of them can be None if no up/low lim
            if low_lim=up_lim the variable will also be fixed

    *Inversion for covariance:
        - None for normal inversion
        - "diag" to ignore off diagonal covariance elements in inversion
        - dict(num_reset) to reset a certain number of smallest *covariance* eigvals 
        - dict(rat_reset) to reset *correlation* eigvals smaller than rat_reset * (bigger eigval)
    """
    # print(xdata,ydata)

    # check for limits and fixed
    lim_fix = dict(lims = [], fixes = [])
    num_lim = 0
    num_fixed = 0
    for limit in limits:
        if limit == None:
            lim_fix['lims'].append(limit)
            lim_fix['fixes'].append(False)

        elif len(limit) == 1 or limit[0] == limit[1]:
            lim_fix['lims'].append(None)
            lim_fix['fixes'].append(True)
            num_fixed += 1
            num_lim += 0

        else:
            lim_fix['lims'].append(limit)
            lim_fix['fixes'].append(False)
            num_lim +=1

    if num_lim:
        if num_fixed:
            fitfun = dict(fitfun="do_fit_limits_fixed",limits=lim_fix['lims'],fixed=lim_fix['fixes'])
        else:
            fitfun = dict(fitfun="do_fit_limits",limits=lim_fix['lims'])

    if num_fixed:
        fitfun = dict(fitfun="do_fit_fixed",fixed=lim_fix['fixes'])        
    else:
        fitfun = dict(fitfun="do_fit")


    # the average is used for the initial fit, and to compute the covariance
    m_ydata = meanense(ydata)

    # check for how to invert covariance
    if inversion == "diag":
        cov_ydata = np.diag(np.diag(covmatense(ydata, m_ydata)))
        inv_cov = invert_cov(cov_ydata)

    else:
        cov_ydata = covmatense(ydata, m_ydata)
        inv_cov = invert_cov(cov_ydata, inversion)


    return fit_data_input_cov(xdata, (ydata, m_ydata), fun, inv_cov, fitfun=fitfun, verb=verb, **kwargs)


def fit_data_boot(xdata, ydata, fun, inversion=None, **kwargs):
    """
    Perform a bootstrap fit to ydata=fun(xdata) using the ydata ensemble to calculate covariance
    provide the initial value of the fit parameters as kwargs
    
    * Limits can be provided in a list, the list needs to be the same length as the variables
        - None for no limit
        - use 1 value to fix variable
        - [low_lim, up_lim], one of them can be None if no up/low lim
            if low_lim=up_lim the variable will also be fixed

    *Inversion for covariance:
        - None for normal inversion
        - "diag" to ignore off diagonal covariance elements in inversion
        - dict(num_reset) to reset a certain number of smallest *covariance* eigvals 
        - dict(rat_reset) to reset *correlation* eigvals smaller than rat_reset * (bigger eigval)
    """

    m_ydata = meanense(ydata)

    cov_ydata = np.cov(ydata, rowvar=False, ddof=1)

    # check for how to invert covariance
    if inversion == "diag":
        cov_ydata = np.diag(np.diag(cov_ydata))
        inv_cov_ydata = invert_cov(cov_ydata)

    else:
        inv_cov_ydata = invert_cov(cov_ydata, inversion)

    
    mean_fit = do_fit(fun, xdata, m_ydata, inv_cov_ydata, **kwargs)

    priors = mean_fit.values
    priordict = {}
    for nn, ele in enumerate(mean_fit.parameters):
        priordict[ele] = priors[nn] 
    
    chi2 = []
    params = []
    failed = 0

    for nn, b_y in enumerate(ydata):
        m = do_fit(fun, xdata, b_y, inv_cov_ydata, **priordict)

        if m.valid:
            chi2.append(m.fval)
            params.append(m.values)
        else:
            failed += 1
            print("The fit {0} did not converge".format(nn))
            print(b_y)
            print(m)
            print("-------")


    mean_chi2 = meanense(chi2)


    return mean_fit, params, mean_chi2


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
        if value == 0: # just in case
            out.update(vsd = value, vpr = 0, esd = 0, epr = 0, sci= False)
            return out
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
        
        # sh, keep at most 3 decimals, or 4sd
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

def value_error_rounding_str(value, error):
    ve_dict = value_error_rounding(value, error)
    return ve_dict2string(ve_dict)

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

def boot_summarize_result(fit_data_result, pretty_vars = [],  svd_reset=None):
    mean_fit = fit_data_result[0]
    params = fit_data_result[1]
    mean_chi2 = fit_data_result[2]
    xdata = range(int(fit_data_result[0].ndof + fit_data_result[0].nfit))


    mean_pars = meanense(params)
    cov_pars = np.cov(params, rowvar=False, ddof=1)
    
    if len(mean_pars) == 1:
        cov_pars = [[cov_pars]]
    
    cor_pars = cov2cor(cov_pars)

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

    names = pretty_vars
    fit_info = ['']
    n_fix = 0

    for p, v, e in zip(names, mean_pars, np.diag(cov_pars)**(1/2)):
        if e == 0:
            n_fix += 1
            p = p + "\\mathrm{\\ [fixed]}"

        fit_info.append(add_fit_info_ve(p, v, e))

    # print(n_fix)

    fit_info[0] = f"$\\chi^2/\\text{{dof}} = {mean_chi2:.2f} / ({lenxdata} - {ffit-n_fix}) = {mean_chi2/(lenxdata - ffit + n_fix):.2f}$"




    return [mean_pars, np.diag(cov_pars)**(1/2)], cor_pars, fit_info    

def remove_fixed_variables_cor(input_cor, input_names):
    fix_ind = []
    new_names = []
    for nn, col in enumerate(input_cor):
        if sum(np.abs(col)) == 1:
            fix_ind.append(nn)
        else:
            new_names.append(input_names[nn])

    new_len = len(input_cor) - len(fix_ind)

    new_cor = np.ones((new_len, new_len))

    nn_new = 0
    for nn in range(len(input_cor)):
        mm_new = nn_new + 1
        if nn in fix_ind:
            continue
        for mm in range(nn+1, len(input_cor)):
            if mm in fix_ind:
                continue
            new_cor[nn_new, mm_new] = input_cor[nn, mm]
            new_cor[mm_new, nn_new] = input_cor[mm, nn]
            mm_new += 1

        nn_new += 1

    return new_cor, new_names


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

""" Example of a fit

def exp(t, m, C):
    return np.exp(-t*m) * C

ncfgs, TL, xd, yd_check,a,b = mJK.maskdata("check_ortho_state0_state0.jack", mask=range(9))

fit_check_corr = mJK.fit_data(xd, yd_check, exp, verb=0, m=.045, C=1)
fit_check_corr[0]

sum_check = mJK.summarize_fit_result(fit_check_corr)

mass_exp = mJK.calc(fit_check_corr[1])[0,0]

def lead_exp(t):    
    return np.exp(mass_exp*(t))

yd_check_res = mJK.td_ensemble_op(lead_exp, lambda x,y: x*y, xd, yd_check)

yd_check_res_me = mJK.calc(yd_check_res)


%matplotlib inline
f, ax = plt.subplots(num=11, layout='constrained')

mJK.plot_data(ax, xd, yd_check_res_me[:,0], yd_check_res_me[:,1], 0, label='data')

mJK.plot_line_model(ax, 'fit', 1, xd, exp, fit_check_corr[1], lead_exp)


ax.text(0.05,.05, "\n".join(sum_check[2]), transform=ax.transAxes,  bbox=dict(fc="w"))
ax.text(0.05,.05, "\n".join(sum_check[2]), transform=ax.transAxes,  bbox=dict(alpha=0))
ax.legend()

"""

