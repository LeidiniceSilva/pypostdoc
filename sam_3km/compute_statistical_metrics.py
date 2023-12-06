# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "Statistical metrics for model evaluation"

import numpy as np
import scipy.stats as st

from scipy.stats import norm


def compute_corr(model, obs):

    """
    The input arrays must have the same dimentions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Pearson Linear Correlation
    """
   
    corr = np.corrcoef(obs, model)[0][1]
    
    return corr


def compute_r2(model, obs):

	"""
	The input arrays must have the same dimensions
	:Param model: Numpy array with model data
	:Param obs: Numpy array with obs data
	:Return: R-squared
	"""

	corr = np.corrcoef(obs, model)[0][1]
	r2 = corr ** 2

	return r2

    
def compute_mae(model, obs):

    """
    The input arrays must have the same dimensions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Mean Absoluty Error
    """

    mae = np.mean(np.abs(np.array(model) - np.array(obs)))
    
    return mae
    

def compute_rmse(model, obs):

    """
    The input arrays must have the same dimensions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Root Mean Square Error
    """

    rmse = np.sqrt(((np.array(model) - np.array(obs)) ** 2).mean()) 
    
    return rmse
    
     
def compute_mbe(model, obs):

    """
    The input arrays must have the same dimensions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Mean Bias Error
    """

    mbe = np.nanmean(np.array(model) - np.array(obs))
    
    return mbe


def compute_pbias(model, obs):

    """
    The input arrays must have the same dimensions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Percentage Bias
    """

    pbias = 100.0 * sum(np.array(model) - np.array(obs)) / sum(np.array(obs))
    
    return pbias
        
    
def compute_apb(model, obs):

    """
    The input arrays must have the same dimensions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Absolute Percent Bias
    """

    apb = 100.0 * sum(np.abs(np.array(model), np.array(obs))) / sum(np.array(obs))
    
    return apb


def compute_nse(model, obs):

    """
    The input arrays must have the same dimensions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Nash–Sutcliffe Efficient Coefficient
    """
    
    p1 = sum((model - obs) ** 2)
    p2 = sum((obs - np.mean(obs)) ** 2)
    nse = 1 - p1 / p2
    
    return nse
    
        
def compute_diso(model, obs):

	"""
	The input arrays must have the same dimensions
	:Param model: Numpy array with model data
	:Param obs: Numpy array with obs data
	:Return: Distance between Indices of Simulation and Observation
	"""

	p0 = np.abs(np.array(model).mean())
	p1 = np.corrcoef(model, obs)[0][1]
	p2 = np.mean(np.abs(np.array(model) - np.array(obs))) / p0
	p3 = np.sqrt(((np.array(model) - np.array(obs)) ** 2).mean()) / p0
	diso = np.sqrt((p1)** 2 + (p2)** 2 + (p3)** 2)

	return diso
    
    
def compute_ioa(model, obs):

    """
    The input arrays must have the same dimensions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Index of Agreement
    """

    p1 = (model - obs)**2
    p2 = np.abs(model - np.mean(obs))
    p3 = np.abs(obs - np.mean(obs))
    ioa = 1 - sum(p1) / sum((p2 + p3)**2)
    
    return ioa
    

def compute_av(gcm, rcm, obs):

    """
    The input arrays must have the same dimensions
    :Param rcm: Numpy array with regional model data
    :Param gcm: Numpy array with global model data
    :Param obs: Numpy array with obs data
    :Return: Added Value Index
    """

    p1 = (gcm - obs)**2
    p2 = (rcm - obs)**2
    p3 = p1 - p2
    p4 = np.max([p1, p2], axis=0)  
    av = p3 / p4
   
    return av
    
    
def compute_ivs(obs, model):

    """
    The input arrays must have the same dimensions
    :Param rcm: Numpy array with regional model data
    :Param gcm: Numpy array with global model data
    :Param obs: Numpy array with obs data
    :Return: Interannual Variability Skill Score
    """

    p1 = np.std(obs)
    p2 = np.std(model)
    p3 = p2 / p1
    p4 = p1 / p2
    ivs = (p3 - p4)**2  
    
    return ivs    
    
    
def compute_cdf(data):

	"""
	The input arrays must have the same dimensions
	:Param data: Numpy array with model or obs data
	:Return: Cumulative Density Function
	"""

	x = np.linspace(np.min(data), np.max(data))
	y = np.nanmean(x)
	z = np.nanstd(x)
	cdf = norm.cdf(x,y,z)

	return x, cdf


def compute_pdf(data):

	"""
	The input arrays must have the same dimensions
	:Param data: Numpy array with model or obs data
	:Return: Cumulative Density Function
	"""

	x = np.linspace(np.min(data), np.max(data))
	y = np.nanmean(x)
	z = np.nanstd(x)
	pdf = norm.pdf(x,y,z)

	return x, pdf
	

def compute_relative_change(rcp, hist):

    """
    The input arrays must have the same dimensions
    :Param rcp: Numpy array with rcp period model
    :Param hist: Numpy array with hist period model
    :Return: Relative change
    """

    p1 = rcp 
    p2 = hist
    p3 = p1 - p2
    p4 = p3 / p2
    rc = p4 * 100
   
    return rc
    
    	
def compute_anomaly(model, fcst):

    """
    The input arrays must have the same dimensions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Anomaly and Standard Anomaly
    """

    clim_mean = np.nanmean(model, axis=0)
    clim_std = np.nanstd(model, axis=0)
    anomaly = fcst - clim_mean
    standard_anomaly = (fcst - clim_mean)/clim_std
    
    return anomaly, standard_anomaly
    

def compute_fcst_correct(model, obs, fcst):

    """
    The input arrays must have the same dimensions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Forecast Data Correction
    """

    sim = np.sort(model)
    alpha_mod, loc_mod, beta_mod = ss.gamma.fit(sim, loc=0)
    obs = np.sort(obs)
    alpha_obs, loc_obs, beta_obs = ss.gamma.fit(obs, loc=0)

    fcst_fcst_correc = []
    for i in fcst:
        prob = ss.gamma.cdf(i, alpha_mod, scale=beta_mod)
        fcst_correc.append(ss.gamma.ppf(prob, alpha_obs, scale=beta_obs))
        
    return fcst_correct
   

def compute_kge(obs, model):

	"""
	The input arrays must have the same dimensions
	Param model: Numpy array with model data
	Param obs: Numpy array with obs data
	Return: Kling-Gupta Efficiency
	"""

	p1 = np.corrcoef(obs, model)[0][1]
	p2 = np.nanmean(model) / np.nanmean(obs)
	p3 = np.nanstd(model, ddof=0) / np.nanmean(model) * 100 
	p4 = np.nanstd(obs, ddof=0) / np.nanmean(obs) * 100 
	p5 = p3/p4
	kge = 1 - np.sqrt((1 - p1)**2 + (1 - p2)**2 + (1 - p5)**2)

	return kge
