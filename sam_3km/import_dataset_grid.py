# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script import grid datasets"

import os
import netCDF4
import numpy as np
import xarray as xr

path = '/marconi/home/userexternal/mdasilva/user/mdasilva/sam_3km/postproc/'

	
def import_obs_srf(param, domain, dataset):

	arq = xr.open_dataset('{0}'.fotmat(path) + '{0}_{1}_{2}_mon_2018-2021_lonlat.nc'.format(param, domain, dataset))
	data = arq[param]
	time = data.sel(time=slice('2018-01-01','2021-12-31'))
	var = time.groupby('time.season').mean(dim='time')
	lat = var.lat
	lon = var.lon
	
	return lat, lon, var
	
		
def import_rcm_srf(param, domain, dataset):

	arq = xr.open_dataset('{0}'.fotmat(path) + '{0}_{1}_{2}_mon_2018-2021_lonlat.nc'.format(param, domain, dataset))
	data = arq[param]
	time = data.sel(time=slice('2018-01-01','2021-12-31'))
	var = time.groupby('time.season').mean(dim='time')
	lat = var.lat
	lon = var.lon
	
	return lat, lon, var
