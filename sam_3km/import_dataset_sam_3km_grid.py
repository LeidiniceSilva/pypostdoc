# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script import grid datasets"

import os
import netCDF4
import numpy as np

path = '/marconi/home/userexternal/mdasilva/user/mdasilva/sam_3km/postproc/'

	
def import_obs(param, dataset):

	arq = xr.open_dataset('{0}'.fotmat(path) + '{0}_SAM-3km_{1}_mon_2018-2021_lonlat.nc'.format(param, dataset))
	data = arq[param]
	time = data.sel(time=slice('2018-01-01','2021-12-31'))
	var = time.groupby('time.season').mean(dim='time')
	lat = var.lat
	lon = var.lon
	
	return var, lat, lon
	
		
def import_rcm(param, dataset):

	arq = xr.open_dataset('{0}'.fotmat(path) + '{0}_SAM-3km_{1}_mon_2018-2021_lonlat.nc'.format(param, dataset))
	data = arq[param]
	time = data.sel(time=slice('2018-01-01','2021-12-31'))
	var = time.groupby('time.season').mean(dim='time')
	lat = var.lat
	lon = var.lon
	
	return var, lat, lon
