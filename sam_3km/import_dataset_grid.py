# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script import grid datasets"

import xarray as xr

path = '/marconi/home/userexternal/mdasilva/user/mdasilva/sam_3km/postproc/'

	
def import_obs_srf(param, domain, dataset):

	arq = xr.open_dataset('{0}'.format(path) + '{0}_{1}_{2}_mon_2018-2021_lonlat.nc'.format(param, domain, dataset))
	data = arq[param]
	time = data.sel(time=slice('2018-01-01','2021-12-31'))
	mean = time.groupby('time.season').mean(dim='time')
	lat = mean.lat
	lon = mean.lon
	var = mean.values

	return lat, lon, var
	
		
def import_rcm_srf(param, domain, dataset):

	arq = xr.open_dataset('{0}'.format(path) + '{0}_{1}_{2}_mon_2018-2021_lonlat.nc'.format(param, domain, dataset))
	data = arq[param]
	time = data.sel(time=slice('2018-01-01','2021-12-31'))
	mean = time.groupby('time.season').mean(dim='time')
	lat = mean.lat
	lon = mean.lon
	var = mean.values
	
	return lat, lon, var
