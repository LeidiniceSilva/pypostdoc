# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script import in-situ datasets"

import xarray as xr

from dict_inmet_stations import inmet

path = '/marconi/home/userexternal/mdasilva'

skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,7793,100,105,106,107,112,117,124,135,137,139,
149,152,155,158,168,174,177,183,186,199,204,210,212,224,226,239,240,248,249,253,254,276,277,280,293,298,
303,305,306,308,319,334,335,341,343,359,362,364,384,393,396,398,399,400,402,413,416,417,422,423,426,427,
443,444,446,451,453,457,458,467,474,479,483,488,489,490,495,505,509,513,514,516,529,534,544,559,566]
	
	
def import_obs_situ(param):
	
	iy, ix, mean_i = [], [], []
	for station in range(1, 567):
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		iy.append(inmet[station][2])
		ix.append(inmet[station][3])

		arq = xr.open_dataset('{0}/OBS/BDMET/database/nc/hourly/{1}/'.format(path, param) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(param, inmet[station][0]))
		data = arq[param]
		time = data.sel(time=slice('2018-01-01','2018-12-31'))
		var = time.groupby('time.season').mean(dim='time')
		
		if param == 'pr':
			mean_i.append(var.values*24)
		else:
			mean_i.append(var.values)

	return iy, ix, mean_i
	
	
def import_rcm_situ(param, domain, dataset):
	
	iy, ix, mean_i = [], [], []
	for station in range(1, 567):
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		iy.append(inmet[station][2])
		ix.append(inmet[station][3])

		arq = xr.open_dataset('{0}/user/mdasilva/sam_3km/postproc/'.format(path) + '{0}_{1}_{2}_mon_2018-2021_lonlat.nc'.format(param, domain, dataset))
		data = arq[param]
		time = data.sel(time=slice('2018-01-01','2018-12-31'))
		var = time.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		mean = var.groupby('time.season').mean(dim='time')
		mean_i.append(mean.values)
				
	return iy, ix, mean_i
	
	
