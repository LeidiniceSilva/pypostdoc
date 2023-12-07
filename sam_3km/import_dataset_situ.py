# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script import in-situ datasets"

import xarray as xr

from dict_inmet_stations import inmet

path = '/marconi/home/userexternal/mdasilva'

	
def import_obs_situ(param):
	
	iy, ix, mean_i = [], [], []
	for station in range(1, 567):

		if station == 2:
			continue
		if station == 15:
			continue
		if station == 19:
			continue
		if station == 23:
			continue
		if station == 35:
			continue
		if station == 47:
			continue
		if station == 59:
			continue
		if station == 64:
			continue
		if station == 93:
			continue
		if station == 96:
			continue
		if station == 100:
			continue
		if station == 105:
			continue
		if station == 112:
			continue
		if station == 117:
			continue
		if station == 124:
			continue
		if station == 137:
			continue
		if station == 149:
			continue
		if station == 152:
			continue
		if station == 155:
			continue
		if station == 158:
			continue
		if station == 174:
			continue
		if station == 183:
			continue
		if station == 210:
			continue
		if station == 212:
			continue
		if station == 240:
			continue
		if station == 248:
			continue
		if station == 253:
			continue
		if station == 303:
			continue
		if station == 305:
			continue
		if station == 308:
			continue
		if station == 335:
			continue
		if station == 343:
			continue
		if station == 359:
			continue
		if station == 393:
			continue
		if station == 398:
			continue
		if station == 399:
			continue
		if station == 413:
			continue
		if station == 417:
			continue
		if station == 422:
			continue
		if station == 426:
			continue
		if station == 427:
			continue
		if station == 444:
			continue
		if station == 453:
			continue
		if station == 457:
			continue
		if station == 458:
			continue
		if station == 479:
			continue
		if station == 490:
			continue
		if station == 495:
			continue
		if station == 505:
			continue
		if station == 514:
			continue
		if station == 529:
			continue
		if station == 566:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		iy.append(inmet[station][2])
		ix.append(inmet[station][3])

		arq = xr.open_dataset('{0}/OBS/BDMET/nc/hourly/pre/'.format(path) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(param, inmet[station][0]))
		data = arq[param]
		time = data.sel(time=slice('2018-01-01','2021-12-31'))
		var = time.groupby('time.season').mean(dim='time')
		mean_i.append(var.values*24)
				
	return iy, ix, mean_i
	
	
def import_rcm_situ(param):
	
	iy, ix, mean_i = [], [], []
	for station in range(1, 567):

		if station == 2:
			continue
		if station == 15:
			continue
		if station == 19:
			continue
		if station == 23:
			continue
		if station == 35:
			continue
		if station == 47:
			continue
		if station == 59:
			continue
		if station == 64:
			continue
		if station == 93:
			continue
		if station == 96:
			continue
		if station == 100:
			continue
		if station == 105:
			continue
		if station == 112:
			continue
		if station == 117:
			continue
		if station == 124:
			continue
		if station == 137:
			continue
		if station == 149:
			continue
		if station == 152:
			continue
		if station == 155:
			continue
		if station == 158:
			continue
		if station == 174:
			continue
		if station == 183:
			continue
		if station == 210:
			continue
		if station == 212:
			continue
		if station == 240:
			continue
		if station == 248:
			continue
		if station == 253:
			continue
		if station == 303:
			continue
		if station == 305:
			continue
		if station == 308:
			continue
		if station == 335:
			continue
		if station == 343:
			continue
		if station == 359:
			continue
		if station == 393:
			continue
		if station == 398:
			continue
		if station == 399:
			continue
		if station == 413:
			continue
		if station == 417:
			continue
		if station == 422:
			continue
		if station == 426:
			continue
		if station == 427:
			continue
		if station == 444:
			continue
		if station == 453:
			continue
		if station == 457:
			continue
		if station == 458:
			continue
		if station == 479:
			continue
		if station == 490:
			continue
		if station == 495:
			continue
		if station == 505:
			continue
		if station == 514:
			continue
		if station == 529:
			continue
		if station == 566:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		iy.append(inmet[station][2])
		ix.append(inmet[station][3])

		arq = xr.open_dataset('{0}/user/mdasilva/sam_3km/postproc/'.format(path) + '{0}_SAM-3km_RegCM5_mon_2018-2021_lonlat.nc'.format(param))
		data = arq[param]
		time = data.sel(time=slice('2018-01-01','2021-12-31'))
		var = time.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		mean = var.groupby('time.season').mean(dim='time')
		mean_i.append(mean.values)
				
	return iy, ix, mean_i
	
	
