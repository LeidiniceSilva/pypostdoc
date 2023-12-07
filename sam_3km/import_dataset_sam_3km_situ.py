# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script import in-situ datasets"

import os
import netCDF4
import numpy as np

path = '/marconi/home/userexternal/mdasilva'

	
def import_stations(param, dataset):
	
	iy, ix = [], []
	mean_i, mean_ii = [], []

	# Select lat and lon 
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
		if inmet[station][2] <= 
			
		yy = inmet[station][2]
		xx = inmet[station][3]
		name = inmet[station][0]
		
		print('Reading weather station:', station, inmet[station][0])		
		# reading regcm usp 
		d_i = xr.open_dataset('{0}/user/mdasilva/sam_3km/postproc/'.format(path) + '{0}_SAM-3km_RegCM5_mon_2018-2021_lonlat.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		d_i = d_i.groupby('time.year').mean('time')
		d_i = np.nanmean(d_i.values)
		mean_i.append(d_i)

		# Reading inmet 
		d_vi = xr.open_dataset('{0}/OBS/BDMET/nc/hourly/pre/'.format(path) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.groupby('time.year').mean('time')
		d_vi = np.nanmean(d_vi.values)
		mean_vi.append(d_vi*24)
		
				
	return iy, ix, mean_i, mean_ii
