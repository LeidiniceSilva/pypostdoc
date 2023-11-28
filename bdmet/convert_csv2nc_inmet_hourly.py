# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 05, 2023"
__description__ = "This script convert .csv to .nc"

import os
import numpy as np
import pandas as pd

from netCDF4 import Dataset
from dict_inmet_stations import inmet

path = '/marconi/home/userexternal/mdasilva/REF/bdmet'

# choose variable: 2, 7, 8 and 21
idx=2 

if idx == 2:
	nc_var = 'pre'
	unit_var = 'mm'
	name_var = 'Hourly total of precipitation'
	std_var = 'precipitation'
elif idx == 7:
	nc_var = 'rad'
	unit_var = 'kJ.m**-2'
	name_var = 'Hourly mean of solar radiation'
	std_var = 'temperature'
elif idx == 8:
	nc_var = 'tmp'
	unit_var = 'degrees C'
	name_var = 'Hourly mean of air temperature'
	std_var = 'temperature'
else:
	nc_var = 'uv'
	unit_var = 'm.s**-1'
	name_var = 'Hourly mean of wind speed'
	std_var = 'wind'

# create date list
dt = pd.date_range('2018-01-01','2022-01-01', freq='H')
dt = dt[:-1]

for station in range(1, 567):

	if station == 15:
		continue
	if station == 23:
		continue
	if station == 47:
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
	if station == 158:
		continue
	if station == 174:
		continue
	if station == 183:
		continue
	if station == 335:
		continue
	if station == 343:
		continue
	if station == 359:
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
	if station == 529:
		continue
	if station == 566:
		continue
														
	print('Reading inmet station:', station, inmet[station][0])
	# Reading smn station
	data = pd.read_csv(os.path.join('{0}/csv/hourly/'.format(path), 'dados_{0}_H_2018-01-01_2021-12-31.csv'.format(inmet[station][0])), skiprows=9, encoding='ISO-8859-1', decimal=',', delimiter=';')
	data_i = data.iloc[:, idx]
	data_ii = data_i.replace(-999., np.nan)
	data_values = np.array(data_ii, dtype=float)
			
	data_dates  = []
	for i in range(len(dt)):
			
		data_dates.append('{0}'.format(dt[i]))
		print('Date organized:', data_dates[i], data_values[i])
		
	nc_output = '{0}/nc/hourly/{1}/{1}_{2}_H_2018-01-01_2021-12-31.nc'.format(path, nc_var, inmet[station][0])

	# create netcdf
	ds = Dataset(nc_output, mode='w', format='NETCDF4_CLASSIC')

	ds.Conventions 	= 'CF-1.6'
	ds.title 	= 'Automatic Weather stations.'
	ds.institution 	= 'Instituto Nacional de Meteorologia, INMET.'
	ds.source 	= 'INMET Meteorological Database.'
	ds.history 	= 'Rewrote via python script.'
	ds.references 	= 'https://bdmep.inmet.gov.br/.'
	ds.comment 	= 'This script convert .csv to .nc of weather station'
		
	ds.createDimension('time', None)

	time 		= ds.createVariable('time', float, ('time'))
	time.axis 	= 'L'
	time.calendar 	= 'standard'
	time.units	= 'Hours since {}'.format(data_dates[0])
	time[:]		= range(len(data_dates))

	var 		= ds.createVariable(nc_var,  float, ('time'))
	var.units 	= unit_var
	var.long_name 	= name_var
	var.standard_name = std_var
	var.missing_value = -999
	var[:] 		= data_values[:]
		
	ds.close()
	
	if os.path.exists(nc_output): 
		print('Done -->', nc_output)
exit()
