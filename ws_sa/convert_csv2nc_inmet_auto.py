# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 20, 2023"
__description__ = "This script convert .csv to .nc"

import os
import numpy as np
import pandas as pd

from netCDF4 import Dataset
from dict_inmet_auto_stations import inmet_auto

freq='daily'
path = '/marconi/home/userexternal/mdasilva/user/mdasilva/WS-SA/INMET'

skip_list = [2,15,19,23,35,47,59,64,93,96,100,105,112,117,124,137,149,152,
155,158,174,177,183,210,212,240,248,253,303,305,308,335,343,359,393,398,399,
413,417,422,426,427,444,453,457,458,479,490,495,505,514,516,529,566]

if freq == 'hourly':	
	freq_i='H'
	idx=2 # 2, 7, 8, 10 
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
else:
	freq_i='D'
	idx=10 # 1, 5, 10 
	if idx == 1:
		nc_var = 'pre'
		unit_var = 'mm'
		name_var = 'Daily total of precipitation'
		std_var = 'precipitation'
	elif idx == 5:
		nc_var = 'tmp'
		unit_var = 'degrees C'
		name_var = 'Daily mean of air temperature'
		std_var = 'temperature'
	else:
		nc_var = 'uv'
		unit_var = 'm.s**-1'
		name_var = 'Daily mean of wind speed'
		std_var = 'wind'

	
# create date list
dt = pd.date_range('2018-01-01','2022-01-01', freq='{0}'.format(freq_i))
dt = dt[:-1]

for station in range(1, 567):
	if station in skip_list:
		continue
		
	# Reading inmet station												
	print('Reading inmet station:', station, inmet_auto[station][0])
	data = pd.read_csv(os.path.join('{0}/automatic/csv/{1}/'.format(path, freq), 'dados_{0}_{1}_2018-01-01_2021-12-31.csv'.format(inmet_auto[station][0], freq_i)), skiprows=9, encoding='ISO-8859-1', decimal=',', delimiter=';')
	data_i = data.iloc[:, idx]
	data_ii = data_i.replace(-999., np.nan)
	data_values = np.array(data_ii, dtype=float)
			
	data_dates  = []
	for i in range(len(dt)):
			
		data_dates.append('{0}'.format(dt[i]))
		print('Date organized:', data_dates[i], data_values[i])
		
	nc_output = '{0}/automatic/nc/{1}/{2}/{2}_{3}_{4}_2018-01-01_2021-12-31.nc'.format(path, freq, nc_var, inmet_auto[station][0], freq_i)

	# create netcdf
	ds = Dataset(nc_output, mode='w', format='NETCDF4_CLASSIC')

	ds.Conventions 	  = 'CF-1.6'
	ds.title 	  = 'Automatic Weather stations.'
	ds.institution 	  = 'Instituto Nacional de Meteorologia, INMET.'
	ds.source 	  = 'INMET Meteorological Database.'
	ds.history 	  = 'Rewrote via python script.'
	ds.references 	  = 'https://bdmep.inmet.gov.br/.'
	ds.comment 	  = 'This script convert .csv to .nc of weather station'
		
	ds.createDimension('time', None)

	time 		  = ds.createVariable('time', float, ('time'))
	time.axis 	  = 'L'
	time.calendar 	  = 'standard'
	time.units	  = 'Hours since {}'.format(data_dates[0])
	time[:]		  = range(len(data_dates))

	var 		  = ds.createVariable(nc_var,  float, ('time'))
	var.units 	  = unit_var
	var.long_name 	  = name_var
	var.standard_name = std_var
	var.missing_value = -999
	var[:] 		  = data_values[:]
		
	ds.close()
	
	if os.path.exists(nc_output): 
		print('Done -->', nc_output)
exit()
