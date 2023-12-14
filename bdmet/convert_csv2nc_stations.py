# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 20, 2023"
__description__ = "This script convert .csv to .nc"

import os
import numpy as np
import pandas as pd

from netCDF4 import Dataset
from dict_inmet_stations import inmet

path = '/marconi/home/userexternal/mdasilva/OBS/BDMET/database'

skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,7793,100,105,106,107,112,117,124,135,137,139,
149,152,155,158,168,174,177,183,186,199,204,210,212,224,226,239,240,248,249,253,254,276,277,280,293,298,
303,305,306,308,319,334,335,341,343,359,362,364,384,393,396,398,399,400,402,413,416,417,422,423,426,427,
443,444,446,451,453,457,458,467,474,479,483,488,489,490,495,505,509,513,514,516,529,534,544,559,566]

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
	if station in skip_list:
		continue
	if inmet[station][2] >= -11.25235:
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
