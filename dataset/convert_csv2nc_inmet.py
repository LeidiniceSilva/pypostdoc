# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 05, 2023"
__description__ = "This script convert .csv to .nc"

import os
import numpy as np
import pandas as pd

from netCDF4 import Dataset

idx=1
if idx == 1:
	nc_var = 'pre'
	unit_var = 'mm'
	name_var = 'Daily total of precipitation'
	std_var = 'precipitation'

elif idx == 6:
	nc_var = 'tmp'
	unit_var = 'C'
	name_var = 'Daily mean of air temperature'
	std_var = 'temperature'
else:
	nc_var = 'uv'
	unit_var = 'm s**-1'
	name_var = 'Daily mean of wind speed'
	std_var = 'wind'

path = '/home/nice/Documentos/FPS_SESA/database/obs/inmet/BDMET/daily'

# create date list
dt = pd.date_range('2018-01-01','2021-12-31', freq='D').strftime('%Y-%m-%d').tolist()
	
for station in range(1, 2):
																																																						
	print('Reading inmet station:', station)
	# Reading smn station
	data = pd.read_csv(os.path.join('{0}'.format(path), 'dados_A001_D_2018-01-01_2021-12-31.csv'), skiprows=9, encoding='ISO-8859-1', decimal=',', delimiter=';')
	data = data.iloc[:, idx]
	data_values = data.replace(-999., np.nan)
	data_values = np.array(data_values, dtype=float)

	data_dates  = []
	for i in range(len(dt)):
			
		data_dates.append('{0}'.format(dt[i]))
		print('Date organized:', data_dates[i], data_values[i])
		
	nc_output = '{0}/{1}_A001_D_2018-01-01_2021-12-31.nc'.format(path, nc_var)

	# create netcdf
	ds = Dataset(nc_output, mode='w', format='NETCDF4_CLASSIC')

	ds.Conventions 	= 'CF-1.6'
	ds.title 		= 'Weather stations.'
	ds.institution 	= 'INMET.'
	ds.source 		= 'Automatic weather station.'
	ds.history 		= 'Rewrote via python script.'
	ds.references 	= 'https://bdmep.inmet.gov.br/.'
	ds.comment 		= 'This script convert .csv to .nc of weather station'
		
	ds.createDimension('time', None)

	time 				= ds.createVariable('time', float, ('time'))
	time.axis 			= 'L'
	time.calendar 		= 'standard'
	time.units			= 'days since {}'.format(data_dates[0])
	time[:]				= range(len(data_dates))

	var 				= ds.createVariable(nc_var,  float, ('time'))
	var.units 			= unit_var
	var.long_name 		= name_var
	var.standard_name 	= std_var
	var.missing_value 	= -999
	var[:] 				= data_values[:]
		
	ds.close()
	
	if os.path.exists(nc_output): 
		print('Done -->', nc_output)

exit()
