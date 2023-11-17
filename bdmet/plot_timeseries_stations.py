# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 16, 2023"
__description__ = "This script plot weather stations timeseries"

import os
import math
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i

stations = 'smn_i' # 100 or 73
path = '/afs/ictp.it/home/m/mda_silv/Documents/FPS_SESA'

# Import latitude, longitude and database
for i in range(1, 73):
	
	if stations == 'inmet':
		yy = inmet[i][2]
		xx = inmet[i][3]
		name = inmet[i][0]
	else:
		yy = smn_i[i][1]
		xx = smn_i[i][2]
		name = smn_i[i][0]

	print('Reading weather station:', stations, name)
	if stations == 'inmet':
		d_i = xr.open_dataset('{0}/database/obs/inmet/inmet_nc_sesa/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))	
		d_i = d_i.values	
	else:
		d_i = xr.open_dataset('{0}/database/obs/smn_i/smn_nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[i][0]))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))	
		d_i = d_i.values
		
	# Plot figure
	fig = plt.figure()
	
	ts = pd.date_range(start="20180101", end="20220101", freq="H")
	ts_h = pd.Series(data=d_i, index=ts[:-1])
	
	x = sum(math.isnan(x) for x in d_i)
	percent_missing = x * 100 / len(d_i)

	ax = fig.add_subplot(1, 1, 1)
	plt.plot(ts_h, linewidth=1, color='blue')
	plt.title(u'{0} NaN={1}% \n lat: {2} lon: {3}'.format(name, round(percent_missing), yy, xx), fontweight='bold')
	plt.xlabel(u'2018-01-01 00h - 2021-12-31 23h', fontweight='bold')
	plt.ylabel(u'Precipitation (mm h⁻¹)', fontweight='bold')
	plt.ylim(0, 200)
	plt.grid()
					
	# Path out to save figure
	path_out = '{0}/figs/sesa_v2'.format(path)
	name_out = 'pyplt_ts_station_{0}_{1}_pr.png'.format(stations, name)
	plt.savefig(os.path.join(path_out, name_out), dpi=100, bbox_inches='tight')
	plt.close('all')
	plt.cla()
exit()

