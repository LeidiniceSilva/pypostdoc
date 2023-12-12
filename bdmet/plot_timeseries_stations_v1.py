# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 20, 2023"
__description__ = "This script plot weather station timeseries"

import os
import math
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet

path='/afs/ictp.it/home/m/mda_silv/Documents/BDMET'

# Import latitude, longitude and database
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
		
	yy = inmet[station][2]
	xx = inmet[station][3]
	name = inmet[station][0]

	print('Reading weather station:', station, name)
	d_i = xr.open_dataset('{0}/database/nc/hourly/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[station][0]))
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
	path_out = '{0}/figs/ts_v1'.format(path)
	name_out = 'pyplt_ts_pr_inmet_{0}_2018-2021.png'.format(name)
	plt.savefig(os.path.join(path_out, name_out), dpi=100, bbox_inches='tight')
	plt.close('all')
	plt.cla()
exit()
