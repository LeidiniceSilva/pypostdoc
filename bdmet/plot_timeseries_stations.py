# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 05, 2023"
__description__ = "This script plot weather stations timeseries"

import os
import math
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet

path = '/marconi/home/userexternal/mdasilva'

# Import latitude, longitude and database
for i in range(1, 567):

	if i == 15:
		continue
	if i == 23:
		continue
	if i == 47:
		continue
	if i == 105:
		continue
	if i == 112:
		continue
	if i == 117:
		continue
	if i == 124:
		continue
	if i == 137:
		continue
	if i == 149:
		continue
	if i == 158:
		continue
	if i == 174:
		continue
	if i == 183:
		continue
	if i == 335:
		continue
	if i == 343:
		continue
	if i == 359:
		continue
	if i == 398:
		continue
	if i == 399:
		continue
	if i == 413:
		continue
	if i == 417:
		continue
	if i == 422:
		continue
	if i == 426:
		continue
	if i == 444:
		continue
	if i == 453:
		continue
	if i == 457:
		continue
	if i == 458:
		continue
	if i == 479:
		continue
	if i == 490:
		continue
	if i == 495:
		continue
	if i == 505:
		continue
	if i == 529:
		continue
	if i == 566:
		continue
		
	yy = inmet[i][2]
	xx = inmet[i][3]
	name = inmet[i][0]

	print('Reading weather station:', i, name)
	d_i = xr.open_dataset('{0}/REF/bdmet/nc/hourly/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
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
	path_out = '{0}/figs/bdmet'.format(path)
	name_out = 'pyplt_ts_pr_inmet_{0}_2018-2021.png'.format(name)
	plt.savefig(os.path.join(path_out, name_out), dpi=100, bbox_inches='tight')
	plt.close('all')
	plt.cla()
exit()

