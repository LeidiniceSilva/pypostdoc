# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/09/2023"
__description__ = "This script check nan in automatic weather station"

import os
import math
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from dict_sesa_inmet_stations import inmet
from dict_urug_smn_stations import urug_smn
from dict_arg_emas_stations import arg_emas

dt = 'H_2018-01-01_2021-12-31'

for i in range(1, 101):
	
	print('Reading weather station:', i, inmet[i][0], inmet[i][1])

	# Reading inmet 
	d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'pre_{0}_{1}.nc'.format(inmet[i][0], dt))
	d_i = d_i.pre.sel(time=slice('2019-01-01','2021-12-31'))
	d_i = d_i.values
	x = sum(math.isnan(x) for x in d_i)
	percent_missing = x * 100 / len(d_i)
		
for j in range(1, 72):

	print('Reading Uruguai weather station:', j, urug_smn[j][0])	
	
	# Reading Uruguai weather stations
	d_j = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/urug_smn/urug_smn_nc/' + 'pre_{0}_{1}.nc'.format(urug_smn[j][0], dt))
	d_j = d_j.pre.sel(time=slice('2019-01-01','2021-12-31'))
	d_j = d_j.values	
	x = sum(math.isnan(x) for x in d_j)
	percent_missing = x * 100 / len(d_j)

for k in range(1, 88):

	if k == 2:
		continue
	if k == 4:
		continue
	if k == 9:
		continue
	if k == 17:
		continue
	if k == 24:
		continue
	if k == 26:
		continue
	if k == 28:
		continue
	if k == 31:
		continue
	if k == 32:
		continue
	if k == 40:
		continue
	if k == 43:
		continue
	if k == 57:
		continue
	if k == 67:
		continue
	if k == 72:
		continue
	if k == 76:
		continue
	if k == 81:
		continue

	print('Reading Argentina weather station:', k, arg_emas[k][0])	
	
	# Reading Argentina weather stations
	d_k = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/arg_emas/arg_emas_nc/' + 'precip_{0}_{1}.nc'.format(arg_emas[k][0], dt))
	d_k = d_k.precip.sel(time=slice('2019-01-01','2021-12-31'))
	d_k = d_k.values
	print(len(d_k))
	x = sum(math.isnan(x) for x in d_k)
	print(x)
	percent_missing = x * 100 / len(d_k)
	print(np.around(percent_missing),'%')


