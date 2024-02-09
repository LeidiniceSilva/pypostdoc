# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot bias maps"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from import_climate_tools import compute_pdf

var = 'pr'
domain = 'SESA-3km'
dt = '2018-2021'  
path = '/marconi/home/userexternal/mdasilva'


def import_obs(param, domain, dataset, season):

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/obs/{1}_{2}_{3}_{4}_2018-2021_lonlat.nc'.format(path, param, domain, dataset, season)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:]
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	one_d_list = []
	for i in range(len(value)):
		for j in range(len(value[i])):
			for k in range(len(value[i][j])):
				one_d_list.append(value[i][j][k])
	return value


def import_cp_3km(param, domain, dataset, season):

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/{1}_{2}_{3}_{4}_2018-2021_lonlat.nc'.format(path, param, domain, dataset, season)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:]
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][0:30,:,:]

	one_d_list = []
	for i in range(len(value)):
		for j in range(len(value[i])):
			for k in range(len(value[i][j])):
				one_d_list.append(value[i][j][k])
	
	return value	

# Import model and obs dataset
if var == 'pr':
	dict_var = {'pr': ['pre', 'pre', 'precip', 'sat_gauge_precip', 'tp']}

	cpc = import_obs(dict_var[var][2], domain, 'CPC', 'day')
	regcm = import_cp_3km(var, domain, 'RegCM5', 'day')
	
	# Round values
	round_cpc = np.round(cpc, 0)
	round_regcm = np.round(regcm, 0)
	
	# Filter 0 mm/day
	filter_cpc = round_cpc[cpc > 0.]
	filter_regcm = round_regcm[regcm > 0.]

	# Compute pdf function
	x_pdf_cpc, pdf_cpc = np.unique(filter_cpc, return_counts=True)
	x_pdf_regcm, pdf_regcm = np.unique(filter_regcm, return_counts=True)
else:
	dict_var = {'tas': ['tmax', 'tmin']}

	cpc_tmax = import_obs(dict_var[var][0], domain, 'CPC', 'mon')
	cpc_tmin = import_obs(dict_var[var][1], domain, 'CPC', 'mon')
	cpc = (cpc_tmax+cpc_tmin)/2

	regcm = import_cp_3km(var, domain, 'RegCM5', 'mon')
	regcm = np.nanmean(regcm, axis=1)
	
	# import pdf function
	x_pdf_cpc, pdf_cpc = compute_pdf(cpc)
	x_pdf_regcm, pdf_regcm = compute_pdf(regcm)

# Plot figure
fig = plt.figure()
font_size = 8

if var == 'pr':
	plt.plot(x_pdf_cpc, pdf_cpc, marker='.', markersize=4, mfc='black', mec='black', linestyle='None', label='CPC')
	plt.plot(x_pdf_regcm, pdf_regcm, marker='.', markersize=4, mfc='red', mec='red', linestyle='None', label='RegCM5')
	plt.xlabel('Intensity of  precipitation (> 0 mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
	plt.ylabel('Frequency', fontsize=font_size, fontweight='bold')
	plt.yscale('log')
	plt.xlim(0,200)
	plt.grid(linestyle='--')
	plt.legend(loc=1, ncol=1)
else:
	plt.plot(x_pdf_cpc, pdf_cpc, color='black', linestyle='--', label='CPC')
	plt.plot(x_pdf_regcm, pdf_regcm, color='red', linestyle='-', label='RegCM5')
	plt.xlabel('Temperature (°C)', fontsize=font_size, fontweight='bold')
	plt.ylabel('Frequency', fontsize=font_size, fontweight='bold')
	plt.xlim(10,30)
	plt.ylim(0,0.2)
	plt.grid(linestyle='--')
	plt.legend(loc=1, ncol=1)

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_frequency_{0}_CP-RegCM5_{1}_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
