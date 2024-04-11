# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot pdf function"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

var = 'pr'
domain = 'SESA-3km'
idt, fdt = '2018', '2021'
dt = '{0}-{1}'.format(idt, fdt)
 
path = '/marconi/home/userexternal/mdasilva'


def import_obs(param, dataset, season):

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/obs/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:]
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	value_ = value.flatten()

	return value_


def import_rcm(param, dataset, season):

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:]
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	value_ = value.flatten()

	return value_	


# Import model and obs dataset
if var == 'pr':
	obs = import_obs('precip', 'CPC', 'day')
	regcm = import_rcm(var, 'RegCM5', 'day')
else:
	obs = import_obs('tmp', 'CRU', 'mon')
	regcm = import_rcm(var, 'RegCM5', 'mon')

# Plot figure
fig = plt.figure()
font_size = 8

if var == 'pr':
	plt.plot(x_pdf_obs, pdf_obs, marker='.', markersize=4, mfc='black', mec='black', linestyle='None', label='CPC')
	plt.plot(x_pdf_regcm, pdf_regcm, marker='.', markersize=4, mfc='red', mec='red', linestyle='None', label='RegCM5')
	plt.xlabel('Precipitation (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
	plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
	plt.xscale('log')
	plt.yscale('log')
	plt.grid(linestyle='--')
	plt.legend(loc=1, ncol=1)
else:
	plt.plot(x_pdf_obs, pdf_obs, color='black', linestyle='--', label='CPC')
	plt.plot(x_pdf_regcm, pdf_regcm, color='red', linestyle='-', label='RegCM5')
	plt.xlabel('Temperature (°C)', fontsize=font_size, fontweight='bold')
	plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
	plt.yscale('log')
	plt.grid(linestyle='--')
	plt.legend(loc=1, ncol=1)

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_pdf_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
