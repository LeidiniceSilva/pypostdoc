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

var = 'tas'
domain = 'SAM-3km'
dt = '2018-2021'  
path = '/marconi/home/userexternal/mdasilva'


def import_obs(param, dataset):

	arq   = '{0}/user/mdasilva/SAM-3km/NoTo-SAM/pdfs_v2/{1}_{2}_SAM_2018-2021_pdf.nc'.format(path, param, dataset)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:]
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = np.squeeze(var[:][0,:,0,0])

	return value


def import_cp_3km(param, dataset):

	arq   = '{0}/user/mdasilva/SAM-3km/NoTo-SAM/pdfs_v2/{1}_{2}_SAM_2018-2021_pdf.nc'.format(path, param, dataset)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:]
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = np.squeeze(var[:][0,:,0,0])
		
	return value	


# Import model and obs dataset
if var == 'pr':
	cpc = import_obs(var, 'CPC')
	regcm = import_cp_3km(var, 'RegCM5')
else:
	cpc = import_obs(var, 'CRU')
	regcm = import_cp_3km(var, 'RegCM5')
	
# Plot figure
fig = plt.figure()
font_size = 8

if var == 'pr':
	plt.plot(cpc, marker='.', markersize=4, mfc='black', mec='black', linestyle='None', label='CPC')
	plt.plot(regcm, marker='.', markersize=4, mfc='red', mec='red', linestyle='None', label='RegCM5')
	plt.xlabel('Precipitation (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
	plt.ylabel('Frequency', fontsize=font_size, fontweight='bold')
	plt.xscale('log')
	plt.yscale('log')
	plt.grid(linestyle='--')
	plt.legend(loc=1, ncol=1)
else:
	plt.plot(cpc, color='black', linestyle='--', label='CPC')
	plt.plot(regcm, color='red', linestyle='-', label='RegCM5')
	plt.xlabel('Temperature (°C)', fontsize=font_size, fontweight='bold')
	plt.ylabel('Frequency', fontsize=font_size, fontweight='bold')
	plt.grid(linestyle='--')
	plt.legend(loc=1, ncol=1)

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_pdf_{0}_{1}_RegCM5_2018-2021.png'.format(var, domain)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
