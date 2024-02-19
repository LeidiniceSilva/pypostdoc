# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot diurnal cycle"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

var = 'pr'
domain = 'SESA-3km'
path = '/marconi/home/userexternal/mdasilva'


def import_obs(param, domain, dataset, period):

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/obs/{1}_{2}_{3}_{4}_2018-2021_lonlat.nc'.format(path, param, domain, dataset, period)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean  = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	return mean
	

def import_cp_3km(param, domain, dataset, period):

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/{1}_{2}_{3}_{4}_2018-2021_lonlat.nc'.format(path, param, domain, dataset, period)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]

	if param == 'pr':
		mean  = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	else:
		mean  = np.nanmean(np.nanmean(var[:][:,0,:,:], axis=1), axis=1)
	
	return mean
		
	
# Import model and obs dataset
if var == 'pr':
	dict_var = {'pr': ['pre', 'precip', 'sat_gauge_precip', 'tp']}
	cp_3km = import_cp_3km(var, domain, 'RegCM5', 'dc')

else:
	dict_var = {'tas': ['t2m']}
	cp_3km = import_cp_3km(var, domain, 'RegCM5', 'dc')
	
# Plot figure
fig = plt.figure()
time = np.arange(0.5, 24 + 0.5)
font_size = 8

if var == 'pr':
	plt1 = plt.plot(time, cp_3km)
	plt.title(u'a) SESA', loc='left', fontweight='bold', fontsize=8)
	l1 = plt1
	plt.setp(l1, linewidth=1., linestyle='-', markersize=3, marker='o', markerfacecolor='white', color='red')
	plt.ylim(0, 8)
	plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
	plt.yticks(np.arange(0, 9, 1), fontsize=font_size)
	plt.xlabel('Hours', fontsize=font_size, fontweight='bold')
	plt.ylabel('Precipitation (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
	plt.grid(linestyle='--')
	plt.legend(plt1, ['RegCM5'], fontsize=font_size, ncol=1, loc=1, shadow=True)
else:
	plt1 = plt.plot(time, cp_3km)
	plt.title(u'a) SESA', loc='left', fontweight='bold', fontsize=8)
	l1 = plt1
	plt.setp(l1, linewidth=1., linestyle='-', markersize=3, marker='o', markerfacecolor='white', color='red')
	plt.ylim(8, 30)
	plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
	plt.yticks(np.arange(8, 32, 2), fontsize=font_size)
	plt.xlabel('Hours', fontsize=font_size, fontweight='bold')
	plt.ylabel('Temperature (°C)', fontsize=font_size, fontweight='bold')
	plt.grid(linestyle='--')
	plt.legend(plt1, ['RegCM5'], fontsize=font_size, ncol=1, loc=1, shadow=True)

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km_v1/figs/evaluate'.format(path)
name_out = 'pyplt_diurnal_cycle_{0}_{1}_RegCM5_2018-2021.png'.format(var, domain)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
