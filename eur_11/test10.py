# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot annual cycle"

import os
import netCDF4
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt

var = 'pr'
domain = 'EUR-11'
dt = '1970-1973'
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11'


def import_obs(param, dataset):

	arq   = '{0}/postproc/obs/{1}_{2}_FPS_{3}_mon_{4}_lonlat.nc'.format(path, param, domain, dataset, dt) 
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean = np.nanmean(np.nanmean(value, axis=1), axis=1)

	clim = []
	for mon in range(0, 11 + 1):
		mean_ = np.nanmean(mean[mon::12], axis=0)
		clim.append(mean_)

	return clim


def import_rcm(param, dataset):

	arq   = '{0}/postproc/rcm/{1}_{2}_FPS_{3}_mon_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean = np.nanmean(np.nanmean(value, axis=1), axis=1)

	clim = []
	for mon in range(0, 11 + 1):
		mean_ = np.nanmean(mean[mon::12], axis=0)
		clim.append(mean_)

	return clim
	
	
# Import model and obs dataset
dict_var = {'pr': ['precip', 'rr', 'tp']}

eobs = import_obs(dict_var[var][1], 'EOBS')
noto_v1 = import_rcm(var, 'NoTo-Europe_RegCM5')
noto_v2 = import_rcm(var, 'NoTo-Europe_cordex5_RegCM5')
wdm7 = import_rcm(var, 'WDM7-Europe_RegCM5')
wsm7 = import_rcm(var, 'WSM7-Europe_RegCM5')
wsm5 = import_rcm(var, 'WSM5-Europe_RegCM5')

# Plot figure
fig = plt.figure(figsize=(10, 5))
time = np.arange(1, 12 + 1)
font_size = 10

ax = fig.add_subplot(1, 1, 1)  
plt.plot(time, eobs, label='EOBS', color='black', linewidth=1., marker='.', linestyle='-')
plt.plot(time, noto_v1, label='NoTo (CORE)', color='blue', linewidth=1., marker='.', linestyle='-')
plt.plot(time, noto_v2, label='NoTo (CORDEX5)', color='gray', linewidth=1., marker='.', linestyle='-')
plt.plot(time, wdm7, label='WDM7', color='green', linewidth=1., marker='.', linestyle='-')
plt.plot(time, wsm7, label='WSM7', color='magenta', linewidth=1., marker='.', linestyle='-')
plt.plot(time, wsm5, label='WSM5', color='red', linewidth=1., marker='.', linestyle='-')

plt.title('(a)', loc='left', fontsize=font_size, fontweight='bold') 
plt.xlabel('Months (1970-1973)', fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.xlim(0, 13)
plt.ylim(0, 8)
plt.xticks(np.arange(0, 13, 1), fontsize=font_size)
plt.yticks(np.arange(0, 8.5, 0.5), fontsize=font_size)
plt.grid(linestyle='--')
plt.legend(loc=2, ncol=6, fontsize=font_size, shadow=True)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_graph_annual_cycle_{0}_{1}_RegCM5_{2}_v2.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

