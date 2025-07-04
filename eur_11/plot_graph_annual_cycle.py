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
dt = '2000-2009'
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11'


def import_obs(param, dataset):

	arq   = '{0}/postproc/obs/{1}_{2}_{3}-FPS_mon_{4}_lonlat.nc'.format(path, param, dataset, domain, dt) 
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

	arq   = '{0}/postproc/rcm/{1}_{2}-FPS_mon_{3}_lonlat.nc'.format(path, param, dataset, dt)	
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
noto = import_rcm(var, 'RegCM5_NoTo-EUR')
wdm7 = import_rcm(var, 'RegCM5_WDM7-EUR')
wsm7 = import_rcm(var, 'RegCM5_WSM7-EUR')
wsm5 = import_rcm(var, 'RegCM5_WSM5-EUR')

# Plot figure
fig = plt.figure(figsize=(10, 5))
time = np.arange(1, 12 + 1)
font_size = 10

ax = fig.add_subplot(1, 1, 1)  
plt.plot(time, eobs, label='EOBS', color='black', linewidth=1., marker='.', linestyle='-')
plt.plot(time, noto, label='NoTo', color='blue', linewidth=1., marker='.', linestyle='-')
plt.plot(time, wdm7, label='WDM7', color='green', linewidth=1., marker='.', linestyle='-')
plt.plot(time, wsm7, label='WSM7', color='magenta', linewidth=1., marker='.', linestyle='-')
plt.plot(time, wsm5, label='WSM5', color='red', linewidth=1., marker='.', linestyle='-')

plt.title('(a)', loc='left', fontsize=font_size, fontweight='bold') 
plt.xlabel('Annual cycle ({0})'.format(dt), fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.xlim(0, 13)
plt.ylim(0, 7)
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'), fontsize=font_size)
plt.yticks(np.arange(0, 7.5, 0.5), fontsize=font_size)
plt.grid(linestyle='--')
plt.legend(loc=2, ncol=5, fontsize=font_size, shadow=True)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_graph_annual_cycle_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

