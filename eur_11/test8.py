# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot clim maps"

import os
import netCDF4
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap

var = 'pr'
dt = '20000101'
domain = 'EUR-11'
path = '/home/mda_silv/scratch/EUR-11/postproc'


def import_obs(param, dataset):

	arq   = '{0}/obs/{1}_{2}_FPS_{3}_{4}_lonlat.nc'.format(path, param, domain, dataset, dt) 
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return mean


def import_rcm(param, dataset):

	arq   = '{0}/rcm/{1}_{2}_FPS_{3}_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:-1,:,:]
	mean = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return mean
	
	
# Import model and obs dataset
dict_var = {'pr': ['precip', 'rr', 'tp']}

cpc_jan = import_obs(dict_var[var][0], 'CPC')
noto_jan = import_rcm(var, 'NoTo-Europe')
wsm5_jan = import_rcm(var, 'WSM5-Europe')
wsm7_jan = import_rcm(var, 'WSM7-Europe')
wdm7_jan = import_rcm(var, 'WDM7-Europe')

# Plot figure
fig = plt.figure(figsize=(10, 5))
time = np.arange(1, 31 + 1)
font_size = 10

ax = fig.add_subplot(1, 1, 1)  
plt.plot(time, cpc_jan, label='CPC', color='black', linewidth=1., marker='.', linestyle='-')
plt.plot(time, noto_jan, label='NoTo', color='blue', linewidth=1., marker='.', linestyle='-')
plt.plot(time, wsm5_jan, label='WSM5', color='red', linewidth=1., marker='.', linestyle='-')
plt.plot(time, wsm7_jan, label='WSM7', color='magenta', linewidth=1., marker='.', linestyle='-')
plt.plot(time, wdm7_jan, label='WDM7', color='green', linewidth=1., marker='.', linestyle='-')

plt.title('(a)', loc='left', fontsize=font_size, fontweight='bold') 
plt.xlabel('Jan 2000', fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.xlim(0, 32)
plt.ylim(0, 11)
plt.xticks(np.arange(0, 32, 1), fontsize=font_size)
plt.yticks(np.arange(0, 11.5, 0.5), fontsize=font_size)
plt.grid(linestyle='--')
plt.legend(loc=2, ncol=1, fontsize=font_size, shadow=True)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_timesieries_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

