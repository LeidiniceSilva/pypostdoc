# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot diurnal cycle"

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

	arq   = '{0}/postproc/obs/{1}_{2}_{3}-FPS_diurnal_cycle_{4}_lonlat.nc'.format(path, param, dataset, domain, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return mean


def import_rcm(param, dataset):

	arq   = '{0}/postproc/rcm/{1}_{2}-FPS_diurnal_cycle_{3}_lonlat.nc'.format(path, param, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return mean
	
	
# Import model and obs dataset
dict_var = {'pr': ['precip', 'rr', 'tp']}

era5 = import_obs(dict_var[var][2], 'ERA5')
noto = import_rcm(var, 'RegCM5_NoTo-EUR')
wdm7 = import_rcm(var, 'RegCM5_WDM7-EUR')
wsm7 = import_rcm(var, 'RegCM5_WSM7-EUR')
wsm5 = import_rcm(var, 'RegCM5_WSM5-EUR')

# Plot figure
fig = plt.figure(figsize=(10, 5))
font_size = 10
time = np.arange(1, 24 + 1)

ax = fig.add_subplot(1, 1, 1)  
plt.plot(time, era5*24, label='ERA5', color='black', linewidth=1., marker='.', linestyle='-')
plt.plot(time, noto*24, label='NoTo', color='blue', linewidth=1., marker='.', linestyle='-')
plt.plot(time, wsm5*24, label='WSM5', color='red', linewidth=1., marker='.', linestyle='-')
plt.plot(time, wsm7*24, label='WSM7', color='magenta', linewidth=1., marker='.', linestyle='-')
plt.plot(time, wdm7*24, label='WDM7', color='green', linewidth=1., marker='.', linestyle='-')

plt.title('(b)', loc='left', fontsize=font_size, fontweight='bold') 
plt.xlabel('Diurnal cycle ({0})'.format(dt), fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.ylim(1, 3.5)
plt.xlim(0, 25)
plt.yticks(np.arange(1, 3.75, 0.25), fontsize=font_size)
plt.xticks(time, ('00hr', '', '02hr', '', '04hr', '', '06hr', '', '08hr', '', '10hr', '', '12hr', '', '14hr', '', '16hr', '', '18hr', '', '20hr', '', '22hr', ''), fontsize=font_size)
plt.grid(linestyle='--')
plt.legend(loc=2, ncol=3, fontsize=font_size, shadow=True)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_graph_diurnal_cycle_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()




