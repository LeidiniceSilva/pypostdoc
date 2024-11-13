# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot clim maps"

import os
import netCDF4
import numpy as np
import xarray as xr
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap

var = 'pr'
dt = '2000010100'
domain = 'EUR-11'
path = '/home/mda_silv/scratch/EUR-11/postproc'


def import_obs(param, dataset):
	
	arq  = xr.open_dataset('{0}/obs/{1}_{2}_FPS_{3}_{4}_lonlat.nc'.format(path, param, domain, dataset, dt))
	data = arq[param]
	time = data.sel(time=slice('2000-01-01','2000-01-31'))
	var  = time.groupby('time.hour').mean()
	mean = np.nanmean(np.nanmean(var.values, axis=1), axis=1)

	return mean


def import_rcm(param, dataset):

	arq  = xr.open_dataset('{0}/rcm/{1}_{2}_FPS_{3}_{4}_lonlat.nc'.format(path, param, domain, dataset, dt))
	data = arq[param]
	time = data.sel(time=slice('2000-01-01','2000-01-31'))
	var  = time.groupby('time.hour').mean()
	mean = np.nanmean(np.nanmean(var.values, axis=1), axis=1)

	return mean
	
	
# Import model and obs dataset
dict_var = {'pr': ['precip', 'rr', 'tp']}

era5_jan = import_obs(dict_var[var][2], 'ERA5')
noto_jan = import_rcm(var, 'NoTo-Europe')
wsm5_jan = import_rcm(var, 'WSM5-Europe')
wsm7_jan = import_rcm(var, 'WSM7-Europe')
wdm7_jan = import_rcm(var, 'WDM7-Europe')

# Plot figure
fig = plt.figure(figsize=(10, 5))
font_size = 10
time = np.arange(1, 24 + 1)

ax = fig.add_subplot(1, 1, 1)  
plt.plot(time, era5_jan*24, label='ERA5', color='black', linewidth=1., marker='.', linestyle='-')
plt.plot(time, noto_jan*24, label='NoTo', color='blue', linewidth=1., marker='.', linestyle='-')
plt.plot(time, wsm5_jan*24, label='WSM5', color='red', linewidth=1., marker='.', linestyle='-')
plt.plot(time, wsm7_jan*24, label='WSM7', color='magenta', linewidth=1., marker='.', linestyle='-')
plt.plot(time, wdm7_jan*24, label='WDM7', color='green', linewidth=1., marker='.', linestyle='-')

plt.title('(b)', loc='left', fontsize=font_size, fontweight='bold') 
plt.xlabel('Hours', fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.ylim(1, 3.5)
plt.xlim(0, 25)
plt.yticks(np.arange(1, 3.75, 0.25), fontsize=font_size)
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.grid(linestyle='--')
plt.legend(loc=10, ncol=5, fontsize=font_size, shadow=True)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_diurnal_cycle_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()




