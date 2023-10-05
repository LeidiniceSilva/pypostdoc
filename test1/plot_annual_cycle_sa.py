# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "September 16, 2023"
__description__ = "This script plot annual cycle"

import os
import netCDF4
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from mpl_toolkits.basemap import Basemap


def import_ref(param, ref, min_lon, min_lat, max_lon, max_lat):

	arq = xr.open_dataset('/home/nice/Downloads/' + '{0}_{1}_mon_2000_lonlat.nc'.format(param, ref))
	data = arq[param]
	data = data.sel(lat=slice(min_lat,max_lat), lon=slice(min_lon,max_lon))
	time = data.sel(time=slice('2000-01-01','2000-12-31'))	
	var = time.values
	mean = np.nanmean(np.nanmean(var, axis=1), axis=1)

	return mean
	
	
def import_rcm(param, min_lon, min_lat, max_lon, max_lat):

	arq = xr.open_dataset('/home/nice/Downloads/' + '{0}_SAM-22_Reg5_mon_2000_lonlat.nc'.format(param))
	data = arq[param]
	data = data.sel(lat=slice(min_lat,max_lat), lon=slice(min_lon,max_lon))
	time = data.sel(time=slice('2000-01-01','2000-12-31'))
	var = time.values
	mean = np.nanmean(np.nanmean(var, axis=1), axis=1)
	# mean = np.nanmean(mean, axis=1)
	
	return mean
	

# Import cmip models and obs database 
var_obs = 'pre'
var_rea = 'tp'
var_rcm = 'pr'

samz_min_lon, samz_min_lat, samz_max_lon, samz_max_lat = -70,-12.5,-50,-5
lpb_min_lon, lpb_min_lat, lpb_max_lon, lpb_max_lat = -63,-32.5,-49,-20
neb_min_lon, neb_min_lat, neb_max_lon, neb_max_lat = -46,-15,-35,-3

clim_samz_obs = import_ref(var_obs, 'cru_ts4.07', samz_min_lon, samz_min_lat, samz_max_lon, samz_max_lat)
clim_lpb_obs = import_ref(var_obs, 'cru_ts4.07', lpb_min_lon, lpb_min_lat, lpb_max_lon, lpb_max_lat)
clim_neb_obs = import_ref(var_obs, 'cru_ts4.07', neb_min_lon, neb_min_lat, neb_max_lon, neb_max_lat)

clim_samz_rea = import_ref(var_rea, 'era5', samz_min_lon, samz_min_lat, samz_max_lon, samz_max_lat)
clim_lpb_rea = import_ref(var_rea, 'era5', lpb_min_lon, lpb_min_lat, lpb_max_lon, lpb_max_lat)
clim_neb_rea = import_ref(var_rea, 'era5', neb_min_lon, neb_min_lat, neb_max_lon, neb_max_lat)

clim_samz_rcm = import_rcm(var_rcm, samz_min_lon, samz_min_lat, samz_max_lon, samz_max_lat)
clim_lpb_rcm = import_rcm(var_rcm, lpb_min_lon, lpb_min_lat, lpb_max_lon, lpb_max_lat)
clim_neb_rcm = import_rcm(var_rcm, neb_min_lon, neb_min_lat, neb_max_lon, neb_max_lat)

# Plot cmip models and obs database
fig = plt.figure()
time = np.arange(0.5, 12+0.5)

ax = fig.add_subplot(3, 1, 1)
annual_cycle = ax.plot(time, clim_samz_rcm, time, clim_samz_rea, time, clim_samz_obs)
plt.title(u'(a) AMZ', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
if var_rcm == 'pr':
	plt.yticks(np.arange(0, 18, 2), fontsize=8)
	plt.ylim(0, 16)
else:
	plt.yticks(np.arange(24, 28.5, 0.5), fontsize=8)
	plt.ylim(24, 28)
plt.grid(linestyle='--')
plt.setp(ax.get_xticklabels(), visible=False)
l1, l2, l3 = annual_cycle
plt.setp(l1, color='blue')
plt.setp(l2, color='red')
plt.setp(l3, color='black')

legend = ['RegCM5', 'ERA5', 'CRU']
plt.legend(annual_cycle, legend, ncol=6, loc=1, fontsize=8)

ax = fig.add_subplot(3, 1, 2)
annual_cycle = ax.plot(time, clim_lpb_rcm, time, clim_lpb_rea, time, clim_lpb_obs)
plt.title(u'(b) LPB', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
if var_rcm == 'pr':
	plt.ylabel('Precipitation (mm d⁻¹)', fontweight='bold', fontsize=8)
	plt.yticks(np.arange(0, 16, 2), fontsize=8)
	plt.ylim(0, 14)
else:
	plt.ylabel('Temperature (°C)', fontweight='bold', fontsize=8)
	plt.yticks(np.arange(13, 28, 2), fontsize=8)
	plt.ylim(13, 27)
plt.grid(linestyle='--')
plt.setp(ax.get_xticklabels(), visible=False)
l1, l2, l3 = annual_cycle
plt.setp(l1, color='blue')
plt.setp(l2, color='red')
plt.setp(l3, color='black')

ax = fig.add_subplot(3, 1, 3)
annual_cycle = ax.plot(time, clim_neb_rcm, time, clim_neb_rea, time, clim_neb_obs)
plt.title(u'(c) NEB', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.xlabel('Months', fontsize=8, fontweight='bold')
if var_rcm == 'pr':
	plt.yticks(np.arange(0, 16, 2), fontsize=8)
	plt.ylim(0, 14)
else:
	plt.yticks(np.arange(24, 28.5, 0.5), fontsize=8)
	plt.ylim(24, 28)
plt.grid(linestyle='--')
l1, l2, l3 = annual_cycle
plt.setp(l1, color='blue')
plt.setp(l2, color='red')
plt.setp(l3, color='black')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_annual_cycle_{0}.png'.format(var_rcm)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


