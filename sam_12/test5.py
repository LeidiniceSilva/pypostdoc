# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot bias maps"

import os
import netCDF4
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat

from cartopy import config
from matplotlib.colors import LinearSegmentedColormap
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

var = 'pr'
domain = 'CSAM-3'
dt = '200001'
font_size = 10

subdomains = {'LPB': {'lat': (-31.7, -21.7), 'lon': (-63.8, -53.8)}}

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5'

def import_obs(param, domain, dataset):

	arq   = '{0}/postproc/evaluate/test_bin/obs/{1}_{2}_{3}_day_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	bounds = subdomains['LPB']
	lat_inds = np.where((lat >= bounds['lat'][0]) & (lat <= bounds['lat'][1]))[0]
	lon_inds = np.where((lon >= bounds['lon'][0]) & (lon <= bounds['lon'][1]))[0]
	value_ = value[:, lat_inds[:, None], lon_inds]

	mean = np.nanmean(np.nanmean(value_, axis=1), axis=1)

	return mean


def import_rcm_i(param, domain):

	arq   = '{0}/postproc/evaluate/test_bin/era5-csam-3/{1}_{2}_day_{3}_lonlat.nc'.format(path, param, domain, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	bounds = subdomains['LPB']
	lat_inds = np.where((lat >= bounds['lat'][0]) & (lat <= bounds['lat'][1]))[0]
	lon_inds = np.where((lon >= bounds['lon'][0]) & (lon <= bounds['lon'][1]))[0]
	value_ = value[:, lat_inds[:, None], lon_inds]

	mean = np.nanmean(np.nanmean(value_, axis=1), axis=1)

	return mean


def import_rcm_ii(param, domain):

	arq   = '{0}/postproc/evaluate/test_bin/era5-csam/{1}_{2}_day_{3}_lonlat.nc'.format(path, param, domain, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	bounds = subdomains['LPB']
	lat_inds = np.where((lat >= bounds['lat'][0]) & (lat <= bounds['lat'][1]))[0]
	lon_inds = np.where((lon >= bounds['lon'][0]) & (lon <= bounds['lon'][1]))[0]
	value_ = value[:, lat_inds[:, None], lon_inds]

	mean = np.nanmean(np.nanmean(value_, axis=1), axis=1)

	return mean


# Import model and obs dataset
mean_cpc = import_obs('precip', domain, 'CPC')
mean_regcm_i = import_rcm_i('pr', domain)
mean_regcm_ii = import_rcm_ii('pr', domain)

# Plot figure   
fig, ax = plt.subplots(figsize=(12, 6))

ax.plot(mean_cpc, color='black', label='CPC', linewidth=1., marker='.', linestyle='-')
ax.plot(mean_regcm_i[:-1], color='red', label='RegCM5 (old_bin)', linewidth=1., marker='.', linestyle='-')
ax.plot(mean_regcm_ii[:-1], color='blue', label='RegCM5 (new_bin)', linewidth=1., marker='.', linestyle='-')
ax.set_title('(a) LPB (Jan/2000)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_ylabel('Precipitation (mm d$^{-1}$)', fontsize=font_size, fontweight='bold')
ax.set_xlabel('Daily timeseries', fontsize=font_size, fontweight='bold')
ax.grid(True, linestyle='--')
plt.legend(loc=1, ncol=1, fontsize=font_size)

# Path out to save figure
path_out = '{0}/postproc/evaluate/test_bin/figs'.format(path)
name_out = 'pyplt_graph_timeseries_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
