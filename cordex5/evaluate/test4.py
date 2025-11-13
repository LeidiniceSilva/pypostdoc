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

	step=1
	rain_min=0.1
	rain_max=500.0
	bins = np.arange(rain_min, rain_max + step, step)
	rain_hist = np.zeros(len(bins) - 1)
	for t in range(value_.shape[0]):
		values = value_[t, :, :]
		wet = values[values >= rain_min]
		hist, _ = np.histogram(wet, bins=bins)
		rain_hist += hist

	data.close()

	total = np.sum(rain_hist)
	rain_hist[rain_hist < 1] = np.nan  
	pdf = rain_hist / (total * step)

	bin_centers = (bins[:-1] + bins[1:]) / 2
	
	return pdf


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

	step=1
	rain_min=0.1
	rain_max=500.0
	bins = np.arange(rain_min, rain_max + step, step)
	rain_hist = np.zeros(len(bins) - 1)
	for t in range(value_.shape[0]):
		values = value_[t, :, :]
		wet = values[values >= rain_min]
		hist, _ = np.histogram(wet, bins=bins)
		rain_hist += hist

	data.close()

	total = np.sum(rain_hist)
	rain_hist[rain_hist < 1] = np.nan  
	pdf = rain_hist / (total * step)

	bin_centers = (bins[:-1] + bins[1:]) / 2

	return pdf


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

	step=1
	rain_min=0.1
	rain_max=500.0
	bins = np.arange(rain_min, rain_max + step, step)
	rain_hist = np.zeros(len(bins) - 1)
	for t in range(value_.shape[0]):
		values = value_[t, :, :]
		wet = values[values >= rain_min]
		hist, _ = np.histogram(wet, bins=bins)
		rain_hist += hist

	data.close()

	total = np.sum(rain_hist)
	rain_hist[rain_hist < 1] = np.nan  
	pdf = rain_hist / (total * step)

	bin_centers = (bins[:-1] + bins[1:]) / 2

	return pdf


# Import model and obs dataset
pdf_cpc = import_obs('precip', domain, 'CPC')
pdf_regcm_i = import_rcm_i('pr', domain)
pdf_regcm_ii = import_rcm_ii('pr', domain)

# Plot figure   
fig, ax = plt.subplots()

ax.plot(pdf_cpc, marker='o', markersize=4, mfc='black', mec='black', alpha=0.75, linestyle='None', label='CPC')
ax.plot(pdf_regcm_i, marker='o', markersize=4, mfc='red', mec='red', alpha=0.75, linestyle='None', label='RegCM5 (old_bin)')
ax.plot(pdf_regcm_ii, marker='o', markersize=4, mfc='blue', mec='b', alpha=0.75, linestyle='None', label='RegCM5 (new_bin)')
ax.set_title('(a) LPB (Jan/2000)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_ylabel('Frequency (Log)', fontsize=font_size, fontweight='bold')
ax.set_xlabel('Precipitation (mm d$^{-1}$)', fontsize=font_size, fontweight='bold')
ax.set_yscale('log')
ax.grid(True, linestyle='--')
plt.legend(loc=1, ncol=1, fontsize=font_size)

# Path out to save figure
path_out = '{0}/postproc/evaluate/test_bin/figs'.format(path)
name_out = 'pyplt_graph_pdf_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
