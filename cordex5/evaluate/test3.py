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
stats='int'
domain = 'CSAM-3'
dt = '200001'
font_size = 8

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5'

def import_obs(param, domain, dataset):

	arq   = '{0}/postproc/evaluate/test_bin/obs/{1}_{2}_{3}_{4}_day_{5}_th1_lonlat.nc'.format(path, param, stats, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,0,:,:]
	
	return lat, lon, mean


def import_rcm_i(param, domain):

	arq   = '{0}/postproc/evaluate/test_bin/era5-csam-3/{1}_{2}_{3}_day_{4}_th1_lonlat.nc'.format(path, param, stats, domain, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,0,:,:]

	return lat, lon, mean


def import_rcm_ii(param, domain):

	arq   = '{0}/postproc/evaluate/test_bin/era5-csam/{1}_{2}_{3}_day_{4}_th1_lonlat.nc'.format(path, param, stats, domain, dt)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,0,:,:]

	return lat, lon, mean


# Import model and obs dataset
lat, lon, cpc = import_obs('precip', domain, 'CPC')
lat, lon, regcm_i = import_rcm_i('pr', domain)
lat, lon, regcm_ii = import_rcm_ii('pr', domain)

regcm_i_cpc = regcm_i - cpc
regcm_ii_cpc = regcm_ii - cpc
regcm_ii_regcm_i = regcm_ii - regcm_i

# Plot figure   
def configure_subplot(ax):

	lon_min = np.round(np.min(lon), 1)
	lon_max = np.round(np.max(lon), 1)
	lat_min = np.round(np.min(lat), 1)
	lat_max = np.round(np.max(lat), 1)
	ax.set_extent([np.min(lon), np.max(lon), np.min(lat), np.max(lat)], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(lon_min,lon_max,10), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(lat_min,lat_max,5), crs=ccrs.PlateCarree())

	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(font_size)

	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.add_feature(cfeat.BORDERS, linewidth=0.5)
	ax.coastlines(linewidth=0.5)


# Plot figure   
fig, axes = plt.subplots(1, 3, figsize=(12, 5), subplot_kw={'projection': ccrs.PlateCarree()})

if stats == 'freq':
	levs = np.arange(-44, 46, 2)
	legend = 'Frequency (%)'
else:
	levs = np.arange(-22, 23, 1)
	legend = 'Intensity (mm d$^-$$^1$)'


ax1 = axes[0]
plt_map1 = ax1.contourf(lon, lat, regcm_i_cpc, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='both') 
ax1.set_title(u'(a) RegCM5 (old_bin) - CPC', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

ax2 = axes[1]
plt_map2 = ax2.contourf(lon, lat, regcm_ii_cpc, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='both') 
ax2.set_title(u'(b) RegCM5 (new_bin) - CPC', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

ax3 = axes[2] 
plt_map3 = ax3.contourf(lon, lat, regcm_ii_regcm_i, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='both') 
ax3.set_title(u'(c) new_bin - old_bin', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

# Set colobar
cbar1_ax = fig.add_axes([0.92, 0.25, 0.02, 0.5]) 
cbar1 = fig.colorbar(plt_map1, cax=cbar1_ax, orientation='vertical')
cbar1.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar1.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/postproc/evaluate/test_bin/figs'.format(path)
name_out = 'pyplt_maps_clim_{0}_{1}_{2}_RegCM5_{3}.png'.format(var, stats, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

