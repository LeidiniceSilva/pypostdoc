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
font_size = 8

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5'

def import_obs(param, domain, dataset):

	arq   = '{0}/postproc/evaluate/test_bin/obs/{1}_{2}_{3}_mon_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]
	
	return lat, lon, mean


def import_rcm_i(param, domain):

	arq   = '{0}/postproc/evaluate/test_bin/era5-csam-3/{1}_{2}_mon_{3}_lonlat.nc'.format(path, param, domain, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

	return lat, lon, mean


def import_rcm_ii(param, domain):

	arq   = '{0}/postproc/evaluate/test_bin/era5-csam/{1}_{2}_mon_{3}_lonlat.nc'.format(path, param, domain, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

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
fig, axes = plt.subplots(2, 3, figsize=(12, 5), subplot_kw={'projection': ccrs.PlateCarree()})
pr_colormap=['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']

ax1 = axes[0, 0]
plt_map1 = ax1.contourf(lon, lat, cpc, transform=ccrs.PlateCarree(), levels=np.arange(0, 18, 1), cmap=matplotlib.colors.ListedColormap(pr_colormap), extend='max') 
ax1.set_title(u'(a) CPC Jan2000', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

ax2 = axes[0, 1]
plt_map2 = ax2.contourf(lon, lat, regcm_i, transform=ccrs.PlateCarree(), levels=np.arange(0, 18, 1), cmap=matplotlib.colors.ListedColormap(pr_colormap), extend='max') 
ax2.set_title(u'(b) RegCM5 (old_bin) Jan2000', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

ax3 = axes[0, 2] 
plt_map3 = ax3.contourf(lon, lat, regcm_ii, transform=ccrs.PlateCarree(), levels=np.arange(0, 18, 1), cmap=matplotlib.colors.ListedColormap(pr_colormap), extend='max') 
ax3.set_title(u'(c) RegCM5 (new_bin) Jan2000', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

ax4 = axes[1, 0]
plt_map4 = ax4.contourf(lon, lat, regcm_i_cpc, transform=ccrs.PlateCarree(), levels=np.arange(-10, 11, 1), cmap=cm.BrBG, extend='both') 
ax4.set_title(u'(d) RegCM5 (old_bin) - CPC', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax4)

ax5 = axes[1, 1]
plt_map5 = ax5.contourf(lon, lat, regcm_ii_cpc, transform=ccrs.PlateCarree(), levels=np.arange(-10, 11, 1), cmap=cm.BrBG, extend='both') 
ax5.set_title(u'(e) RegCM5 (new_bin) - CPC', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax5)

ax6 = axes[1, 2]
plt_map6 = ax6.contourf(lon, lat, regcm_ii_regcm_i, transform=ccrs.PlateCarree(), levels=np.arange(-10, 11, 1), cmap=cm.BrBG, extend='both') 
ax6.set_title(u'(f) new_bin - old_bin', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax6)

# Set colobar
cbar1_ax = fig.add_axes([0.92, 0.55, 0.02, 0.35]) 
cbar1 = fig.colorbar(plt_map1, cax=cbar1_ax, orientation='vertical')
cbar1.set_label('Precipitation (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
cbar1.ax.tick_params(labelsize=font_size)

cbar2_ax = fig.add_axes([0.92, 0.1, 0.02, 0.35])  
cbar2 = fig.colorbar(plt_map4, cax=cbar2_ax, orientation='vertical')
cbar2.set_label('Bias of  precipitation (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
cbar2.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/postproc/evaluate/test_bin/figs'.format(path)
name_out = 'pyplt_maps_clim_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
