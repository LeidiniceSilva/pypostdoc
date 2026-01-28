# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot bias maps"

import os
import netCDF4
import argparse
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat

from cartopy import config
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from import_climate_tools import compute_mbe

parser = argparse.ArgumentParser()
parser.add_argument('--var', required=True, help='Variable name')
parser.add_argument('--freq', required=True, help='Frequency')
parser.add_argument('--domain', required=True, help='Domain name')
parser.add_argument('--idt', required=True, help='Initial year')
parser.add_argument('--fdt', required=True, help='Final year')
args = parser.parse_args()

var = args.var
freq = args.freq
domain = args.domain
idt = args.idt
fdt = args.fdt
font_size = 6

if freq == 'hourly':
	dt = '1hr_{0}-{1}'.format(idt, fdt)
else:
	dt = '{0}-{1}'.format(idt, fdt)

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5'


def import_obs(param, domain, dataset):
		
	arq = '{0}/postproc/evaluate/obs/p99_{1}_{2}_{3}_lonlat.nc'.format(path, domain, dataset, dt)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]
	
	return lat, lon, mean
	
		
def import_rcm(param, domain, dataset):
		
	arq = '{0}/postproc/evaluate/rcm/p99_{1}_{2}_{3}_lonlat.nc'.format(path, domain, dataset, dt)		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]
	
	return lat, lon, mean


# Import model and obs dataset
if freq == 'hourly':
	lat, lon, cmorph = import_obs('cmorph', domain, 'CMORPH')
	lat, lon, era5   = import_obs('tp', domain, 'ERA5')
	lat, lon, regcm  = import_rcm('pr', domain, 'RegCM5')

	mbe_regcm_cmorph = regcm - cmorph
	mbe_regcm_era5   = regcm - era5
else:
	lat, lon, cpc    = import_obs('precip', domain, 'CPC')
	lat, lon, cmorph = import_obs('cmorph', domain, 'CMORPH')
	lat, lon, mswep  = import_obs('precipitation', domain, 'MSWEP')
	lat, lon, era5   = import_obs('tp', domain, 'ERA5')
	lat, lon, regcm  = import_rcm('pr', domain, 'RegCM5')
	
	mbe_regcm_cpc    = regcm - cpc
	mbe_regcm_cmorph = regcm - cmorph
	mbe_regcm_mswep  = regcm - mswep
	mbe_regcm_era5   = regcm - era5

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
		label.set_fontsize(6)

	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.add_feature(cfeat.BORDERS, linewidth=0.5)
	ax.coastlines(linewidth=0.5)
		
if freq == 'hourly':
	fig, axes = plt.subplots(1, 2, figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree()})
	levs = np.arange(-5, 5.5, 0.5)
	legend = 'Hourly p99 (mm h$^-$$^1$)'

	ax1 = axes[0]
	plt_map = ax1.contourf(lon, lat, mbe_regcm_cmorph, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax1.set_title(u'(a) RegCM5-CMORPH', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax1)

	ax2 = axes[1]
	plt_map = ax2.contourf(lon, lat, mbe_regcm_era5, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax2.set_title(u'(b) RegCM5-ERA5', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax2)
else:
	fig, axes = plt.subplots(1, 4, figsize=(12, 4), subplot_kw={'projection': ccrs.PlateCarree()})
	levs = np.arange(-60, 70, 10)
	legend = 'Daily p99 (mm d$^-$$^1$)'

	ax1 = axes[0]
	plt_map = ax1.contourf(lon, lat, mbe_regcm_cpc, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax1.set_title(u'(a) RegCM5-CPC', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax1)

	ax2 = axes[1]
	plt_map = ax2.contourf(lon, lat, mbe_regcm_cmorph, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax2.set_title(u'(b) RegCM5-CMORPH', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax2)

	ax3 = axes[2]
	plt_map = ax3.contourf(lon, lat, mbe_regcm_mswep, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax3.set_title(u'(c) RegCM5-MSWEP', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax3)

	ax4 = axes[3]
	plt_map = ax4.contourf(lon, lat, mbe_regcm_era5, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax4.set_title(u'(d) RegCM5-ERA5', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax4)

# Set colobar
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.91, 0.3, 0.01, 0.4]))
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/figs/evaluate'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
#plt.show()
exit()
