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
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from import_climate_tools import compute_mbe

var = 'pr'
stats = 'int'
freq = 'hourly'
domain = 'SAM-3km'
idt, fdt = '2018', '2021'

if freq == 'hourly':
	dataset = 'CMORPH'
else:
	dataset = 'CPC'

path = '/leonardo/home/userexternal/mdasilva/leonardo_work'


def import_obs(param, domain, dataset, season):

	if freq == 'hourly':
		dt = '1hr_{0}_{1}-{2}_th0.5'.format(season, idt, fdt)
	else:
		dt = '{0}_{1}-{2}'.format(season, idt, fdt)

	arq   = '{0}/SAM-3km/postproc/evaluate/obs/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, stats, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,0,:,:]
	
	return lat, lon, mean
	
		
def import_rcm(param, domain, dataset, season):

	if freq == 'hourly':
		dt = '1hr_{0}_{1}-{2}_th0.5'.format(season, idt, fdt)
	else:
		dt = '{0}_{1}-{2}'.format(season, idt, fdt)

	arq   = '{0}/SAM-3km/postproc/evaluate/rcm/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, stats, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,0,:,:]
	
	return lat, lon, mean
	

def configure_subplot(ax):

	ax.set_extent([-80, -34, -38, -8], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(-80,-34,12), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(-38,-8,6), crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.tick_params(labelsize=font_size)
	ax.add_feature(cfeat.BORDERS)
	ax.coastlines()	
	ax.grid(c='k', ls='--', alpha=0.4)
	
	
# Import model and obs dataset
lat, lon, regcm_djf = import_rcm('pr', domain, 'RegCM5', 'DJF')
lat, lon, regcm_mam = import_rcm('pr', domain, 'RegCM5', 'MAM')
lat, lon, regcm_jja = import_rcm('pr', domain, 'RegCM5', 'JJA')
lat, lon, regcm_son = import_rcm('pr', domain, 'RegCM5', 'SON')

if freq == 'hourly':
	lat, lon, obs_djf = import_obs('cmorph', domain, dataset, 'DJF')
	lat, lon, obs_mam = import_obs('cmorph', domain, dataset, 'MAM')
	lat, lon, obs_jja = import_obs('cmorph', domain, dataset, 'JJA')
	lat, lon, obs_son = import_obs('cmorph', domain, dataset, 'SON')
else:
	lat, lon, obs_djf = import_obs('precip', domain, dataset, 'DJF')
	lat, lon, obs_mam = import_obs('precip', domain, dataset, 'MAM')
	lat, lon, obs_jja = import_obs('precip', domain, dataset, 'JJA')
	lat, lon, obs_son = import_obs('precip', domain, dataset, 'SON')

mbe_regcm_obs_djf = regcm_djf - obs_djf
mbe_regcm_obs_mam = regcm_mam - obs_mam
mbe_regcm_obs_jja = regcm_jja - obs_jja
mbe_regcm_obs_son = regcm_son - obs_son

# Plot figure
fig, axes = plt.subplots(4, 1, figsize=(5, 12), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()
font_size = 8

if stats == 'freq':
	if freq == 'hourly':
		levs = np.arange(-22, 23, 1)
		legend = 'Frequency (%)'
		dt = '1hr_{0}-{1}_th0.5'.format(idt, fdt)
	else:
		levs = np.arange(-44, 46, 2)
		legend = 'Frequency (%)'
		dt = '{0}-{1}'.format(idt, fdt)
else:
	if freq == 'hourly':
		levs = np.arange(-6, 6.5, 0.5)
		legend = 'Intensity (mm h$^-$$^1$)'
		dt = '1hr_{0}-{1}_th0.5'.format(idt, fdt)
	else:
		levs = np.arange(-22, 23, 1)
		legend = 'Intensity (mm d$^-$$^1$)'
		dt = '{0}-{1}'.format(idt, fdt)

ax1 = axes[0]
plt_map = ax1.contourf(lon, lat, mbe_regcm_obs_djf, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
ax1.set_title(u'(a) RegCM5-{0}'.format(dataset), loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

ax2 = axes[1]
plt_map = ax2.contourf(lon, lat, mbe_regcm_obs_mam, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
ax2.set_title(u'(b) RegCM5-{0}'.format(dataset), loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

ax3 = axes[2]
plt_map = ax3.contourf(lon, lat, mbe_regcm_obs_jja, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
ax3.set_title(u'(c) RegCM5-{0}'.format(dataset), loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

ax4 = axes[3]
plt_map = ax4.contourf(lon, lat, mbe_regcm_obs_son, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
ax4.set_title(u'(d) RegCM5-{0}'.format(dataset), loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax4)

cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.93, 0.3, 0.03, 0.4]))
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_{2}_RegCM5_{3}.png'.format(var, stats, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
exit()
