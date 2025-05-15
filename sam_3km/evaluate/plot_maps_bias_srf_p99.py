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

var = 'p99'
freq = 'hourly'
domain = 'SAM-3km'
idt, fdt = '2018', '2021'

if freq == 'hourly':
	dt = '1hr_{0}-{1}'.format(idt, fdt)
	dataset = 'CMORPH'
else:
	dt = '{0}-{1}'.format(idt, fdt)
	dataset = 'CPC'

path = '/leonardo/home/userexternal/mdasilva/leonardo_work'


def import_obs(param, domain, dataset):

	arq   = '{0}/SAM-3km/postproc/evaluate/obs/p99_{1}_{2}_{3}_lonlat.nc'.format(path, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]
	
	return lat, lon, mean
	
		
def import_rcm(param, domain, dataset):

	arq   = '{0}/SAM-3km/postproc/evaluate/rcm/p99_{1}_{2}_{3}_lonlat.nc'.format(path, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]
	
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
lat, lon, regcm = import_rcm('pr', domain, 'RegCM5')

if freq == 'hourly':
	lat, lon, obs = import_obs('cmorph', domain, dataset)
else:
	lat, lon, obs = import_obs('precip', domain, dataset)

mbe_regcm_obs = regcm - obs

# Plot figure
fig, ax = plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()})
font_size = 8
	
if freq == 'hourly':
	levs = np.arange(-5, 5.5, 0.5)
	legend = 'Hourly p99 (mm h$^-$$^1$)'
else:
	levs = np.arange(-60, 70, 10)
	legend = 'Daily p99 (mm d$^-$$^1$)'

plt_map = ax.contourf(lon, lat, mbe_regcm_obs, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
ax.set_title(u'(a) RegCM5-{0}'.format(dataset), loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax)

cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.93, 0.3, 0.03, 0.4]))
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
exit()
