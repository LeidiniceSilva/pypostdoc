# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
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
dt = '2000-2009'
domain = 'EUR-11'
dataset = 'EOBS'
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11'


def import_obs(param, dataset):

	arq   = '{0}/postproc/obs/p99_{1}_{2}_{3}_lonlat.nc'.format(path, dataset, domain, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	
		
def import_rcm(param, dataset):

	arq   = '{0}/postproc/rcm/p99_{1}_{2}_lonlat.nc'.format(path, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean


# Import model and obs dataset
lat, lon, obs = import_obs('rr', dataset)
lat, lon, noto = import_rcm('pr', 'RegCM5_NoTo-EUR')
lat, lon, wdm7 = import_rcm('pr', 'RegCM5_WDM7-EUR')
lat, lon, wsm7 = import_rcm('pr', 'RegCM5_WSM7-EUR')
lat, lon, wsm5 = import_rcm('pr', 'RegCM5_WSM5-EUR')
	
mbe_noto_obs = compute_mbe(noto, obs)	
mbe_wdm7_obs = compute_mbe(wdm7, obs)	
mbe_wsm7_obs = compute_mbe(wsm7, obs)	
mbe_wsm5_obs = compute_mbe(wsm5, obs)  

# Plot figure   
def configure_subplot(ax):

	lon_min = np.round(np.min(lon), 1)
	lon_max = np.round(np.max(lon), 1)
	lat_min = np.round(np.min(lat), 1)
	lat_max = np.round(np.max(lat), 1)
	ax.set_extent([np.min(lon), np.max(lon), np.min(lat), np.max(lat)], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(lon_min,lon_max,20), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(lat_min,lat_max,15), crs=ccrs.PlateCarree())

	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(6)

	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.coastlines(linewidth=0.6)

fig, axes = plt.subplots(1, 4, figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()

dict_plot = {'p99': ['Daily p99 (mm d$^-$$^1$)', np.arange(-50, 55, 5), cm.BrBG]}
font_size = 8

plot_data = {'Plot 1': {'data': mbe_noto_obs[0], 'title': '(a) NoTo-{0}'.format(dataset)},
'Plot 2': {'data': mbe_wdm7_obs[0], 'title': '(b) WDM7-{0}'.format(dataset)},
'Plot 3': {'data': mbe_wsm7_obs[0], 'title': '(c) WSM7-{0}'.format(dataset)},
'Plot 4': {'data': mbe_wsm5_obs[0], 'title': '(d) WSM5-{0}'.format(dataset)}}

for ax, (key, value) in zip(axes, plot_data.items()):
    data = value['data']
    title = value['title']
    
    contour = ax.contourf(lon, lat, data, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
    ax.set_title(title, loc='left', fontsize=font_size, fontweight='bold')
    configure_subplot(ax)

# Set colobar
cbar = fig.colorbar(contour, ax=fig.axes, orientation='horizontal', pad=0.1, aspect=50)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
exit()

