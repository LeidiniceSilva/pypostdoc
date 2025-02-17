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
dt = '2000-2001'
domain = 'EUR-11'
dataset = 'EOBS'
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11'


def import_obs(param, domain, dataset):

	arq   = '{0}/postproc/obs/p99_{1}_{2}_{3}_lonlat.nc'.format(path, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	
		
def import_rcm(param, domain, dataset):

	arq   = '{0}/postproc/rcm/p99_{1}_{2}_{3}_lonlat.nc'.format(path, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean


def configure_subplot(ax):

    ax.set_xticks(np.arange(-40,65,25), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(30,85,15), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.tick_params(axis='x', labelsize=6, labelcolor='black')
    ax.tick_params(axis='y', labelsize=6, labelcolor='black')
    ax.grid(c='k', ls='--', alpha=0.4)
    ax.coastlines()
	
	
# Import model and obs dataset
lat, lon, obs = import_obs('rr', domain, dataset)
lat, lon, noto = import_rcm('pr', domain, 'NoTo-Europe_RegCM5')
lat, lon, wdm7 = import_rcm('pr', domain, 'WDM7-Europe_RegCM5')
lat, lon, wsm7 = import_rcm('pr', domain, 'WSM7-Europe_RegCM5')
lat, lon, wsm5 = import_rcm('pr', domain, 'WSM5-Europe_RegCM5')
	
mbe_noto_obs = compute_mbe(noto, obs)	
mbe_wdm7_obs = compute_mbe(wdm7, obs)	
mbe_wsm7_obs = compute_mbe(wsm7, obs)	
mbe_wsm5_obs = compute_mbe(wsm5, obs)  

# Plot figure
fig, axes = plt.subplots(1, 4, figsize=(9, 4), subplot_kw={'projection': ccrs.PlateCarree()})
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
path_out = '{0}/figs/totc'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

