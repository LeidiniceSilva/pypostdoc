# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot bias maps"

import os
import netCDF4
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap
from import_climate_tools import compute_mbe

var = 'pr_int'
domain = 'EUR-11'
dt = '20000101'

path = '/marconi/home/userexternal/mdasilva'


def import_obs(param, domain, dataset):

	arq   = '{0}/user/mdasilva/EUR-11/post_evaluate/obs/obs_v2/{1}_freq_{2}_{3}_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	
		
def import_rcm(exp, param, domain, dataset):

	arq   = '{0}/user/mdasilva/EUR-11/post_evaluate/rcm/{1}/{2}_freq_{3}_{4}_{5}_lonlat.nc'.format(path, exp, param, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	
	
def basemap(lat, lon):

	lat_start, lat_end, lon_start, lon_end = 15, 75, -45, 65
	
	map = Basemap(projection='cyl', llcrnrlon=lon_start, llcrnrlat=lat_start, urcrnrlon=lon_end,urcrnrlat=lat_end, resolution='c')
	map.drawmeridians(np.arange(lon_start, lon_end, 20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(lat_start, lat_end, 10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	map.drawcoastlines(linewidth=0.5, color='black')
	
	lons, lats = np.meshgrid(lon, lat)
	xx, yy = map(lons,lats)
	
	return map, xx, yy
	

# Import model and obs dataset
lat, lon, cpc_jan = import_obs('precip', domain, 'CPC')
lat, lon, wdm7_jan_v1 = import_rcm('wdm7-Europe_v1', 'pr', domain, 'RegCM5')
lat, lon, wdm7_jan_v2 = import_rcm('wdm7-Europe_v2', 'pr', domain, 'RegCM5')

mbe_wdm7_jan_v1_cpc = compute_mbe(wdm7_jan_v1, cpc_jan)
mbe_wdm7_jan_v2_cpc = compute_mbe(wdm7_jan_v2, cpc_jan)

# Plot figure
fig = plt.figure()   
font_size = 8
	
levs = np.arange(-35, 38, 3)
legend = 'Intensity (mm d$^-$$^1$)'

ax = fig.add_subplot(1, 2, 1)
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_jan_v1_cpc[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(a) WDM7_v1 - CPC Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(1, 2, 2)
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_jan_v2_cpc[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(b) WDM7_v2 - CPC Jan', loc='left', fontsize=font_size, fontweight='bold')

cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.3, 0.018, 0.4]))
cbar.set_label('Frequency (%)'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/user/mdasilva/EUR-11/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_WDM7_v1-v2_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
