# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jan 02, 2024"
__description__ = "This script plot cyclone tracking"

import os
import cmocean
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import maskoceans

var = 'pr'
cyclone = 'cyclone-ii'
path = '/marconi/home/userexternal/mdasilva'

if cyclone == 'cyclone-i':
	datetime = 'Jun 14-22, 2023'
else:
	datetime = 'Jul 12-16, 2023'


def import_era5(param, domain, dataset, freq, dt):

	arq   = '{0}/user/mdasilva/SAM-3km-cyclone/post/obs/era5/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, freq, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]

	return lat, lon, mean
	
	
def import_persiann(param, dataset, freq, dt):

	arq   = '{0}/user/mdasilva/SAM-3km-cyclone/post/obs/persiann/{1}_{2}_{3}_{4}_lonlat.nc'.format(path, param, dataset, freq, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]

	return lat, lon, mean
	
	
def import_regcm5(param, domain, exp, dataset, freq, dt):

	arq   = '{0}/user/mdasilva/SAM-3km-cyclone/post/rcm/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, param, domain, exp, dataset, freq, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]

	return lat, lon, mean


def basemap(lat, lon):

	map = Basemap(projection='cyl', llcrnrlon=-82., llcrnrlat=-50., urcrnrlon=-32.,urcrnrlat=-10., resolution='c')
	map.drawmeridians(np.arange(-82., -32., 10.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-50., -10., 5.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')
	map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black')
	
	lons, lats = np.meshgrid(lon, lat)
	xx, yy = map(lons,lats)
	
	return map, xx, yy


# Import model and obs dataset
lat, lon, era5 = import_era5('tp', 'SAM-3km-cyclone', 'ERA5', '1hr', cyclone)
era5_sum = np.sum(era5, axis=0)

lat, lon, persiann = import_persiann('prec', 'persiann', '1hr', cyclone)
persiann_sum = np.sum(persiann, axis=0)

lat, lon, regcm = import_regcm5('pr', 'SAM-3km-cyclone', 'ECMWF-ERA5_evaluation_r1i1p1f1', 'ICTP-RegCM5', '3hr', cyclone)
regcm_sum = np.sum(regcm, axis=0)

# Plot figure
fig = plt.figure(figsize=(12, 4))
font_size = 8
colors=["#ffffffff","#d7f0fcff","#ade0f7ff","#86c4ebff","#60a5d6ff","#4794b3ff","#49a67cff","#55b848ff","#9ecf51ff","#ebe359ff","#f7be4aff","#f58433ff","#ed5a28ff","#de3728ff","#cc1f27ff","#b01a1fff","#911419ff"]

ax = fig.add_subplot(1, 3, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, era5_sum*1000, levels=np.arange(0, 200, 10), cmap=matplotlib.colors.ListedColormap(colors), extend='max')
plt.title(u'(a) ERA5 {0}'.format(datetime), loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=20, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=30, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(1, 3, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, persiann_sum, levels=np.arange(0, 200, 10), cmap=matplotlib.colors.ListedColormap(colors), extend='max') 
plt.title(u'(b) PERSIANN {0}'.format(datetime), loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=20, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(1, 3, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, regcm_sum*10800, levels=np.arange(0, 200, 10), cmap=matplotlib.colors.ListedColormap(colors), extend='max') 
plt.title(u'(c) CP-RegCM5 {0}'.format(datetime), loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=20, fontsize=font_size, fontweight='bold')
cbar = map.colorbar(plt_map, ax=ax)
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km-cyclone/figs'.format(path)
name_out = 'pyplt_maps_tracking_cyclone_{0}_{1}.png'.format(var, cyclone)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
