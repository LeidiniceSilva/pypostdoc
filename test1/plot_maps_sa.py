# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "September 16, 2023"
__description__ = "This script plot annual cycle"

import os
import netCDF4
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from mpl_toolkits.basemap import Basemap


def import_ref(param, ref):

	arq = xr.open_dataset('/home/nice/Downloads/' + '{0}_{1}_mon_2000_lonlat.nc'.format(param, ref))
	data = arq[param]
	time = data.sel(time=slice('2000-01-01','2000-12-31'))
	var = time.groupby('time.year').mean('time')
	lat = var.lat
	lon = var.lon
	mean = np.nanmean(var.values, axis=0)
	
	return lat, lon, mean
	
	
def import_rcm(param):

	arq = xr.open_dataset('/home/nice/Downloads/' + '{0}_SAM-22_Reg5_mon_2000_lonlat.nc'.format(param))
	data = arq[param]
	time = data.sel(time=slice('2000-01-01','2000-12-31'))
	var = time.groupby('time.year').mean('time')
	lat = var.lat
	lon = var.lon
	mean = np.nanmean(var.values, axis=0)
	
	return lat, lon, mean
	
	
def basemap(lat, lon):
	
	aux_lon1 = []
	aux_lon2 = []
	for l in lon:
		if l <= 180:
			aux_lon1.append(l)
		else:
			aux_lon2.append(l-360)
		
	lon = np.array(aux_lon1[::-1] + aux_lon2[::-1])
	new_lat = lat
	new_lon = lon[::-1]

	map = Basemap(projection='cyl', llcrnrlon=-105., llcrnrlat=-57., urcrnrlon=-16.,urcrnrlat=18., resolution='c')
	map.drawmeridians(np.arange(-105., -16., 20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-57., 18., 15.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
		
	lons, lats = np.meshgrid(new_lon, new_lat)
	xx, yy = map(lons,lats)

	# Import shapefile 	
	path = '/home/nice/Documentos/github_projects/shp'
	map.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black')
	
	return map, xx, yy
	
	
# Import model and obs database 
var_obs = 'tmp'
var_rea = 't2m'
var_rcm = 'tas'

lat, lon, clim_obs = import_ref(var_obs, 'cru_ts4.07')
lat, lon, clim_rea = import_ref(var_rea, 'era5')
lat, lon, clim_rcm = import_rcm(var_rcm)

bias_rea = clim_rcm[0] - clim_rea 
bias_obs = clim_rcm[0] - clim_obs 

# Plot figure   
fig = plt.figure(figsize=(8, 5))

if var_rcm == 'pr':
	levs0 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
	color0 = cm.Blues
	levs = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]
	color = cm.BrBG
else:
	levs0 = [14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34]
	color0 = cm.Reds
	levs = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]
	color = cm.bwr

ax = fig.add_subplot(2, 3, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, clim_rcm[0], levels=levs0, latlon=True, cmap=color0, extend='max') 
plt.title(u'(a) RegCM5', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=8, fontweight='bold')

ax = fig.add_subplot(2, 3, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, clim_rea, levels=levs0, latlon=True, cmap=color0, extend='max') 
plt.title(u'(b) ERA5', loc='left', fontsize=8, fontweight='bold')

ax = fig.add_subplot(2, 3, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, clim_obs, levels=levs0, latlon=True, cmap=color0, extend='max') 
plt.title(u'(c) CRU', loc='left', fontsize=8, fontweight='bold')
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.91, 0.3, 0.02, 0.4]))
cbar.ax.tick_params(labelsize=8)

ax = fig.add_subplot(2, 3, 4)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_rea, levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(d) RegCM5 - ERA5', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=8, fontweight='bold')

ax = fig.add_subplot(2, 3, 5)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_obs, levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(e) RegCM5 - CRU', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=8, fontweight='bold')

cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.98, 0.3, 0.02, 0.4]))
cbar.ax.tick_params(labelsize=8)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_{0}.png'.format(var_rcm)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


