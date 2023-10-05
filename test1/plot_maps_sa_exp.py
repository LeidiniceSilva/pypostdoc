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


def import_ref(param, exp):

	arq = xr.open_dataset('/home/mda_silv/scratch/test1/posp/' + '{0}_{1}_2000010100_lonlat.nc'.format(param, exp))
	data = arq[param]
	time = data.sel(time=slice('2000-01-01','2000-01-31'))
	var = time.groupby('time.month').mean('time')
	lat = var.lat
	lon = var.lon
	mean = np.nanmean(var.values, axis=0)
	
	return lat, lon, mean
	
	
def import_rcm(param, exp):

	arq = xr.open_dataset('/home/mda_silv/scratch/test1/posp/' + '{0}_SAM-22_Reg5_{1}_mon_2000010100_lonlat.nc'.format(param, exp))
	data = arq[param]
	time = data.sel(time=slice('2000-01-01','2000-01-31'))
	var = time.groupby('time.month').mean('time')
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
	path = '/home/mda_silv/github_projects/shp'
	map.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black')
	
	return map, xx, yy
	
	
# Import model and obs database 
var_rcm = 'rsnscl'

lat, lon, clim_rcm_exp1 = import_rcm(var_rcm, 'exp1')
lat, lon, clim_rcm_exp2 = import_rcm(var_rcm, 'exp2')

diff = clim_rcm_exp2 - clim_rcm_exp1

# Plot figure   
fig = plt.figure()

if var_rcm == 'pr':
	levs0 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16]
	color0 = cm.Blues
	levs1 = [-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3]
	levs2 = [-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10]
	color = cm.BrBG
else:
	levs0 = [175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425]
	color0 = cm.rainbow
	levs1 = [-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3]
	color1 = cm.PuOr

ax = fig.add_subplot(1, 3, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, clim_rcm_exp1, levels=levs0, latlon=True, cmap=color0, extend='max') 
plt.title(u'(a) RegCM5 exp1', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=10, fontsize=8, fontweight='bold')

ax = fig.add_subplot(1, 3, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, clim_rcm_exp2, levels=levs0, latlon=True, cmap=color0, extend='max') 
plt.title(u'(b) RegCM5 exp2', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=10, fontsize=8, fontweight='bold')
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.91, 0.3, 0.02, 0.4]))
cbar.ax.tick_params(labelsize=8)

ax = fig.add_subplot(1, 3, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, diff, levels=levs1, latlon=True, cmap=color1, extend='both') 
plt.title(u'(c) Exp 2 - Exp 1', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=10, fontsize=8, fontweight='bold')
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.98, 0.3, 0.02, 0.4]))
cbar.ax.tick_params(labelsize=8)

# Path out to save figure
path_out = '/home/mda_silv/figs/test1'
name_out = 'pyplt_maps_{0}_exp.png'.format(var_rcm)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


