# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot bias maps"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap
from import_climate_tools import compute_mbe

var = 'p99'
domain = 'SAM-3km'
dt = '2018-2021'
path = '/marconi/home/userexternal/mdasilva'


def import_obs(param, domain, dataset, season):

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/obs/p99_{1}_{2}_{3}_lonlat.nc'.format(path, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	
		
def import_cp_3km(param, domain, dataset, season):

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/p99_{1}_{2}_{3}_lonlat.nc'.format(path, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean


def basemap(lat, lon):
	
	map = Basemap(projection='cyl', llcrnrlon=-80., llcrnrlat=-38., urcrnrlon=-34.,urcrnrlat=-8., resolution='c')
	map.drawmeridians(np.arange(-80., -34., 12.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-38., -8., 6.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black')
	
	lons, lats = np.meshgrid(lon, lat)
	xx, yy = map(lons,lats)
	
	return map, xx, yy
	
	
# Import model and obs dataset
lat, lon, cpc = import_obs('precip', domain, 'CPC', 'DJF')

lat, lon, era5 = import_obs('tp', domain, 'ERA5', 'DJF')

lat, lon, regcm = import_cp_3km('pr', domain, 'RegCM5', 'DJF')

# Plot figure
fig = plt.figure()   
lev=np.arange(0, 100, 5)
cor=cm.Blues
font_size = 8
	
ax = fig.add_subplot(1, 3, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, cpc[0], levels=lev, cmap=cor, extend='both') 
plt.title(u'(a) CPC', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(1, 3, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, era5[0], levels=lev, cmap=cor, extend='both') 
plt.title(u'(a) ERA5', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(1, 3, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, regcm[0], levels=lev, cmap=cor, extend='both') 
plt.title(u'(c) RegCM5', loc='left', fontsize=font_size, fontweight='bold')

cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.84, 0.3, 0.03, 0.4]))
cbar.set_label('Extreme wet precipitation (mm $^-$$^1$)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_maps_bias_{0}_CP-RegCM5_{1}_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
