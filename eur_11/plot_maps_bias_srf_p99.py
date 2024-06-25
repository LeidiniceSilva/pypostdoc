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

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import maskoceans
from import_climate_tools import compute_mbe

var = 'p99'
freq = 'daily'
domain = 'EUR-11'
path = '/marconi/home/userexternal/mdasilva'

if freq == 'hourly':
	dt = '1hr_2000-2000'
else:
	dt = '2000-2000'


def import_obs(param, domain, dataset):

	arq   = '{0}/user/mdasilva/EUR-11/post_evaluate/obs/p99_{1}_{2}_{3}_lonlat.nc'.format(path, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	
		
def import_rcm(param, domain, dataset):

	arq   = '{0}/user/mdasilva/EUR-11/post_evaluate/rcm/{2}/p99_{1}_{2}_RegCM5_{3}_lonlat.nc'.format(path, domain, dataset, dt)	
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
if freq == 'hourly':
	lat, lon, eobs = import_obs('rr', domain, 'EOBS')
	lat, lon, noto = import_rcm('pr', domain, 'Noto-Europe')
	lat, lon, wdm7 = import_rcm('pr', domain, 'wdm7-Europe')
	lat, lon, wsm7 = import_rcm('pr', domain, 'wsm7-Europe')
	lat, lon, wsm5 = import_rcm('pr', domain, 'wsm5-Europe')
	
	mbe_noto_eobs = compute_mbe(noto, eobs)	
	mbe_wdm7_eobs = compute_mbe(wdm7, eobs)	
	mbe_wsm7_eobs = compute_mbe(wsm7, eobs)	
	mbe_wsm5_eobs = compute_mbe(wsm5, eobs)  
else:
	lat, lon, eobs = import_obs('rr', domain, 'EOBS')
	lat, lon, noto = import_rcm('pr', domain, 'Noto-Europe')
	lat, lon, wdm7 = import_rcm('pr', domain, 'wdm7-Europe')
	lat, lon, wsm7 = import_rcm('pr', domain, 'wsm7-Europe')
	lat, lon, wsm5 = import_rcm('pr', domain, 'wsm5-Europe')
	
	mbe_noto_eobs = compute_mbe(noto, eobs)	
	mbe_wdm7_eobs = compute_mbe(wdm7, eobs)	
	mbe_wsm7_eobs = compute_mbe(wsm7, eobs)	
	mbe_wsm5_eobs = compute_mbe(wsm5, eobs)  
	
# Plot figure
fig = plt.figure(figsize=(8, 3))   
font_size = 8

if freq == 'hourly':
	levs = np.arange(-5, 5.5, 0.5)
	legend = 'Hourly p99 (mm h$^-$$^1$)'
else:
	levs = np.arange(-50, 55, 5)
	legend = 'Daily p99 (mm d$^-$$^1$)'

ax = fig.add_subplot(1, 4, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_noto_eobs[0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(a) NoTo-EOBS', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(1, 4, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_eobs[0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(b) WDM7-EOBS', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(1, 4, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wsm7_eobs[0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(c) WSM7-EOBS', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(1, 4, 4)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wsm5_eobs[0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(d) WSM5-EOBS', loc='left', fontsize=font_size, fontweight='bold')
	
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.3, 0.01, 0.4]))
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/user/mdasilva/EUR-11/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
