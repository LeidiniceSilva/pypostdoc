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

var = 'p99'
domain = 'EUR-11'
dt = '20000101'
path = '/home/mda_silv/scratch/EUR-11/postproc'

	
def import_obs(param, dataset):

	arq   = '{0}/obs/p99_{1}_{2}_{3}_lonlat.nc'.format(path, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean


def import_rcm(param, dataset):

	arq   = '{0}/rcm/p99_{1}_{2}_{3}_lonlat.nc'.format(path, domain, dataset, dt)	
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
lat, lon, cpc_jan = import_obs('precip', 'CPC')

lat, lon, noto_jan = import_rcm('pr', 'NoTo-Europe')
lat, lon, wsm5_jan = import_rcm('pr', 'WSM5-Europe')
lat, lon, wsm7_jan = import_rcm('pr', 'WSM7-Europe')
lat, lon, wdm7_jan = import_rcm('pr', 'WDM7-Europe')

mbe_noto_cpc_jan = compute_mbe(noto_jan, cpc_jan)
mbe_wsm5_cpc_jan = compute_mbe(wsm5_jan, cpc_jan)
mbe_wsm7_cpc_jan = compute_mbe(wsm7_jan, cpc_jan)
mbe_wdm7_cpc_jan = compute_mbe(wdm7_jan, cpc_jan)

# Plot figure
fig = plt.figure(figsize=(10, 6))   
font_size = 8
	
levs = np.arange(-70, 75, 5)
legend = 'Daily p99 (mm d$^-$$^1$)'

ax = fig.add_subplot(2, 2, 1)
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_noto_cpc_jan[0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(a) NoTo - CPC Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 2, 2)
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wsm5_cpc_jan[0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(b) WSM5 - CPC Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 2, 3)
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wsm7_cpc_jan[0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(c) WSM7 - CPC Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 2, 4)
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_cpc_jan[0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(d) WDM7 - CPC Jan', loc='left', fontsize=font_size, fontweight='bold')

cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.3, 0.018, 0.4]))
cbar.set_label('Daily p99 (mm d$^-$$^1$)'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
