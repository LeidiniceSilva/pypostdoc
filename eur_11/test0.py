# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot clim maps"

import os
import netCDF4
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import maskoceans

var = 'pr'
dt = '200001'
domain = 'EUR-11'
path = '/marconi/home/userexternal/mdasilva'
	
			
def import_obs(param, dataset):

	if param == 'precip':
		param_ = 'precip'
	elif param == 'rr':
		param_ = 'rr'
	else:
		param_ = 'tp'
	
	arq   = '{0}/user/mdasilva/EUR-11/postproc/obs/{1}_{2}_{3}_{4}_lonlat.nc'.format(path, param, domain, dataset, dt) 
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param_][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

	return lat, lon, mean


def import_rcm(param, dataset):

	arq   = '{0}/user/mdasilva/EUR-11/postproc/rcm/{1}_{2}_{3}_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

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
dict_var = {'pr': ['precip', 'rr', 'pr']}

lat, lon, cpc_jan = import_obs(dict_var[var][0], 'CPC')
lat, lon, eobs_jan = import_obs(dict_var[var][1], 'EOBS')
lat, lon, era5_jan = import_obs(dict_var[var][2], 'ERA5')

lat, lon, noto_jan = import_rcm(var, 'NoTo-Europe')
lat, lon, wsm5_jan = import_rcm(var, 'WSM5-Europe')
lat, lon, wsm7_jan = import_rcm(var, 'WSM7-Europe')
lat, lon, wdm7_jan = import_rcm(var, 'WDM7-Europe')

# Plot figure
fig = plt.figure(figsize=(10, 6))

color = ['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']
dict_plot = {'pr': ['Precipitation (mm d$^-$$^1$)', np.arange(0, 18, 1), matplotlib.colors.ListedColormap(color)]}
font_size = 8
	
ax = fig.add_subplot(3, 3, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, cpc_jan, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(a) CPC Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, eobs_jan, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(b) EOBS Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, era5_jan, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(c) ERA5 Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 4)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, noto_jan, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(d) NoTo Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 5)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, wsm5_jan, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(e) WSM5 Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 6)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, wsm7_jan, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(f) WSM7 Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 7)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, wdm7_jan, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(g) WDM7 Jan', loc='left', fontsize=font_size, fontweight='bold')

# Set colobar
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.3, 0.01, 0.4]))
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/user/mdasilva/EUR-11/figs'.format(path)
name_out = 'pyplt_maps_clim_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

