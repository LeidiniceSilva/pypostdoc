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

var = 'pr'
dt = '20000101'
domain = 'EUR-11'
path = '/marconi/home/userexternal/mdasilva'


def import_rcm(exp, param, dataset):

	arq   = '{0}/user/mdasilva/EUR-11/post_evaluate/rcm/{1}/{2}_{3}_{4}_{5}_lonlat.nc'.format(path, exp, param, domain, dataset, dt)	
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
dict_var = {'pr': ['rr', 'precipitation', 'precip', ]}

lat, lon, wdm7_jan_v1 = import_rcm('wdm7-Europe_v1', var, 'RegCM5')
lat, lon, wdm7_jan_v2 = import_rcm('wdm7-Europe_v2', var, 'RegCM5')
lat, lon, wdm7_jan_v3 = import_rcm('wdm7-Europe_v3', var, 'RegCM5')
lat, lon, wdm7_jan_v4 = import_rcm('wdm7-Europe_v4', var, 'RegCM5')

mbe_wdm7_jan_v4_v1 = wdm7_jan_v4 - wdm7_jan_v1
mbe_wdm7_jan_v4_v2 = wdm7_jan_v4 - wdm7_jan_v2
mbe_wdm7_jan_v4_v3 = wdm7_jan_v4 - wdm7_jan_v3

# Plot figure
fig = plt.figure(figsize=(10, 6))
dict_plot = {'pr': ['Diff of precipitation (mm d$^-$$^1$)', np.arange(-4, 4.5, 0.5), cm.PiYG]}
font_size = 8

ax = fig.add_subplot(2, 2, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_jan_v4_v1, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(a) WDM7_v4(less ccn2) - WDM7_v1(ctrl) Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 2, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_jan_v4_v2, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(a) WDM7_v4(less ccn2) - WDM7_v2(more ccn) Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 2, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_jan_v4_v3, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(a) WDM7_v4(less ccn2) - WDM7_v3(less ccn) Jan', loc='left', fontsize=font_size, fontweight='bold')

# Set colobar
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.3, 0.01, 0.4]))
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/user/mdasilva/EUR-11/figs'.format(path)
name_out = 'pyplt_maps_diff_{0}_{1}_RegCM5_WDM7_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

