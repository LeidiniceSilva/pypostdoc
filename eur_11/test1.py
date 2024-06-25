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
from import_climate_tools import compute_mbe

var = 'pr'
dt = '20000101'
domain = 'EUR-11'
path = '/marconi/home/userexternal/mdasilva'
	
			
def import_obs(param, dataset):

	arq   = '{0}/user/mdasilva/EUR-11/post_evaluate/obs/{1}_{2}_{3}_{4}_lonlat.nc'.format(path, param, domain, dataset, dt) 
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

	return lat, lon, mean


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

lat, lon, eobs_jan = import_obs(dict_var[var][0], 'EOBS')
lat, lon, mswep_jan = import_obs(dict_var[var][1], 'MSWEP')
lat, lon, cpc_jan = import_obs(dict_var[var][2], 'CPC')
lat, lon, wdm7_jan_v1 = import_rcm('wdm7-Europe_v1', var, 'RegCM5')
lat, lon, wdm7_jan_v2 = import_rcm('wdm7-Europe_v3', var, 'RegCM5')
lat, lon, wdm7_jan_v3 = import_rcm('wdm7-Europe_v4', var, 'RegCM5')
lat, lon, wdm7_jan_v4 = import_rcm('wdm7-Europe_v5', var, 'RegCM5')

mbe_wdm7_jan_v1_eobs = compute_mbe(wdm7_jan_v1, eobs_jan)
mbe_wdm7_jan_v1_mswep = compute_mbe(wdm7_jan_v1, mswep_jan)
mbe_wdm7_jan_v1_cpc = compute_mbe(wdm7_jan_v1, cpc_jan)

mbe_wdm7_jan_v2_eobs = compute_mbe(wdm7_jan_v2, eobs_jan)
mbe_wdm7_jan_v2_mswep = compute_mbe(wdm7_jan_v2, mswep_jan)
mbe_wdm7_jan_v2_cpc = compute_mbe(wdm7_jan_v2, cpc_jan)

mbe_wdm7_jan_v3_eobs = compute_mbe(wdm7_jan_v3, eobs_jan)
mbe_wdm7_jan_v3_mswep = compute_mbe(wdm7_jan_v3, mswep_jan)
mbe_wdm7_jan_v3_cpc = compute_mbe(wdm7_jan_v3, cpc_jan)

mbe_wdm7_jan_v4_eobs = compute_mbe(wdm7_jan_v4, eobs_jan)
mbe_wdm7_jan_v4_mswep = compute_mbe(wdm7_jan_v4, mswep_jan)
mbe_wdm7_jan_v4_cpc = compute_mbe(wdm7_jan_v4, cpc_jan)

# Plot figure
fig = plt.figure(figsize=(10, 8))
dict_plot = {'pr': ['Bias of precipitation (mm d$^-$$^1$)', np.arange(-10, 11, 1), cm.BrBG]}
font_size = 8

ax = fig.add_subplot(4, 3, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_jan_v1_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(a) WDM7_v1 - EOBS Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_jan_v1_mswep, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(b) WDM7_v1 - MSWEP Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_jan_v1_cpc, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(c) WDM7_v1 - CPC Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 4)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_jan_v2_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(d) WDM7_v2 (fix bugs) - EOBS Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 5)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_jan_v2_mswep, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(e) WDM7_v2 (fix bugs) - MSWEP Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 6)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_jan_v2_cpc, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(f) WDM7_v2 (fix bugs) - CPC Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 7)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_jan_v3_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(g) WDM7_v3 (more ccn) - EOBS Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 8)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_jan_v3_mswep, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(h) WDM7_v3 (more ccn) - MSWEP Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 9)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_jan_v3_cpc, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(i) WDM7_v3 (more ccn) - CPC Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 10)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_jan_v4_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(j) WDM7_v4 (less ccn) - EOBS Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 11)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_jan_v4_mswep, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(k) WDM7_v4 (less ccn) - MSWEP Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 12)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wdm7_jan_v4_cpc, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(l) WDM7_v4 (less ccn) - CPC Jan', loc='left', fontsize=font_size, fontweight='bold')

# Set colobar
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.3, 0.01, 0.4]))
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/user/mdasilva/EUR-11/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_WDM7_v1-v2-v3-v4_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

