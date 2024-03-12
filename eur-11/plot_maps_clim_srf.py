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
domain = 'EUR-11'
dt = '2000-2000'
path = '/marconi/home/userexternal/mdasilva'
	
			
def import_obs(param, dataset, season):

	arq   = '{0}/user/mdasilva/EUR-11/post_evaluate/obs/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

	return lat, lon, mean


def import_rcm(param, dataset, season):

	arq   = '{0}/user/mdasilva/EUR-11/post_evaluate/rcm/{1}/{2}_{3}_{4}_RegCM5_{5}_{6}_lonlat.nc'.format(path, dataset, param, domain, dataset, season, dt)	
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

lat, lon, eobs_djf = import_obs(dict_var[var][0], 'EOBS', 'DJF')
lat, lon, eobs_mam = import_obs(dict_var[var][0], 'EOBS', 'MAM')
lat, lon, eobs_jja = import_obs(dict_var[var][0], 'EOBS', 'JJA')
lat, lon, eobs_son = import_obs(dict_var[var][0], 'EOBS', 'SON')

lat, lon, mswep_djf = import_obs(dict_var[var][1], 'MSWEP', 'DJF')
lat, lon, mswep_mam = import_obs(dict_var[var][1], 'MSWEP', 'MAM')
lat, lon, mswep_jja = import_obs(dict_var[var][1], 'MSWEP', 'JJA')
lat, lon, mswep_son = import_obs(dict_var[var][1], 'MSWEP', 'SON')

lat, lon, cpc_djf = import_obs(dict_var[var][2], 'CPC', 'DJF')
lat, lon, cpc_mam = import_obs(dict_var[var][2], 'CPC', 'MAM')
lat, lon, cpc_jja = import_obs(dict_var[var][2], 'CPC', 'JJA')
lat, lon, cpc_son = import_obs(dict_var[var][2], 'CPC', 'SON')

lat, lon, noto_djf = import_rcm(var, 'Noto-Europe', 'DJF')
lat, lon, noto_mam = import_rcm(var, 'Noto-Europe', 'MAM')
lat, lon, noto_jja = import_rcm(var, 'Noto-Europe', 'JJA')
lat, lon, noto_son = import_rcm(var, 'Noto-Europe', 'SON')

lat, lon, wdm7_djf = import_rcm(var, 'wdm7-Europe', 'DJF')
lat, lon, wdm7_mam = import_rcm(var, 'wdm7-Europe', 'MAM')
lat, lon, wdm7_jja = import_rcm(var, 'wdm7-Europe', 'JJA')
lat, lon, wdm7_son = import_rcm(var, 'wdm7-Europe', 'SON')

lat, lon, wsm7_djf = import_rcm(var, 'wsm7-Europe', 'DJF')
lat, lon, wsm7_mam = import_rcm(var, 'wsm7-Europe', 'MAM')
lat, lon, wsm7_jja = import_rcm(var, 'wsm7-Europe', 'JJA')
lat, lon, wsm7_son = import_rcm(var, 'wsm7-Europe', 'SON')

lat, lon, wsm5_djf = import_rcm(var, 'wsm5-Europe', 'DJF')
lat, lon, wsm5_mam = import_rcm(var, 'wsm5-Europe', 'MAM')
lat, lon, wsm5_jja = import_rcm(var, 'wsm5-Europe', 'JJA')
lat, lon, wsm5_son = import_rcm(var, 'wsm5-Europe', 'SON')


# Plot figure
color = ['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']
dict_plot = {'pr': ['Precipitation (mm d$^-$$^1$)', np.arange(0, 18, 1), matplotlib.colors.ListedColormap(color)]}
font_size = 8
	
fig = plt.figure(figsize=(12, 6))

ax = fig.add_subplot(4, 5, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, eobs_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(a) EOBS DJF', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, noto_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(b) NoTo DJF', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, wdm7_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(c) WDM7 DJF', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 4)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, wsm7_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(d) WSM7 DJF', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 5)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, wsm5_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(e) WSM5 DJF', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 6)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, eobs_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(f) EOBS MAM', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 7)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, noto_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(g) NoTo MAM', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 8)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, wdm7_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(h) WDM7 MAM', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 9)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, wsm7_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(i) WSM7 MAM', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 10)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, wsm5_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(j) WSM5 MAM', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 11)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, eobs_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(k) EOBS JJA', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 12)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, noto_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(l) NoTo JJA', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 13)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, wdm7_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(m) WDM7 JJA', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 14)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, wsm7_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(n) WSM7 JJA', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 15)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, wsm5_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(o) WSM5 JJA', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 16)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, eobs_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(p) EOBS SON', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 17)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, noto_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(q) NoTo SON', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 18)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, wdm7_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(r) WDM7 SON', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 19)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, wsm7_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(s) WSM7 SON', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 20)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, wsm5_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(t) WSM5 SON', loc='left', fontsize=font_size, fontweight='bold')

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
