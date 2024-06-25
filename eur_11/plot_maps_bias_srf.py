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

mbe_djf_noto_eobs = compute_mbe(noto_djf, eobs_djf)
mbe_mam_noto_eobs = compute_mbe(noto_mam, eobs_mam)
mbe_jja_noto_eobs = compute_mbe(noto_jja, eobs_jja)
mbe_son_noto_eobs = compute_mbe(noto_son, eobs_son)

mbe_djf_wdm7_eobs = compute_mbe(wdm7_djf, eobs_djf)
mbe_mam_wdm7_eobs = compute_mbe(wdm7_mam, eobs_mam)
mbe_jja_wdm7_eobs = compute_mbe(wdm7_jja, eobs_jja)
mbe_son_wdm7_eobs = compute_mbe(wdm7_son, eobs_son)

mbe_djf_wsm7_eobs = compute_mbe(wsm7_djf, eobs_djf)
mbe_mam_wsm7_eobs = compute_mbe(wsm7_mam, eobs_mam)
mbe_jja_wsm7_eobs = compute_mbe(wsm7_jja, eobs_jja)
mbe_son_wsm7_eobs = compute_mbe(wsm7_son, eobs_son)

mbe_djf_wsm5_eobs = compute_mbe(wsm5_djf, eobs_djf)
mbe_mam_wsm5_eobs = compute_mbe(wsm5_mam, eobs_mam)
mbe_jja_wsm5_eobs = compute_mbe(wsm5_jja, eobs_jja)
mbe_son_wsm5_eobs = compute_mbe(wsm5_son, eobs_son)

# Plot figure
fig = plt.figure(figsize=(10, 6))
dict_plot = {'pr': ['Bias of precipitation (mm d$^-$$^1$)', np.arange(-10, 11, 1), cm.BrBG]}
font_size = 8
	
ax = fig.add_subplot(4, 4, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_djf_noto_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(a) NoTo-EOBS DJF', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 4, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_djf_wdm7_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(b) WDM7-EOBS DJF', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 4, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_djf_wsm7_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(c) WSM7-EOBS DJF', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 4, 4)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_djf_wsm5_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(d) WSM5-EOBS DJF', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 4, 5)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mam_noto_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(e) NoTo-EOBS MAM', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 4, 6)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mam_wdm7_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(f) WDM7-EOBS MAM', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 4, 7)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mam_wsm7_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(g) WSM7-EOBS MAM', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 4, 8)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mam_wsm5_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(h) WSM5-EOBS MAM', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 4, 9)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_jja_noto_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(i) NoTo-EOBS JJA', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 4, 10)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_jja_wdm7_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(j) WDM7-EOBS JJA', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 4, 11)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_jja_wsm7_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(k) WSM7-EOBS JJA', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 4, 12)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_jja_wsm5_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(l) WSM5-EOBS JJA', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 4, 13)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_son_noto_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(m) NoTo-EOBS SON', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 4, 14)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_son_wdm7_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(n) WDM7-EOBS SON', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 4, 15)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_son_wsm7_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(o) WSM7-EOBS SON', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 4, 16)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_son_wsm5_eobs, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(p) WSM5-EOBS SON', loc='left', fontsize=font_size, fontweight='bold')

# Set colobar
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.3, 0.01, 0.4]))
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/user/mdasilva/EUR-11/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
