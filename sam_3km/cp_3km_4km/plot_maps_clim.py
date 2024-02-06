# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot bias maps"

import os
import cmocean
import netCDF4
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import maskoceans

var = 'pr'
dt = '2018-2021'
path = '/marconi/home/userexternal/mdasilva'

	
def import_grid(param, domain, dataset, season, dt):

	arq   = '{0}/user/mdasilva/{1}/post_evaluate/{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, domain, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

	return lat, lon, mean


def basemap(lat, lon):

	map = Basemap(projection='cyl', llcrnrlon=-76., llcrnrlat=-35., urcrnrlon=-38.,urcrnrlat=-12., resolution='c')
	map.drawmeridians(np.arange(-76., -38., 10.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-35., -12., 5.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black')
	
	lons, lats = np.meshgrid(lon, lat)
	xx, yy = map(lons,lats)
	
	return map, xx, yy
	
	
# Import model and obs dataset
dict_var = {
'pr': ['pre', 'precip', 'tp'],
'tas': ['tmp', 't2m'],
'tasmax': ['tmx', 'tmax'],
'tasmin': ['tmn', 'tmin']
}

if var == 'pr':
	lat, lon, cru_djf = import_grid(dict_var[var][0], 'SAM-3km', 'CRU', 'DJF', dt)
	lat, lon, cru_mam = import_grid(dict_var[var][0], 'SAM-3km', 'CRU', 'MAM', dt)
	lat, lon, cru_jja = import_grid(dict_var[var][0], 'SAM-3km', 'CRU', 'JJA', dt)
	lat, lon, cru_son = import_grid(dict_var[var][0], 'SAM-3km', 'CRU', 'SON', dt)

	lat, lon, cpc_djf = import_grid(dict_var[var][1], 'SAM-3km', 'CPC', 'DJF', dt)
	lat, lon, cpc_mam = import_grid(dict_var[var][1], 'SAM-3km', 'CPC', 'MAM', dt)
	lat, lon, cpc_jja = import_grid(dict_var[var][1], 'SAM-3km', 'CPC', 'JJA', dt)
	lat, lon, cpc_son = import_grid(dict_var[var][1], 'SAM-3km', 'CPC', 'SON', dt)

	lat, lon, era5_djf = import_grid(dict_var[var][2], 'SAM-3km', 'ERA5', 'DJF', dt)
	lat, lon, era5_mam = import_grid(dict_var[var][2], 'SAM-3km', 'ERA5', 'MAM', dt)
	lat, lon, era5_jja = import_grid(dict_var[var][2], 'SAM-3km', 'ERA5', 'JJA', dt)
	lat, lon, era5_son = import_grid(dict_var[var][2], 'SAM-3km', 'ERA5', 'SON', dt)

	lat, lon, cp_3km_djf = import_grid(var, 'SAM-3km', 'RegCM5', 'DJF', dt)
	lat, lon, cp_3km_mam = import_grid(var, 'SAM-3km', 'RegCM5', 'MAM', dt)
	lat, lon, cp_3km_jja = import_grid(var, 'SAM-3km', 'RegCM5', 'JJA', dt)
	lat, lon, cp_3km_son = import_grid(var, 'SAM-3km', 'RegCM5', 'SON', dt) 

	lat, lon, cp_4km_djf = import_grid(var, 'CSAM-4', 'RegCM5', 'DJF', dt)
	lat, lon, cp_4km_mam = import_grid(var, 'CSAM-4', 'RegCM5', 'MAM', dt)
	lat, lon, cp_4km_jja = import_grid(var, 'CSAM-4', 'RegCM5', 'JJA', dt)
	lat, lon, cp_4km_son = import_grid(var, 'CSAM-4', 'RegCM5', 'SON', dt)	
		
elif var == 'tas':	
	lat, lon, cru_djf = import_grid(dict_var[var][0], 'SAM-3km', 'CRU', 'DJF', dt)
	lat, lon, cru_mam = import_grid(dict_var[var][0], 'SAM-3km', 'CRU', 'MAM', dt)
	lat, lon, cru_jja = import_grid(dict_var[var][0], 'SAM-3km', 'CRU', 'JJA', dt)
	lat, lon, cru_son = import_grid(dict_var[var][0], 'SAM-3km', 'CRU', 'SON', dt)

	lat, lon, era5_djf = import_grid(dict_var[var][1], 'SAM-3km', 'ERA5', 'DJF', dt)
	lat, lon, era5_mam = import_grid(dict_var[var][1], 'SAM-3km', 'ERA5', 'MAM', dt)
	lat, lon, era5_jja = import_grid(dict_var[var][1], 'SAM-3km', 'ERA5', 'JJA', dt)
	lat, lon, era5_son = import_grid(dict_var[var][1], 'SAM-3km', 'ERA5', 'SON', dt)

	lat, lon, cp_3km_djf = import_grid(var, 'SAM-3km', 'RegCM5', 'DJF', dt)
	lat, lon, cp_3km_mam = import_grid(var, 'SAM-3km', 'RegCM5', 'MAM', dt)
	lat, lon, cp_3km_jja = import_grid(var, 'SAM-3km', 'RegCM5', 'JJA', dt)
	lat, lon, cp_3km_son = import_grid(var, 'SAM-3km', 'RegCM5', 'SON', dt) 

	lat, lon, cp_4km_djf = import_grid(var, 'CSAM-4', 'RegCM5', 'DJF', dt)
	lat, lon, cp_4km_mam = import_grid(var, 'CSAM-4', 'RegCM5', 'MAM', dt)
	lat, lon, cp_4km_jja = import_grid(var, 'CSAM-4', 'RegCM5', 'JJA', dt)
	lat, lon, cp_4km_son = import_grid(var, 'CSAM-4', 'RegCM5', 'SON', dt)

else:
	lat, lon, cru_djf = import_grid(dict_var[var][0], 'SAM-3km', 'CRU', 'DJF', dt)
	lat, lon, cru_mam = import_grid(dict_var[var][0], 'SAM-3km', 'CRU', 'MAM', dt)
	lat, lon, cru_jja = import_grid(dict_var[var][0], 'SAM-3km', 'CRU', 'JJA', dt)
	lat, lon, cru_son = import_grid(dict_var[var][0], 'SAM-3km', 'CRU', 'SON', dt)

	lat, lon, cpc_djf = import_grid(dict_var[var][1], 'SAM-3km', 'CPC', 'DJF', dt)
	lat, lon, cpc_mam = import_grid(dict_var[var][1], 'SAM-3km', 'CPC', 'MAM', dt)
	lat, lon, cpc_jja = import_grid(dict_var[var][1], 'SAM-3km', 'CPC', 'JJA', dt)
	lat, lon, cpc_son = import_grid(dict_var[var][1], 'SAM-3km', 'CPC', 'SON', dt)

	lat, lon, cp_3km_djf = import_grid(var, 'SAM-3km', 'RegCM5', 'DJF', dt)
	lat, lon, cp_3km_mam = import_grid(var, 'SAM-3km', 'RegCM5', 'MAM', dt)
	lat, lon, cp_3km_jja = import_grid(var, 'SAM-3km', 'RegCM5', 'JJA', dt)
	lat, lon, cp_3km_son = import_grid(var, 'SAM-3km', 'RegCM5', 'SON', dt) 

	lat, lon, cp_4km_djf = import_grid(var, 'CSAM-4', 'RegCM5', 'DJF', dt)
	lat, lon, cp_4km_mam = import_grid(var, 'CSAM-4', 'RegCM5', 'MAM', dt)
	lat, lon, cp_4km_jja = import_grid(var, 'CSAM-4', 'RegCM5', 'JJA', dt)
	lat, lon, cp_4km_son = import_grid(var, 'CSAM-4', 'RegCM5', 'SON', dt)

if var == 'pr':
	cp_3km_djf = cp_3km_djf
	cp_3km_mam = cp_3km_mam
	cp_3km_jja = cp_3km_jja
	cp_3km_son = cp_3km_son
	cp_4km_djf = cp_4km_djf
	cp_4km_mam = cp_4km_mam
	cp_4km_jja = cp_4km_jja
	cp_4km_son = cp_4km_son
else:
	cp_3km_djf = cp_3km_djf[0]
	cp_3km_mam = cp_3km_mam[0]
	cp_3km_jja = cp_3km_jja[0]
	cp_3km_son = cp_3km_son[0]
	cp_4km_djf = cp_4km_djf
	cp_4km_mam = cp_4km_mam
	cp_4km_jja = cp_4km_jja
	cp_4km_son = cp_4km_son
	
# Plot figure
color = ["#ffffffff", "#d7f0fcff", "#ade0f7ff", "#86c4ebff", "#60a5d6ff", "#4794b3ff", "#49a67cff", "#55b848ff", 
"#9ecf51ff", "#ebe359ff", "#f7be4aff", "#f58433ff", "#ed5a28ff", "#de3728ff", "#cc1f27ff", "#b01a1fff", "#911419ff"]

dict_plot = {
'pr': ['Precipitation (mm d$^-$$^1$)', np.arange(0, 18, 1), matplotlib.colors.ListedColormap(color)],
'tas': ['Air temperature (°C)', np.arange(0, 34, 2), cm.jet],
'tasmax': ['Maximum air temperature (°C)', np.arange(6, 40, 2), cm.jet],
'tasmin': ['Minimum air temperature (°C)', np.arange(-4, 30, 2), cm.jet]
}

font_size = 8
	
if var == 'pr':
	fig = plt.figure(figsize=(10, 6))

	ax = fig.add_subplot(4, 5, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cru_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(a) CRU DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cpc_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(b) CPC DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(c) ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_3km_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(d) CP-3km DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_4km_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(e) CP-4km DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cru_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(f) CRU MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cpc_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(g) CPC MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(h) ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 9)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_3km_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(i) CP-3km MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 10)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_4km_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(j) CP-4km MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 11)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cru_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(k) CRU JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 12)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cpc_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(l) CPC JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 13)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(m) ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 14)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_3km_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(n) CP-3km JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 15)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_4km_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(o) CP-4km JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 16)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cru_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(p) CRU SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 17)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cpc_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(q) CPC SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 18)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(r) ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 19)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_3km_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(s) CP-3km SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 20)  
	map, xx, yy = basemap(lat, lon)
	#cp_4km_son = maskoceans(xx, yy, cp_4km_son)
	plt_map = map.contourf(xx, yy, cp_4km_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(t) CP-4km SON', loc='left', fontsize=font_size, fontweight='bold')

elif var == 'tas':
	fig = plt.figure(figsize=(8, 6))

	ax = fig.add_subplot(4, 4, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cru_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(a) CRU DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(b) ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_3km_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(c) CP-3km DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_4km_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(d) CP-4km DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cru_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(e) CRU MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(f) ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_3km_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(g) CP-3km MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_4km_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(j) CP-4km MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 9)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cru_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(i) CRU JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 10)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(j) ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 11)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_3km_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(k) CP-3km JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 12)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_4km_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(l) CP-4km JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 13)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cru_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(m) CRU SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 14)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(n) ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 15)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_3km_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(o) CP-3km SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 16)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_4km_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(p) CP-4km SON', loc='left', fontsize=font_size, fontweight='bold')

else:
	fig = plt.figure(figsize=(8, 6))
	
	ax = fig.add_subplot(4, 4, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cru_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(a) CRU DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cpc_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(b) CPC DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_3km_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(c) CP-3km DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_4km_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(d) CP-4km DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cru_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(e) CRU MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cpc_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(f) CPC MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_3km_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(g) CP-3km MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_4km_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(j) CP-4km MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 9)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cru_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(i) CRU JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 10)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cpc_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(j) CPC JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 11)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_3km_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(k) CP-3km JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 12)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_4km_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(l) CP-4km JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 13)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cru_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(m) CRU SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 14)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cpc_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(n) CPC SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 15)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_3km_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(o) CP-3km SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 16)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, cp_4km_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
	plt.title(u'(p) CP-4km SON', loc='left', fontsize=font_size, fontweight='bold')

# Set colobar
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.3, 0.015, 0.4]))
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cp_3km-4km'.format(path)
name_out = 'pyplt_maps_clim_{0}_CP-RegCM5_SAM-3km_{1}.png'.format(var, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
