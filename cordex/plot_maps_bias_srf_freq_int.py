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

var = 'pr'
stats = 'int'
freq = 'daily'
domain = 'CSAM-3'
idt, fdt = '2000', '2009'

path = '/marconi/home/userexternal/mdasilva'


def import_obs(param, domain, dataset, season):

	if freq == 'hourly':
		dt = '1hr_{0}_{1}-{2}_th0.5'.format(season, idt, fdt)
	else:
		dt = '{0}_{1}-{2}'.format(season, idt, fdt)

	arq   = '{0}/user/mdasilva/CORDEX/post_evaluate/obs/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, stats, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	
		
def import_rcm(param, domain, dataset, season):

	if freq == 'hourly':
		dt = '1hr_{0}_{1}-{2}_th0.5'.format(season, 2000, 2005)
	else:
		dt = '{0}_{1}-{2}'.format(season, 2000, 2005)

	arq   = '{0}/user/mdasilva/CORDEX/post_evaluate/rcm/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, stats, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	

def basemap(lat, lon):
	
	map = Basemap(projection='cyl', llcrnrlon=-80., llcrnrlat=-38., urcrnrlon=-34.,urcrnrlat=-10., resolution='c')
	map.drawmeridians(np.arange(-80., -34., 12.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-38., -8., 6.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black')
	
	lons, lats = np.meshgrid(lon, lat)
	xx, yy = map(lons,lats)
	
	return map, xx, yy
	
	
# Import model and obs dataset
if freq == 'hourly':
	lat, lon, era5_djf = import_obs('pr', domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_obs('pr', domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_obs('pr', domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_obs('pr', domain, 'ERA5', 'SON')

	lat, lon, regcm_djf = import_rcm('pr', domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_rcm('pr', domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_rcm('pr', domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_rcm('pr', domain, 'RegCM5', 'SON')

	mbe_regcm_era5_djf = regcm_djf - era5_djf
	mbe_regcm_era5_mam = regcm_mam - era5_mam
	mbe_regcm_era5_jja = regcm_jja - era5_jja
	mbe_regcm_era5_son = regcm_son - era5_son
else:
	lat, lon, cpc_djf = import_obs('precip', domain, 'CPC', 'DJF')
	lat, lon, cpc_mam = import_obs('precip', domain, 'CPC', 'MAM')
	lat, lon, cpc_jja = import_obs('precip', domain, 'CPC', 'JJA')
	lat, lon, cpc_son = import_obs('precip', domain, 'CPC', 'SON')

	lat, lon, mswep_djf = import_obs('precipitation', domain, 'MSWEP', 'DJF')
	lat, lon, mswep_mam = import_obs('precipitation', domain, 'MSWEP', 'MAM')
	lat, lon, mswep_jja = import_obs('precipitation', domain, 'MSWEP', 'JJA')
	lat, lon, mswep_son = import_obs('precipitation', domain, 'MSWEP', 'SON')

	lat, lon, trmm_djf = import_obs('hrf', domain, 'TRMM', 'DJF')
	lat, lon, trmm_mam = import_obs('hrf', domain, 'TRMM', 'MAM')
	lat, lon, trmm_jja = import_obs('hrf', domain, 'TRMM', 'JJA')
	lat, lon, trmm_son = import_obs('hrf', domain, 'TRMM', 'SON')
	
	lat, lon, era5_djf = import_obs('pr', domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_obs('pr', domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_obs('pr', domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_obs('pr', domain, 'ERA5', 'SON')
	
	lat, lon, regcm_djf = import_rcm('pr', domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_rcm('pr', domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_rcm('pr', domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_rcm('pr', domain, 'RegCM5', 'SON')

	mbe_regcm_cpc_djf = regcm_djf - cpc_djf
	mbe_regcm_cpc_mam = regcm_mam - cpc_mam
	mbe_regcm_cpc_jja = regcm_jja - cpc_jja
	mbe_regcm_cpc_son = regcm_son - cpc_son
	
	mbe_regcm_mswep_djf = regcm_djf - mswep_djf
	mbe_regcm_mswep_mam = regcm_mam - mswep_mam
	mbe_regcm_mswep_jja = regcm_jja - mswep_jja
	mbe_regcm_mswep_son = regcm_son - mswep_son

	mbe_regcm_trmm_djf = regcm_djf - trmm_djf
	mbe_regcm_trmm_mam = regcm_mam - trmm_mam
	mbe_regcm_trmm_jja = regcm_jja - trmm_jja
	mbe_regcm_trmm_son = regcm_son - trmm_son
	
	mbe_regcm_era5_djf = regcm_djf - era5_djf
	mbe_regcm_era5_mam = regcm_mam - era5_mam
	mbe_regcm_era5_jja = regcm_jja - era5_jja
	mbe_regcm_era5_son = regcm_son - era5_son
	
# Plot figure
fig = plt.figure(figsize=(8, 6))   
font_size = 8

if stats == 'freq':
	if freq == 'hourly':
		levs = np.arange(-22, 23, 1)
		legend = 'Frequency (%)'
		dt = '1hr_{0}-{1}_th0.5'.format(idt, fdt)
	else:
		levs = np.arange(-44, 46, 2)
		legend = 'Frequency (%)'
		dt = '{0}-{1}'.format(idt, fdt)
else:
	if freq == 'hourly':
		levs = np.arange(-6, 6.5, 0.5)
		legend = 'Intensity (mm h$^-$$^1$)'
		dt = '1hr_{0}-{1}_th0.5'.format(idt, fdt)
	else:
		levs = np.arange(-22, 23, 1)
		legend = 'Intensity (mm d$^-$$^1$)'
		dt = '{0}-{1}'.format(idt, fdt)

if freq == 'hourly':
	ax = fig.add_subplot(4, 4, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_era5_djf[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(a) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'DJF', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_era5_mam[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(b) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'MAM', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_era5_jja[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(c) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'JJA', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_era5_son[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(d) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'SON', labelpad=20, fontsize=font_size, fontweight='bold')

	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.84, 0.3, 0.03, 0.4]))
	cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)
	
else:	
	ax = fig.add_subplot(4, 4, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_cpc_djf[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(a) CPM3 - CPC', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'DJF', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_mswep_djf[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(b) CPM3 - MSWEP', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_trmm_djf[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(c) CPM3 - TRMM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_era5_djf[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(d) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_cpc_mam[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(e) CPM3 - CPC', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'MAM', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_mswep_mam[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(f) CPM3 - MSWEP', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_trmm_mam[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(g) CPM3 - TRMM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_era5_mam[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(h) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 9)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_cpc_jja[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(i) CPM3 - CPC', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'JJA', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 10)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_mswep_jja[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(j) CPM3 - MSWEP', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 11)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_trmm_jja[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(k) CPM3 - TRMM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 12)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_era5_jja[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(l) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 13)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_cpc_son[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(m) CPM3 - CPC', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'SON', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 14)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_mswep_son[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(n) CPM3 - MSWEP', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 15)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_trmm_son[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(o) CPM3 - TRMM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 16)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_era5_son[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(p) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.94, 0.3, 0.03, 0.4]))
	cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)
# Path out to save figure
path_out = '{0}/user/mdasilva/CORDEX/figs/v1'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_{2}_RegCM5_{3}.png'.format(var, stats, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
