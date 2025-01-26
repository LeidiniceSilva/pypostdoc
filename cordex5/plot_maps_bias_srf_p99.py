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
freq = 'hourly'
domain = 'CSAM-3'
idt, fdt = '2000', '2009'

path = '/marconi/home/userexternal/mdasilva'


def import_obs(param, domain, dataset):

	if freq == 'hourly':
		dt = '1hr_{0}-{1}'.format(idt, fdt)
	else:
		dt = '{0}-{1}'.format(idt, fdt)
		
	if param == 'pr':
		param_ = 'tp'
	else:
		param_ = param
		
	arq = '{0}/user/mdasilva/CORDEX/post_evaluate/obs/p99_{1}_{2}_{3}.nc'.format(path, domain, dataset, dt)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param_][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	
		
def import_rcm(param, domain, dataset):

	if freq == 'hourly':
		dt = '1hr_{0}-{1}'.format('2000', '2005')
	else:
		dt = '{0}-{1}'.format('2000', '2005')
		
	arq = '{0}/user/mdasilva/CORDEX/post_evaluate/rcm/p99_{1}_{2}_{3}_lonlat.nc'.format(path, domain, dataset, dt)		
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
	lat, lon, cmorph = import_obs('cmorph', domain, 'CMORPH')
	lat, lon, era5   = import_obs('pr', domain, 'ERA5')
	lat, lon, regcm  = import_rcm('pr', domain, 'RegCM5')

	mbe_regcm_cmorph = regcm - cmorph
	mbe_regcm_era5   = regcm - era5
else:
	lat, lon, cpc    = import_obs('precip', domain, 'CPC')
	lat, lon, cmorph = import_obs('cmorph', domain, 'CMORPH')
	lat, lon, mswep  = import_obs('precipitation', domain, 'MSWEP')
	lat, lon, era5   = import_obs('pr', domain, 'ERA5')
	lat, lon, regcm  = import_rcm('pr', domain, 'RegCM5')
	
	mbe_regcm_cpc    = regcm - cpc
	mbe_regcm_cmorph = regcm - cmorph
	mbe_regcm_mswep  = regcm - mswep
	mbe_regcm_era5   = regcm - era5

# Plot figure
font_size = 8
		
if freq == 'hourly':
	fig = plt.figure(figsize=(10, 4))   
	levs = np.arange(-5, 5.5, 0.5)
	legend = 'Hourly p99 (mm h$^-$$^1$)'
	dt = '1hr_{0}-{1}'.format(idt, fdt)
	
	ax = fig.add_subplot(1, 2, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_cmorph[0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(a) CPM3 - CMORPH', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(1, 2, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_era5[0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(b) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
	
	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.91, 0.2, 0.02, 0.6]))
	cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)

else:
	fig = plt.figure(figsize=(10, 4))   
	levs = np.arange(-60, 70, 10)
	legend = 'Daily p99 (mm d$^-$$^1$)'
	dt = '{0}-{1}'.format(idt, fdt)

	ax = fig.add_subplot(1, 4, 1)
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_cpc[0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(a) CPM3 - CPC', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(1, 4, 2)
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_cmorph[0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(b) CPM3 - CMORPH', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(1, 4, 3)
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_mswep[0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(c) CPM3 - MSWEP', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(1, 4, 4)
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_era5[0], levels=levs, cmap=cm.BrBG, extend='neither') 
	plt.title(u'(d) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.91, 0.35, 0.015, 0.3]))
	cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/user/mdasilva/CORDEX/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
