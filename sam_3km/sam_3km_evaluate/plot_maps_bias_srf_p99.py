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
domain = 'SAM-3km'
path = '/marconi/home/userexternal/mdasilva'

if freq == 'hourly':
	dt = '1hr_2018-2021'
else:
	dt = '2018-2021'


def import_obs(param, domain, dataset):

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/obs/p99_{1}_{2}_{3}_lonlat.nc'.format(path, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	
		
def import_cp_3km(param, domain, dataset):

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/p99_{1}_{2}_{3}_lonlat.nc'.format(path, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean


def basemap(lat, lon):
	
	map = Basemap(projection='cyl', llcrnrlon=-80., llcrnrlat=-38., urcrnrlon=-34.,urcrnrlat=-8., resolution='c')
	map.drawmeridians(np.arange(-80., -34., 12.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-38., -8., 6.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black')
	
	lons, lats = np.meshgrid(lon, lat)
	xx, yy = map(lons,lats)
	
	return map, xx, yy
	
	
# Import model and obs dataset
if freq == 'hourly':
	lat, lon, era5 = import_obs('tp', domain, 'ERA5')
	lat, lon, regcm = import_cp_3km('pr', domain, 'RegCM5')
	
	mbe_regcm_era5 = regcm - era5
else:
	lat, lon, cpc = import_obs('precip', domain, 'CPC')
	lat, lon, era5 = import_obs('tp', domain, 'ERA5')
	lat, lon, regcm = import_cp_3km('pr', domain, 'RegCM5')
	
	mbe_regcm_cpc = regcm - cpc
	mbe_regcm_era5 = regcm - era5

# Plot figure
fig = plt.figure()   
font_size = 8
	
if freq == 'hourly':
	levs = np.arange(-5, 5.5, 0.5)
	legend = 'Hourly p99 (mm h$^-$$^1$)'

	ax = fig.add_subplot(1, 1, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_era5[0], levels=levs, cmap=cm.BrBG, extend='both') 
	plt.title(u'(a) RegCM5 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	cbar = plt.colorbar(plt_map)
	cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)

else:
	levs = np.arange(-60, 70, 10)
	legend = 'Daily p99 (mm d$^-$$^1$)'

	ax = fig.add_subplot(2, 1, 1)
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_cpc[0], levels=levs, cmap=cm.BrBG, extend='both') 
	plt.title(u'(a) RegCM5 - CPC', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(2, 1, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_era5[0], levels=levs, cmap=cm.BrBG, extend='both') 
	plt.title(u'(b) RegCM5 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.75, 0.3, 0.03, 0.4]))
	cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_maps_bias_{0}_RegCM5_{1}_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
