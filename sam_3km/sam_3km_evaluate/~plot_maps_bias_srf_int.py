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

var = 'pr_int'
domain = 'SAM-3km'
freq = 'hourly'
path = '/marconi/home/userexternal/mdasilva'


def import_obs(param, domain, dataset, season):

	if freq == 'hourly':
		dt = '1hr_{0}_2018-2021_th0.5'.format(season)
	else:
		dt = '{0}_2018-2021'.format(season)

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/obs/{1}_i_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	
		
def import_cp_3km(param, domain, dataset, season):

	if freq == 'hourly':
		dt = '1hr_{0}_2018-2021_th0.5'.format(season)
	else:
		dt = '{0}_2018-2021'.format(season)

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/{1}_freq_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
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
	lat, lon, era5_djf = import_obs('tp', domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_obs('tp', domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_obs('tp', domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_obs('tp', domain, 'ERA5', 'SON')

	lat, lon, regcm_djf = import_cp_3km('pr', domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_cp_3km('pr', domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_cp_3km('pr', domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_cp_3km('pr', domain, 'RegCM5', 'SON')

	mbe_regcm_era5_djf = regcm_djf - era5_djf
	mbe_regcm_era5_mam = regcm_mam - era5_mam
	mbe_regcm_era5_jja = regcm_jja - era5_jja
	mbe_regcm_era5_son = regcm_son - era5_son
else:
	lat, lon, cpc_djf = import_obs('precip', domain, 'CPC', 'DJF')
	lat, lon, cpc_mam = import_obs('precip', domain, 'CPC', 'MAM')
	lat, lon, cpc_jja = import_obs('precip', domain, 'CPC', 'JJA')
	lat, lon, cpc_son = import_obs('precip', domain, 'CPC', 'SON')

	lat, lon, regcm_djf = import_cp_3km('pr', domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_cp_3km('pr', domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_cp_3km('pr', domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_cp_3km('pr', domain, 'RegCM5', 'SON')

	mbe_regcm_cpc_djf = regcm_djf - cpc_djf
	mbe_regcm_cpc_mam = regcm_mam - cpc_mam
	mbe_regcm_cpc_jja = regcm_jja - cpc_jja
	mbe_regcm_cpc_son = regcm_son - cpc_son

# Plot figure
fig = plt.figure(figsize=(4, 8))   
font_size = 8

if freq == 'hourly':
	levs = np.arange(-6, 7, 0.5)
	legend = 'Frequency (%)'
	dt = '1hr_{0}_2018-2021_th0.5'.format(season)

	ax = fig.add_subplot(4, 1, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_era5_djf[0][0], levels=np.arange(-45, 50, 5), cmap=cm.BrBG, extend='both') 
	plt.title(u'(a) RegCM5 - ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 1, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_era5_mam[0][0], levels=np.arange(-45, 50, 5), cmap=cm.BrBG, extend='both') 
	plt.title(u'(b) RegCM5 - ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 1, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_era5_jja[0][0], levels=np.arange(-45, 50, 5), cmap=cm.BrBG, extend='both') 
	plt.title(u'(c) RegCM5 - ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 1, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_era5_son[0][0], levels=np.arange(-45, 50, 5), cmap=cm.BrBG, extend='both') 
	plt.title(u'(d) RegCM5 - ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')

	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.84, 0.3, 0.03, 0.4]))
	cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)
	
else:
	levs = np.arange(-24, 28, 4)
	legend = 'Frequency (%)'
	dt = '{0}_2018-2021'.format(season)
			
	ax = fig.add_subplot(4, 1, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_cpc_djf[0][0], levels=np.arange(-45, 50, 5), cmap=cm.BrBG, extend='both') 
	plt.title(u'(a) RegCM5 - CPC DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 1, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_cpc_mam[0][0], levels=np.arange(-45, 50, 5), cmap=cm.BrBG, extend='both') 
	plt.title(u'(b) RegCM5 - CPC MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 1, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_cpc_jja[0][0], levels=np.arange(-45, 50, 5), cmap=cm.BrBG, extend='both') 
	plt.title(u'(c) RegCM5 - CPC JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 1, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_regcm_cpc_son[0][0], levels=np.arange(-45, 50, 5), cmap=cm.BrBG, extend='both') 
	plt.title(u'(d) RegCM5 - CPC SON', loc='left', fontsize=font_size, fontweight='bold')

	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.84, 0.3, 0.03, 0.4]))
	cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_maps_bias_{0}_freq_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()




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
domain = 'SAM-3km'
dt = '2018-2021'
path = '/marconi/home/userexternal/mdasilva'


def import_obs(param, domain, dataset, season):

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/obs/{1}_int_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	
		
def import_cp_3km(param, domain, dataset, season):

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/{1}_int_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
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
lat, lon, cpc_djf = import_obs('precip', domain, 'CPC', 'DJF')
lat, lon, cpc_mam = import_obs('precip', domain, 'CPC', 'MAM')
lat, lon, cpc_jja = import_obs('precip', domain, 'CPC', 'JJA')
lat, lon, cpc_son = import_obs('precip', domain, 'CPC', 'SON')

lat, lon, regcm_djf = import_cp_3km('pr', domain, 'RegCM5', 'DJF')
lat, lon, regcm_mam = import_cp_3km('pr', domain, 'RegCM5', 'MAM')
lat, lon, regcm_jja = import_cp_3km('pr', domain, 'RegCM5', 'JJA')
lat, lon, regcm_son = import_cp_3km('pr', domain, 'RegCM5', 'SON')

mbe_regcm_cpc_djf = regcm_djf - cpc_djf
mbe_regcm_cpc_mam = regcm_mam - cpc_mam
mbe_regcm_cpc_jja = regcm_jja - cpc_jja
mbe_regcm_cpc_son = regcm_son - cpc_son

# Plot figure
fig = plt.figure(figsize=(4, 8))   
font_size = 8
	
ax = fig.add_subplot(4, 1, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_regcm_cpc_djf[0][0], levels=np.arange(-24, 28, 4), cmap=cm.BrBG, extend='both') 
plt.title(u'(a) RegCM5 - CPC DJF', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 1, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_regcm_cpc_mam[0][0], levels=np.arange(-24, 28, 4), cmap=cm.BrBG, extend='both') 
plt.title(u'(b) RegCM5 - CPC MAM', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 1, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_regcm_cpc_jja[0][0], levels=np.arange(-24, 28, 4), cmap=cm.BrBG, extend='both') 
plt.title(u'(c) RegCM5 - CPC JJA', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 1, 4)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_regcm_cpc_son[0][0], levels=np.arange(-24, 28, 4), cmap=cm.BrBG, extend='both') 
plt.title(u'(d) RegCM5 - CPC SON', loc='left', fontsize=font_size, fontweight='bold')

cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.84, 0.3, 0.03, 0.4]))
cbar.set_label('Intensity (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_maps_bias_{0}_int_RegCM5_{1}_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
