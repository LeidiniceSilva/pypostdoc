# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot bias maps"

import os
import netCDF4
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from mpl_toolkits.basemap import Basemap
from import_climate_tools import compute_mbe

var = 'pr'
freq = 'daily'
domain = 'SAM-3km'
idt, fdt = '2018', '2018'
dt = '{0}-{1}'.format(idt, fdt)

path = '/marconi/home/userexternal/mdasilva'


def import_obs(param, domain, dataset, season):

	arq   = '{0}/user/mdasilva/SAM-3km_v5/post/obs/{1}_int_{2}_{3}_{4}_lonlat.nc'.format(path, param, domain, dataset, season)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	
		
def import_rcm(param, domain, dataset, season):

	arq   = '{0}/user/mdasilva/SAM-3km_v5/post/rcm/{1}_int_{2}_{3}_{4}_lonlat.nc'.format(path, param, domain, dataset, season)	
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
lat, lon, obs_jan = import_obs('precip', domain, 'CPC', '20180101')
lat, lon, obs_feb = import_obs('precip', domain, 'CPC', '20180201')
lat, lon, obs_mar = import_obs('precip', domain, 'CPC', '20180301')
lat, lon, obs_apr = import_obs('precip', domain, 'CPC', '20180401')
lat, lon, obs_may = import_obs('precip', domain, 'CPC', '20180501')
lat, lon, obs_jun = import_obs('precip', domain, 'CPC', '20180601')
lat, lon, obs_jul = import_obs('precip', domain, 'CPC', '20180701')
lat, lon, obs_aug = import_obs('precip', domain, 'CPC', '20180801')
lat, lon, obs_sep = import_obs('precip', domain, 'CPC', '20180801')

lat, lon, regcm_jan = import_rcm('pr', domain, 'RegCM5', '20180101')
lat, lon, regcm_feb = import_rcm('pr', domain, 'RegCM5', '20180201')
lat, lon, regcm_mar = import_rcm('pr', domain, 'RegCM5', '20180301')
lat, lon, regcm_apr = import_rcm('pr', domain, 'RegCM5', '20180401')
lat, lon, regcm_may = import_rcm('pr', domain, 'RegCM5', '20180501')
lat, lon, regcm_jun = import_rcm('pr', domain, 'RegCM5', '20180601')
lat, lon, regcm_jul = import_rcm('pr', domain, 'RegCM5', '20180701')
lat, lon, regcm_aug = import_rcm('pr', domain, 'RegCM5', '20180801')
lat, lon, regcm_sep = import_rcm('pr', domain, 'RegCM5', '20180901')
	
mbe_jan = compute_mbe(regcm_jan, obs_jan)
mbe_feb = compute_mbe(regcm_feb, obs_feb)
mbe_mar = compute_mbe(regcm_mar, obs_mar)
mbe_apr = compute_mbe(regcm_apr, obs_apr)
mbe_may = compute_mbe(regcm_may, obs_may)
mbe_jun = compute_mbe(regcm_jun, obs_jun)
mbe_jul = compute_mbe(regcm_jul, obs_jul)
mbe_aug = compute_mbe(regcm_aug, obs_aug)
mbe_sep = compute_mbe(regcm_sep, obs_sep)

# Plot figure
fig = plt.figure(figsize=(10, 7))   
font_size = 8
	
levs = np.arange(-25, 28, 3)
legend = 'Intensity (mm d$^-$$^1$)'

ax = fig.add_subplot(3, 3, 1)
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_jan[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(a) RegCM5-CPC Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 2)
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_feb[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(b) RegCM5-CPC Feb', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 3)
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mar[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(c) RegCM5-CPC Mar', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 4)
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_apr[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(d) RegCM5-CPC Apr', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 5)
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_may[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(e) RegCM5-CPC May', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 6)
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_jun[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(f) RegCM5-CPC Jun', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 7)
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_jul[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(g) RegCM5-CPC Jul', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 8)
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_aug[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(h) RegCM5-CPC Aug', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 9)
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_sep[0][0], levels=levs, cmap=cm.BrBG, extend='neither') 
plt.title(u'(i) RegCM5-CPC Sep', loc='left', fontsize=font_size, fontweight='bold')

cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.3, 0.018, 0.4]))
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km_v5/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_int_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
