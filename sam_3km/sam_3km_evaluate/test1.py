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
domain = 'SAM-3km'
idt, fdt = '2018', '2018'
dt = '{0}-{1}'.format(idt, fdt)

path = '/marconi/home/userexternal/mdasilva'
	
	
def import_obs(param, domain, dataset, season):

	arq   = '{0}/user/mdasilva/SAM-3km_v4/post/obs/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean


def import_rcm(param, domain, dataset, season):

	arq   = '{0}/user/mdasilva/SAM-3km_v4/post/rcm/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
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
dict_var = {
'pr': ['pre', 'pre', 'precip', 'sat_gauge_precip', 'tp'],
'tas': ['tmp', 'tmp', 't2m'],
'tasmax': ['tmx', 'tmax', 'mx2t'],
'tasmin': ['tmn', 'tmin', 'mn2t'],
'clt': ['cld', 'tcc'],
'rsnl': ['msnlwrf'],
}

if var == 'pr':
	database = 'CRU'

	lat, lon, obs_mon = import_obs(dict_var[var][1], domain, database, 'mon')
	lat, lon, regcm_mon = import_rcm(var, domain, 'RegCM5', 'mon')
	mbe_mon = compute_mbe(regcm_mon, obs_mon)
	
elif var == 'tas':
	database = 'CRU'

	lat, lon, obs_mon = import_obs(dict_var[var][1], domain, database, 'mon')
	lat, lon, regcm_mon = import_rcm(var, domain, 'RegCM5', 'mon')
	mbe_mon = compute_mbe(regcm_mon, obs_mon)
	
elif var == 'tasmax':
	database = 'CRU'

	lat, lon, obs_mon = import_obs(dict_var[var][0], domain, database, 'mon')
	lat, lon, regcm_mon = import_rcm(var, domain, 'RegCM5', 'mon')
	mbe_mon = compute_mbe(regcm_mon, obs_mon)

elif var == 'tasmin':
	database = 'CRU'

	lat, lon, obs_mon = import_obs(dict_var[var][0], domain, database, 'mon')
	lat, lon, regcm_mon = import_rcm(var, domain, 'RegCM5', 'mon')
	mbe_mon = compute_mbe(regcm_mon, obs_mon)	
	
elif var == 'clt':
	database = 'CRU'

	lat, lon, obs_mon = import_obs(dict_var[var][0], domain, database, 'mon')
	lat, lon, regcm_mon = import_rcm(var, domain, 'RegCM5', 'mon')
	mbe_mon = compute_mbe(regcm_mon, obs_mon)
	
else:
	database = 'ERA5'
	
	lat, lon, obs_mon = import_obs(dict_var[var][0], domain, database, 'mon')
	lat, lon, regcm_mon = import_rcm(var, domain, 'RegCM5', 'mon')
	mbe_mon = compute_mbe(regcm_mon, obs_mon)

# Plot figure  
fig = plt.figure(figsize=(10, 7))
 
font_size = 8

dict_plot = {
'pr': ['Bias of  precipitation (mm d$^-$$^1$)', np.arange(-10, 11, 1), cm.BrBG],
'tas': ['Bias of air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'tasmax': ['Bias of maximum air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'tasmin': ['Bias of minimum air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'clt': ['Bias of total cloud cover (%)', np.arange(-70, 80, 10), cm.RdGy],
'rsnl': ['Bias of surface net upward longwave flux (W mm$^-$$^2$)', np.arange(-60, 65, 5), cm.coolwarm]
}
	
ax = fig.add_subplot(3, 3, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mon[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(a) RegCM5-{0} Jan'.format(database), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mon[1], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(b) RegCM5-{0} Feb'.format(database), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mon[2], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(c) RegCM5-{0} Mar'.format(database), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 4)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mon[3], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(d) RegCM5-{0} Apr'.format(database), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 5)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mon[4], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(e) RegCM5-{0} May'.format(database), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 6)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mon[5], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(f) RegCM5-{0} Jun'.format(database), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 7)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mon[6], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(g) RegCM5-{0} Jul'.format(database), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 8)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mon[7], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(h) RegCM5-{0} Aug'.format(database), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 9)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mon[8], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(i) RegCM5-{0} Sep'.format(database), loc='left', fontsize=font_size, fontweight='bold')

cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.3, 0.018, 0.4]))
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km_v5/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
