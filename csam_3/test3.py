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
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii
from mpl_toolkits.basemap import Basemap
from import_climate_tools import compute_mbe

var = 'pr_freq'
domain = 'CSAM-3'
dt = '200101'

path = '/marconi/home/userexternal/mdasilva'

			
def import_obs(param, domain, dataset):

	#arq   = '{0}/user/mdasilva/CORDEX/post_evaluate/obs/{1}_{2}_{3}_day_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)
	arq   = '{0}/user/mdasilva/CORDEX/post_evaluate/obs/{1}_{2}_{3}_1hr_{4}_th0.5_lonlat.nc'.format(path, param, domain, dataset, dt)		
	data  = netCDF4.Dataset(arq)
	var   = data.variables['tp'][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]
	
	return lat, lon, mean


def import_rcm(param, domain, dataset):

	#arq   = '{0}/user/mdasilva/CORDEX/post_evaluate/rcm/{1}_{2}_{3}_day_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)
	arq   = '{0}/user/mdasilva/CORDEX/post_evaluate/rcm/{1}_{2}_{3}_1hr_{4}_th0.5_lonlat.nc'.format(path, param, domain, dataset, dt)
	data  = netCDF4.Dataset(arq)
	var   = data.variables['pr'][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]
	
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
lat, lon, era5 = import_obs(var, domain, 'ERA5')
lat, lon, regcm5 = import_rcm(var, domain, 'RegCM5')

mbe_regcm5_era5 = compute_mbe(regcm5, era5)

print(era5.shape)
print(regcm5.shape)
print(mbe_regcm5_era5.shape)

# Plot figure   
#dict_plot = {'pr_freq': ['Frequency (%) Daily', np.arange(-60, 65, 5), cm.BrBG]}
dict_plot = {'pr_freq': ['Frequency (%) Hourly', np.arange(-18, 20, 2), cm.BrBG]}

font_size = 8
	
fig = plt.figure()

ax = fig.add_subplot(1, 1, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_regcm5_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
plt.title(u'(a) RegCM5-ERA5 Jan', loc='left', fontsize=font_size, fontweight='bold')

cbar = plt.colorbar(plt_map)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/user/mdasilva/CORDEX/figs'.format(path)
#name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_1hr_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
