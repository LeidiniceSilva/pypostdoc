# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot pdf"

import os
import netCDF4
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.pyplot as plt

var = 'pr'
domain = 'CSAM-3'
dt = '200101'

path = '/marconi/home/userexternal/mdasilva'

			
def import_obs(param, domain, dataset):

	#arq   = '{0}/user/mdasilva/CORDEX/post_evaluate/obs/{1}_{2}_{3}_day_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)
	arq   = '{0}/user/mdasilva/CORDEX/post_evaluate/obs/{1}_{2}_{3}_1hr_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)		
	data  = netCDF4.Dataset(arq)
	var   = data.variables['tp'][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean


def import_rcm(param, domain, dataset):

	#arq   = '{0}/user/mdasilva/CORDEX/post_evaluate/rcm/{1}_{2}_{3}_day_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)
	arq   = '{0}/user/mdasilva/CORDEX/post_evaluate/rcm/{1}_{2}_{3}_1hr_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)
	data  = netCDF4.Dataset(arq)
	var   = data.variables['pr'][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:-1,:,:]
	
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

# Convert array in list
era5_list = era5.flatten()
regcm5_list = regcm5.flatten()

# Round values
round_era5 = np.round(era5_list,0)
round_regcm5 = np.round(regcm5_list,0)

# Filter 0 mm/day
filter_era5 = round_era5[round_era5 > 0.]
filter_regcm5 = round_regcm5[round_regcm5 > 0.]

# Compute frequency
x_pdf_era5, pdf_era5 = np.unique(filter_era5, return_counts=True) 
x_pdf_regcm5, pdf_regcm5 = np.unique(filter_regcm5, return_counts=True) 

# Plot figure  
fig = plt.figure()
font_size = 8

ax = fig.add_subplot(1, 1, 1)  
plt.plot(x_pdf_era5, pdf_era5, marker='o', markersize=4, mfc='red', mec='red', alpha=0.65, linestyle='None', label='ERA5')
plt.plot(x_pdf_regcm5, pdf_regcm5, marker='o', markersize=4, mfc='black', mec='black', alpha=0.65, linestyle='None', label='RegCM5')
plt.title('(a)', loc='left', fontsize=font_size, fontweight='bold') 
plt.yscale('log')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
#plt.xlabel('Precipitation (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.xlabel('Precipitation (mm h$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.legend(loc=1, ncol=2, fontsize=font_size, shadow=True)

# Path out to save figure
path_out = '{0}/user/mdasilva/CORDEX/figs'.format(path)
#name_out = 'pyplt_pdf_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
name_out = 'pyplt_pdf_{0}_{1}_RegCM5_1hr_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
