# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot map of genesis density"

import os
import netCDF4
import datetime
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeat

from scipy import signal, misc
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/SAM-3km'


def read_dat_file(filename):

	data = []
	with open(filename, 'r') as file:
		lines = file.readlines()
		header = []
		rows = []
	
	# Iterate over lines in the file
	for line in lines:
		line = line.strip().split()
		
		# If the line contains 6 elements, it's considered a header
		if len(line) == 6:
			if rows:  # If we have rows, append them to data
				data.append((header, rows))
				rows = []  # Reset rows
			header = line
		else:
			rows.append(line)
	
	# Append the last header and rows to data
	if header and rows:
		data.append((header, rows))
	
	return data


def open_dat_file(dataset, yr_init, yr_end):

	lat, lon = [], []
	for yr in range(yr_init, yr_end+1):
	
		data = read_dat_file('{0}/postproc/cyclone/{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
		for i, (header, rows) in enumerate(data):
			lat.append(rows[0][1])
			lon.append(rows[0][2])

	return lat, lon
	

def import_data(param, dataset):

	arq   = '{0}/postproc/cyclone/{2}/{1}_SAM-3km_{2}_day_2018-2021_lonlat.nc'.format(path, param, dataset)		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]	
	lon   = data.variables['lon'][:]
	mean = np.nanmean(var[:][:,:,:], axis=0)
		
	return lat, lon, mean
	
	
# Import model and obs dataset
lat_era5, lon_era5 = open_dat_file('ERA5', 2018, 2021)
lat_regcm5, lon_regcm5 = open_dat_file('RegCM5', 2018, 2021)
lat_wrf415, lon_wrf415 = open_dat_file('WRF415', 2018, 2021)

# Import model and obs dataset 
lat_, lon_, data_mean = import_data('tp', 'ERA5')

lat_era5_ls = [float(i) for i in lat_era5]
lon_era5_ls = [float(i) for i in lon_era5]

lat_regcm5_ls = [float(i) for i in lat_regcm5]
lon_regcm5_ls = [float(i) for i in lon_regcm5]

lat_wrf415_ls = [float(i) for i in lat_wrf415]
lon_wrf415_ls = [float(i) for i in lon_wrf415]

# Plot figure
fig, axes = plt.subplots(2,2, figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2), (ax3, ax4) = axes
fig.delaxes(ax4)

font_size = 10
states_provinces = cfeat.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')
color = ['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']
level = np.arange(0,10,1)

ax1.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax1.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax1.xaxis.set_major_formatter(LongitudeFormatter())
ax1.yaxis.set_major_formatter(LatitudeFormatter())
ax1.grid(c='k', ls='--', alpha=0.3)
ax1.add_feature(cfeat.BORDERS)
ax1.add_feature(states_provinces, edgecolor='0.25')
ax1.coastlines()
ax1.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
ax1.set_title('(a) ERA5', loc='left', fontsize=font_size, fontweight='bold')
cf = ax1.contourf(lon_, lat_, data_mean-data_mean, levels=level, cmap=matplotlib.colors.ListedColormap(color), transform=ccrs.PlateCarree())
sc = ax1.scatter(lon_era5_ls, lat_era5_ls, 12, color='blue', edgecolors='black', linewidth=0.5, marker='o') 

ax2.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax2.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax2.xaxis.set_major_formatter(LongitudeFormatter())
ax2.yaxis.set_major_formatter(LatitudeFormatter())
ax2.grid(c='k', ls='--', alpha=0.3)
ax2.add_feature(cfeat.BORDERS)
ax2.add_feature(states_provinces, edgecolor='0.25')
ax2.coastlines()
ax2.set_title('(b) RegCM5', loc='left', fontsize=font_size, fontweight='bold')
cf = ax2.contourf(lon_, lat_, data_mean-data_mean, levels=level, cmap=matplotlib.colors.ListedColormap(color), transform=ccrs.PlateCarree())
sc = ax2.scatter(lon_regcm5_ls, lat_regcm5_ls, 12, color='blue', edgecolors='black', linewidth=0.5, marker='o') 

ax3.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax3.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax3.xaxis.set_major_formatter(LongitudeFormatter())
ax3.yaxis.set_major_formatter(LatitudeFormatter())
ax3.grid(c='k', ls='--', alpha=0.3)
ax3.add_feature(cfeat.BORDERS)
ax3.add_feature(states_provinces, edgecolor='0.25')
ax3.coastlines()
ax3.set_title('(c) WRF415', loc='left', fontsize=font_size, fontweight='bold')
cf = ax3.contourf(lon_, lat_, data_mean-data_mean, levels=level, cmap=matplotlib.colors.ListedColormap(color), transform=ccrs.PlateCarree())
sc = ax3.scatter(lon_wrf415_ls, lat_wrf415_ls, 12, color='blue', edgecolors='black', linewidth=0.5, marker='o') 

# Path out to save figure
path_out = '{0}/figs/cyclone'.format(path)
name_out = 'pyplt_maps_tracking_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
