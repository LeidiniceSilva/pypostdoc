# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot composites"

import os
import netCDF4
import numpy as np
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeat

from datetime import datetime, timedelta
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

path='/marconi/home/userexternal/mdasilva'


def generate_hourly_dates(start_date, end_date):
    
	dates_hr = []
	current_date = start_date
	while current_date <= end_date:
		dates_hr.append(current_date.strftime('%Y%m%d%H'))
		current_date += timedelta(hours=6)
	
	return dates_hr
	

def find_indices_in_date_list(date_list, target_dates):
    
	indices = []
	for target_date in target_dates:
		try:
			index = date_list.index(target_date)
			indices.append(index)
		except ValueError:
			pass  # Date not found in date_list
	
	return indices
	
	
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
			
			# If we have rows, append them to data
			if rows:
				data.append((header, rows))
				rows = []
			header = line
		else:  
			rows.append(line)
	
	# Append the last header and rows to data
	if header and rows:
		data.append((header, rows))
	
	return data


def open_file_dt(dataset):

	dt_hr = []
	for yr in range(2018, 2021+1):
	
		data = read_dat_file('{0}/user/mdasilva/SAM-3km/post_cyclone/ECv2/{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
	
		for i, (header, rows) in enumerate(data):
			if (rows[0][2] < str(-56)):
				dt_hr.append(header[1][:])
				
	return dt_hr
	
	
def import_data(param, dataset, indices):

	if dataset == 'RegCM5':
		arq = '{0}/user/mdasilva/SAM-3km/post_cyclone/regcm5/regcm5/{1}_SAM-3km_{2}_6hr_2018-2021_lonlat.nc'.format(path, param, dataset)
	elif dataset == 'WRF415':
		arq = '{0}/user/mdasilva/SAM-3km/post_cyclone/wrf/wrf/{1}_SAM-3km_{2}_6hr_2018-2021_lonlat.nc'.format(path, param, dataset)
	else:
		arq   = '{0}/user/mdasilva/SAM-3km/post_cyclone/obs/era5/era5/{1}_SAM-25km_{2}_6hr_2018-2021_lonlat.nc'.format(path, param, dataset)		
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]	
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]

	var_i = []
	for idx in indices:
		var_i.append(np.squeeze(mean[idx,:,:]))
	
	mean_i = np.nanmean(var_i, axis=0)
		
	return lat, lon, mean_i
	
	
# Generate list of dates from 2018 to 2021
hourly_dates = generate_hourly_dates(datetime(2018, 1, 1, 0), datetime(2021, 12, 31, 23))

# Import cyclone tracking date 
dt_era5 = open_file_dt('ERA5')
dt_regcm5 = open_file_dt('RegCM5')
dt_wrf415 = open_file_dt('WRF415')

era5_idx = find_indices_in_date_list(hourly_dates, dt_era5)
regcm5_idx = find_indices_in_date_list(hourly_dates, dt_regcm5)
wrf415_idx = find_indices_in_date_list(hourly_dates, dt_wrf415)

# Import model and obs dataset 
lat, lon, msl_era5 = import_data('msl', 'ERA5', era5_idx)
lat, lon, u10_era5 = import_data('u10', 'ERA5', era5_idx)
lat, lon, v10_era5 = import_data('v10', 'ERA5', era5_idx)

lat, lon, msl_regcm5 = import_data('psl', 'RegCM5', regcm5_idx)
lat, lon, u10_regcm5 = import_data('uas', 'RegCM5', regcm5_idx)
lat, lon, v10_regcm5 = import_data('vas', 'RegCM5', regcm5_idx)

lat, lon, msl_wrf415 = import_data('PSL', 'WRF415', wrf415_idx)
lat, lon, u10_wrf415 = import_data('U10e', 'WRF415', wrf415_idx)
lat, lon, v10_wrf415 = import_data('V10e', 'WRF415', wrf415_idx)

# Calculate wind speed
uv10_era5 = np.sqrt(u10_era5**2 + v10_era5**2)
uv10_regcm5 = np.sqrt(u10_regcm5**2 + v10_regcm5**2)
uv10_wrf415 = np.sqrt(u10_wrf415**2 + v10_wrf415**2)

# Plot figure
fig, axes = plt.subplots(1,3, figsize=(14, 3), subplot_kw={"projection": ccrs.PlateCarree()})
ax1, ax2, ax3 = axes
font_size = 10

states_provinces = cfeat.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')
level = np.arange(0,12.5,0.5)

ax1.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax1.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax1.xaxis.set_major_formatter(LongitudeFormatter())
ax1.yaxis.set_major_formatter(LatitudeFormatter())
ax1.grid(c='k', ls='--', alpha=0.3)
ax1.add_feature(cfeat.BORDERS)
ax1.add_feature(states_provinces, edgecolor='0.25')
ax1.coastlines()
ax1.set_title('(a) ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax1.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
cf = ax1.contourf(lon, lat, uv10_era5, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap='jet')
ct = ax1.contour(lon, lat, msl_era5/100, colors='black', linewidths=0.50)
ax1.clabel(ct, inline=1, fontsize=font_size)

ax2.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax2.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax2.xaxis.set_major_formatter(LongitudeFormatter())
ax2.yaxis.set_major_formatter(LatitudeFormatter())
ax2.grid(c='k', ls='--', alpha=0.3)
ax2.add_feature(cfeat.BORDERS)
ax2.add_feature(states_provinces, edgecolor='0.25')
ax2.coastlines()
ax2.set_title('(b) RegCM5', loc='left', fontsize=font_size, fontweight='bold')
ax2.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
cf = ax2.contourf(lon, lat, uv10_regcm5, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap='jet')
ct = ax2.contour(lon, lat, msl_regcm5/100, colors='black', linewidths=0.50)
ax2.clabel(ct, inline=1, fontsize=font_size)

ax3.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax3.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax3.xaxis.set_major_formatter(LongitudeFormatter())
ax3.yaxis.set_major_formatter(LatitudeFormatter())
ax3.grid(c='k', ls='--', alpha=0.3)
ax3.add_feature(cfeat.BORDERS)
ax3.add_feature(states_provinces, edgecolor='0.25')
ax3.coastlines()
ax3.set_title('(c) WRF415', loc='left', fontsize=font_size, fontweight='bold')
ax3.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
cf = ax3.contourf(lon, lat, uv10_wrf415, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap='jet')
ct = ax3.contour(lon, lat, msl_wrf415, colors='black', linewidths=0.50)
ax3.clabel(ct, inline=1, fontsize=font_size)
cb = plt.colorbar(cf, cax=fig.add_axes([0.91, 0.2, 0.015, 0.6]))

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cyclone/paper'.format(path)
name_out = 'pyplt_maps_msl_uv10_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
