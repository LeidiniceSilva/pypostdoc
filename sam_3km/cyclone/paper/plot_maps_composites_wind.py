# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot map of windy speed"

import os
import netCDF4
import datetime
import numpy as np
import matplotlib.cm as cm
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeat

from scipy import signal, misc
from datetime import datetime
from datetime import datetime, timedelta
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

font_size = 10
path='/marconi/home/userexternal/mdasilva'


def generate_hourly_dates(start_date, end_date):
    
	dates = []
	current_date = start_date
	while current_date <= end_date:
		dates.append(current_date.strftime('%Y%m%d%H'))
		current_date += timedelta(hours=6)
	
	return dates


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


def open_dat_file(dataset):

	dt = []
	for yr in range(2018, 2021+1):
	
		data = read_dat_file('{0}/user/mdasilva/SAM-3km/post_cyclone/ECyclone_v2/{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
	
		rows_list = []
		rows_list_i = []
		for i, (header, rows) in enumerate(data):
			if (rows[0][2] < str(-56)):
				rows_list.append(rows[0])
		
		for j  in rows_list:
			dt.append(str(j[0]))
	
	return dt


def import_data(param, dataset, indices):

	if dataset == 'ERA5':
		arq = '{0}/user/mdasilva/SAM-3km/post_cyclone/obs/era5/era5/{1}_SAM-25km_{2}_6hr_2018-2021_lonlat.nc'.format(path, param, dataset) 	 
	elif dataset == 'RegCM5':
		arq = '{0}/user/mdasilva/SAM-3km/post_cyclone/regcm5/regcm5/{1}/{1}_SAM-3km_{2}_6hr_2018-2021_lonlat.nc'.format(path, param, dataset)
	else:
		arq = '{0}/user/mdasilva/SAM-3km/post_cyclone/wrf/wrf/{1}/{1}_SAM-3km_{2}_6hr_2018-2021_lonlat.nc'.format(path, param, dataset)
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]	
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]

	var_i = []
	for idx_i in indices:
		var_i.append(np.squeeze(mean[idx_i,:,:]))
		
	return lat, lon, var_i
	

# Generate list of daily dates from 2018 to 2021
hourly_dates = generate_hourly_dates(datetime(2018, 1, 1, 0), datetime(2021, 12, 31, 23))
print("here")

# Import cyclone tracking date 
dt_era5 = open_dat_file('ERA5')
dt_regcm5 = open_dat_file('RegCM5')
dt_wrf415 = open_dat_file('WRF415')
print("here")

era5_idx_i = find_indices_in_date_list(hourly_dates, dt_era5)
regcm5_idx_i = find_indices_in_date_list(hourly_dates, dt_regcm5)
wrf415_idx_i = find_indices_in_date_list(hourly_dates, dt_wrf415)
print("here")

# Import model and obs dataset 
lat, lon, era5_u10 = import_data('u10', 'ERA5', era5_idx_i)
lat, lon, era5_v10 = import_data('v10', 'ERA5', era5_idx_i)
print("here")

lat, lon, regcm5_uas = import_data('uas', 'RegCM5', regcm5_idx_i)
lat, lon, regcm5_vas = import_data('vas', 'RegCM5', regcm5_idx_i)
print("here")

lat, lon, wrf415_U10e = import_data('U10e', 'WRF415', wrf415_idx_i)
lat, lon, wrf415_V10e = import_data('V10e', 'WRF415', wrf415_idx_i)
print("here")

era5_u10_i = np.nanmean(era5_u10, axis=0)
era5_v10_i = np.nanmean(era5_v10, axis=0)
uv_era5 = np.sqrt(era5_u10_i**2 + era5_v10_i**2)
print("here")

regcm5_uas_i = np.nanmean(regcm5_uas, axis=0)
regcm5_vas_i = np.nanmean(regcm5_vas, axis=0)
uv_regcm5 = np.sqrt(regcm5_uas_i**2 + regcm5_vas_i**2)
print("here")

wrf415_U10e_i = np.nanmean(wrf415_U10e, axis=0)
wrf415_V10e_i = np.nanmean(wrf415_V10e, axis=0)
uv_wrf415 = np.sqrt(wrf415_U10e_i**2 + wrf415_V10e_i**2)
print("here")

# Plot figure
fig, axes = plt.subplots(2,2, figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2), (ax3, ax4) = axes
fig.delaxes(ax4)

colorb = np.arange(0,11,1)
states_provinces = cfeat.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')

ax1.set_xticks(np.arange(-76,38.5,5), crs=ccrs.PlateCarree())
ax1.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax1.xaxis.set_major_formatter(LongitudeFormatter())
ax1.yaxis.set_major_formatter(LatitudeFormatter())
ax1.grid(c='k', ls='--', alpha=0.3)
ax1.add_feature(cfeat.BORDERS)
ax1.add_feature(states_provinces, edgecolor='0.25')
ax1.coastlines()
ax1.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax1.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
ax1.set_title('(a) ERA5', loc='left', fontsize=font_size, fontweight='bold')
cf = ax1.contourf(lon, lat, uv_era5, colorb, transform=ccrs.PlateCarree(), extend='max', cmap='Blues')
qr = ax1.quiver(lon, lat, era5_u10_i, era5_v10_i, color='black', transform=ccrs.PlateCarree())
plt.quiverkey(qr, X=0.9, Y=1.05, U=5, label='5 m/s', labelpos='E', fontproperties={'size': 8})

ax2.set_xticks(np.arange(-76,38.5,5), crs=ccrs.PlateCarree())
ax2.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax2.xaxis.set_major_formatter(LongitudeFormatter())
ax2.yaxis.set_major_formatter(LatitudeFormatter())
ax2.grid(c='k', ls='--', alpha=0.3)
ax2.add_feature(cfeat.BORDERS)
ax2.add_feature(states_provinces, edgecolor='0.25')
ax2.coastlines()
ax2.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax2.set_title('(b) RegCM5', loc='left', fontsize=font_size, fontweight='bold')
cf = ax2.contourf(lon, lat, uv_regcm5, colorb, transform=ccrs.PlateCarree(), extend='max', cmap='Blues')
cb = plt.colorbar(cf, cax=fig.add_axes([0.92, 0.3, 0.015, 0.4]), ticks=colorb)
qr = ax2.quiver(lon, lat, regcm5_uas_i, regcm5_vas_i, color='black', transform=ccrs.PlateCarree())

ax3.set_xticks(np.arange(-76,38.5,5), crs=ccrs.PlateCarree())
ax3.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax3.xaxis.set_major_formatter(LongitudeFormatter())
ax3.yaxis.set_major_formatter(LatitudeFormatter())
ax3.grid(c='k', ls='--', alpha=0.3)
ax3.add_feature(cfeat.BORDERS)
ax3.add_feature(states_provinces, edgecolor='0.25')
ax3.coastlines()
ax3.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax3.set_title('(c) WRF415', loc='left', fontsize=font_size, fontweight='bold')
cf = ax3.contourf(lon, lat, uv_wrf415, colorb, transform=ccrs.PlateCarree(), extend='max', cmap='Blues')
cb = plt.colorbar(cf, cax=fig.add_axes([0.92, 0.3, 0.015, 0.4]), ticks=colorb)
qr = ax3.quiver(lon, lat, wrf415_U10e_i, wrf415_V10e_i, color='black', transform=ccrs.PlateCarree())

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cyclone/paper'.format(path)
name_out = 'pyplt_maps_wind_speed_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
exit()
