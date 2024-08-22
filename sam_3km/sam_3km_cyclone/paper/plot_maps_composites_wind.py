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

path='/marconi/home/userexternal/mdasilva'


def remove_duplicates(date_list):
	
	unique_dates = []
	for date in date_list:
		if date not in unique_dates:
			unique_dates.append(date)
	
	return unique_dates
    

def generate_daily_dates(start_date, end_date):
    
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
	
		data = read_dat_file('{0}/user/mdasilva/SAM-3km/post_cyclone/ECyclone/ECyclone_{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
	
		rows_list = []
		rows_list_i = []
		for i, (header, rows) in enumerate(data):
			if ((str(-20) < rows[0][1] < str(-30)) and (str(-45) < rows[0][2] < str(-55))):
				rows_list.append(rows[0])
		
		for j  in rows_list:
			dt.append(str(j[0]))
	
	return dt


def import_obs(param):

	arq   = '{0}/user/mdasilva/SAM-3km/post_cyclone/obs/postproc/{1}_SAM-25km_ERA5_1hr_20180101-20211201_lonlat.nc'.format(path, param)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]	
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	

def import_rcm(param):

	arq   = '{0}/user/mdasilva/SAM-3km/post_cyclone/rcm/postproc/{1}_SAM-3km_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5_6hr_20180101-20211231.nc'.format(path, param)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]	
	lon   = data.variables['lon'][:]
	mean = var[:][:,0,:,:]
	
	return lat, lon, mean

# Generate list of daily dates from 2018 to 2021
daily_dates = generate_daily_dates(datetime(2018, 1, 1), datetime(2021, 12, 31))

# Import model and obs dataset 
dt_era5 = open_dat_file('ERA5')
dt_regcm5 = open_dat_file('RegCM5')

lat, lon, ua_era5 = import_obs('u')
lat, lon, va_era5 = import_obs('v')

lat, lon, ua_regcm5 = import_rcm('ua')
lat, lon, va_regcm5 = import_rcm('va')

era5_indices = find_indices_in_date_list(dt_era5, daily_dates)
regcm5_indices = find_indices_in_date_list(dt_regcm5, daily_dates)

# Import indices after tracking
era5_idx = remove_duplicates(dt_era5)
era5_idx_i = find_indices_in_date_list(daily_dates, era5_idx)

regcm5_idx = remove_duplicates(dt_regcm5)
regcm5_idx_i = find_indices_in_date_list(daily_dates, regcm5_idx)

ua_era5_i = []
va_era5_i = []
for idx_i in era5_idx_i:
	ua_era5_i.append(ua_era5[idx_i-1,:,:])
	va_era5_i.append(va_era5[idx_i-1,:,:])

ua_era5_ii = np.nanmean(ua_era5_i, axis=0)
va_era5_ii = np.nanmean(va_era5_i, axis=0)
ws_era5 = np.sqrt(ua_era5_ii**2 + va_era5_ii**2)

ua_regcm5_i = []
va_regcm5_i = []
for idx_ii in regcm5_idx_i:
	ua_regcm5_i.append(ua_regcm5[idx_ii-1,:,:])
	va_regcm5_i.append(va_regcm5[idx_ii-1,:,:])

ua_regcm5_ii = np.nanmean(ua_regcm5_i, axis=0)
va_regcm5_ii = np.nanmean(va_regcm5_i, axis=0)
ws_regcm5 = np.sqrt(ua_regcm5_ii**2 + va_regcm5_ii**2)

# Plot figure
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12, 8), subplot_kw={"projection": ccrs.PlateCarree()})
font_size = 10

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
ax1.set_title('a) ERA5', loc='left', fontsize=font_size, fontweight='bold')
cf = ax1.contourf(lon,lat, ws_era5, colorb, transform=ccrs.PlateCarree(), extend='max', cmap='Blues')
qr = ax1.quiver(lon, lat, ua_era5_ii, va_era5_ii, color='black', transform=ccrs.PlateCarree())
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
ax2.set_title('b) RegCM5', loc='left', fontsize=font_size, fontweight='bold')
cf = ax2.contourf(lon, lat, ws_regcm5, colorb, transform=ccrs.PlateCarree(), extend='max', cmap='Blues')
cb = plt.colorbar(cf, cax=fig.add_axes([0.92, 0.3, 0.015, 0.4]), ticks=colorb)
qr = ax2.quiver(lon, lat, ua_regcm5_ii, va_regcm5_ii, color='black', transform=ccrs.PlateCarree())

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cyclone'.format(path)
name_out = 'pyplt_maps_wind_speed_EC_ERA5_RegCM5_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
