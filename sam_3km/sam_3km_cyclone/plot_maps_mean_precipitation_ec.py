# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot map of mean precipitation"

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
        dates.append(current_date.strftime('%Y%m%d'))
        current_date += timedelta(days=1)
    
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
			rows_list.append(rows)
		
		for j  in rows_list:
			for k in j:
				dt.append(str(k[0][:-2]))

	return dt


def import_obs(param):

	arq   = '{0}/user/mdasilva/SAM-3km/post_cyclone/obs/era5/{1}_SAM-25km_ERA5_1hr_2018010100-2021123100.nc'.format(path, param)		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['latitude'][:]	
	lon   = data.variables['longitude'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	

def import_rcm(param):

	arq   = '{0}/user/mdasilva/SAM-3km/post_cyclone/rcm/regcm5/{1}/{1}_SAM-3km_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5_6hr_20180101-20211201_lonlat.nc'.format(path, param)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]	
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	

# Generate list of daily dates from 2018 to 2019
daily_dates = generate_daily_dates(datetime(2018, 1, 1), datetime(2021, 12, 31))

# Import model and obs dataset 
lat, lon, pr_era5 = import_obs('tp')

dt_era5 = open_dat_file('ERA5')
era5_indices = remove_duplicates(dt_era5)
era5_idx = find_indices_in_date_list(daily_dates, era5_indices)

pr_era5_i = []
for idx_i in era5_idx:
	pr_era5_i.append(pr_era5[idx_i,:,:])

print(pr_era5_i)
exit()

pr_era5_ii = np.sum(pr_era5_i*86400, axis=0)

# Plot figure
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12, 8), subplot_kw={"projection": ccrs.PlateCarree()})

color = ['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']

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
cf = ax1.contourf(lon, lat, pr_era5_ii/4, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))

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
cf = ax2.contourf(lon, lat, pr_era5_ii/4, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))
cb = plt.colorbar(cf, cax=fig.add_axes([0.92, 0.3, 0.015, 0.4]))

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cyclone'.format(path)
name_out = 'pyplt_maps_wind_speed_EC_ERA5_RegCM5_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
