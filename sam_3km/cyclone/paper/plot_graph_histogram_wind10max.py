# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot map of precipitation"

import os
import netCDF4
import datetime
import numpy as np
import xarray as xr
import matplotlib.colors
import matplotlib.cm as cm
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeat

from scipy import signal, misc
from datetime import datetime, timedelta
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

font_size = 10
path = '/leonardo/home/userexternal/mdasilva/leonardo_work'


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
		
		if len(line) == 6:
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
	
		data = read_dat_file('{0}/SAM-3km/postproc/cyclone/{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
	
		rows_list = []
		rows_list_i = []
		for i, (header, rows) in enumerate(data):
			if (rows[0][2] < str(-56)):
				rows_list.append(rows)
		
		for j  in rows_list:
			for k in j:
				dt.append(str(k[0][:-2]))

	return dt

	
def import_data(param, dataset, indices):

	if dataset == 'RegCM5':
		arq = '{0}/SAM-3km/postproc/cyclone/RegCM5/{1}_SAM-3km_{2}_day_2018-2021_lonlat.nc'.format(path, param, dataset)
	elif dataset == 'WRF415':
		arq = '{0}/SAM-3km/postproc/cyclone/WRF415/{1}_SAM-3km_{2}_day_2018-2021_lonlat.nc'.format(path, param, dataset)	
	else:
		arq   = '{0}/SAM-3km/postproc/cyclone/ERA5/{1}_SAM-3km_{2}_day_2018-2021_lonlat.nc'.format(path, param, dataset)
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]	
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
		
	return lat, lon, mean

	
def get_max_wind_near_cyclone(lat, lon, wind_data, dataset):
	
	wind_speed = []
	for yr in range(2018, 2021+1):
		data = read_dat_file('{0}/postproc/cyclone/{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
		for i, (header, rows) in enumerate(data):
			lat_cyc = float(row[0][1])
			lon_cyc = float(row[0][1])                   
			lat_mask = (lat >= lat_cyc - 1) & (lat <= lat_cyc + 1)
			lon_mask = (lon >= lon_cyc - 1) & (lon <= lon_cyc + 1)

			if not np.any(lat_mask) or not np.any(lon_mask):
                        	continue

			gust_masked = wind_data[:, lat_mask, lon_mask)]
			wind_speed.append(gust_masked.compressed() if hasattr(gust_masked, 'mask') else gust_masked.flatten())

	return wind_speed


def comp_stats(wind_dataset):

	max_winds = np.max(wind_dataset, axis=0)  
	bin_edges = np.arange(0, 51, 5)
	bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
	counts, _ = np.histogram(max_winds.flatten(), bins=bin_edges)
	normalized_freq = counts / counts.sum()  

	return normalized_freq

# Generate list of daily dates from 2018 to 2021
daily_dates = generate_daily_dates(datetime(2018, 1, 1), datetime(2021, 12, 31))

# Import cyclone tracking date 
dt_era5 = open_dat_file('ERA5')
dt_regcm5 = open_dat_file('RegCM5')
dt_wrf415 = open_dat_file('WRF415')

era5_idx = remove_duplicates(dt_era5)
regcm5_idx = remove_duplicates(dt_regcm5)
wrf415_idx = remove_duplicates(dt_wrf415)

era5_idx_i = find_indices_in_date_list(daily_dates, era5_idx)
regcm5_idx_i = find_indices_in_date_list(daily_dates, regcm5_idx)
wrf415_idx_i = find_indices_in_date_list(daily_dates, wrf415_idx)

# Import model and obs dataset 
lat, lon, era5_u10 = import_data('u10', 'ERA5', era5_idx_i)
lat, lon, era5_v10 = import_data('v10', 'ERA5', era5_idx_i)
lat, lon, uv10_regcm5 = import_data('sfcWindmax', 'RegCM5', regcm5_idx_i)
lat, lon, wrf415_u10 = import_data('U10MAX', 'WRF515', wrf415_idx_i)
lat, lon, wrf415_v10 = import_data('V10MAX', 'WRF515', wrf415_idx_i)

# Calculate wind speed
uv10_era5 = np.sqrt(era5_u10**2 + era5_v10**2)
uv10_regcm5 = np.sqrt(regcm5_uas**2 + regcm5_vas**2)
uv10_wrf415 = np.sqrt(wrf415_u10**2 + wrf415_v10**2)

era5_wind = wind_speed_cyclone(lat, lon, uv10_era5, 'ERA5')
regcm5_wind = wind_speed_cyclone(lat, lon,regcm5_era5, 'RegCM5')
wrf415_wind = wind_speed_cyclone(lat, lon, uv10_wrf415, 'WRF415')

# Plot figure
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))
font_size = 10

ax1.hist(era5_max_wind, bins=40, color=color, edgecolor='black', label='ERA5')
ax1.axvline(era5_max_wind_median, color='red', label=f'Median: {median_val:.1f}')
ax1.axvline(era5_max_wind_mean, color='darkorange', label=f'Mean: {mean_val:.1f}')
ax1.axvline(era5_max_wind_p90, color='blue', label=f'90th percentile: {p90_val:.1f}')
ax1.axvline(era5_max_max, color='purple', label=f'Max: {max_val:.1f}')
ax1.set_xlabel('Max 10 m Wind Gust [m/s]', fontsize=12)
ax1.set_ylabel('Frequency', fontsize=12)
ax1.set_title(label, fontsize=13)
ax1.legend()
ax1.grid(True, linestyle='--', alpha=0.5)

ax2.hist(regcm5_max_wind, bins=40, color=color, edgecolor='black', label='ERA5')
ax2.axvline(regcm5_max_wind_median, color='red', label=f'Median: {median_val:.1f}')
ax2.axvline(regcm5_max_wind_mean, color='darkorange', label=f'Mean: {mean_val:.1f}')
ax2.axvline(regcm5_max_wind_p90, color='blue', label=f'90th percentile: {p90_val:.1f}')
ax2.axvline(regcm5_max_max, color='purple', label=f'Max: {max_val:.1f}')
ax2.set_xlabel('Max 10 m Wind Gust [m/s]', fontsize=12)
ax2.set_ylabel('Frequency', fontsize=12)
ax2.set_title(label, fontsize=13)
ax2.legend()
ax2.grid(True, linestyle='--', alpha=0.5)

ax3.hist(wrf415_max_wind, bins=40, color=color, edgecolor='black', label='ERA5')
ax3.axvline(wrf415_max_wind_median, color='red', label=f'Median: {median_val:.1f}')
ax3.axvline(wrf415_max_wind_mean, color='darkorange', label=f'Mean: {mean_val:.1f}')
ax3.axvline(wrf415_max_wind_p90, color='blue', label=f'90th percentile: {p90_val:.1f}')
ax3.axvline(wrf415_max_max, color='purple', label=f'Max: {max_val:.1f}')
ax3.set_xlabel('Max 10 m Wind Gust [m/s]', fontsize=12)
ax3.set_ylabel('Frequency', fontsize=12)
ax3.set_title(label, fontsize=13)
ax3.legend()
ax3.grid(True, linestyle='--', alpha=0.5)

# Path out to save figure
path_out = '{0}/figs/cyclone'.format(path)
name_out = 'pyplt_graph_histogram_wind_max_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

