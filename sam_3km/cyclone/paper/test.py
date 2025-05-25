# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot wind histogram"

import os
import netCDF4
import datetime
import numpy as np
import xarray as xr
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from scipy import signal, misc
from datetime import datetime, timedelta

font_size = 10
path = '/leonardo/home/userexternal/mdasilva/leonardo_work'


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
			pass  
	
	return indices
	
	
def read_dat_file(filename):

	data = []
	with open(filename, 'r') as file:
		lines = file.readlines()
		header = []
		rows = []
	
	for line in lines:
		line = line.strip().split()
		if len(line) == 6:
			if rows:
				data.append((header, rows))
				rows = []
			header = line
		else:  
			rows.append(line)
	
	if header and rows:
		data.append((header, rows))
	
	return data


def open_file_dt(dataset):

	lat_cyc, lon_cyc, dt_cyc = [], [], []
	for yr in range(2018, 2021+1):
	
		data = read_dat_file('{0}/SAM-3km/postproc/cyclone/{1}/track/resultado_{2}.dat'.format(path, dataset, yr))

		rows_list = []
		for i, (header, rows) in enumerate(data):
			if (rows[0][2] < str(-56)):
				rows_list.append(rows)
	
		for j  in rows_list:
			for k in j:
				lat_cyc.append(k[1][:])
				lon_cyc.append(k[2][:])
				dt_cyc.append(str(k[0][:]))

	return lat_cyc, lon_cyc, dt_hr


def import_data(param, dataset, indices):

	if dataset == 'RegCM5':
		arq = '{0}/SAM-3km/postproc/cyclone/RegCM5/{1}_SAM-3km_{2}_6hr_2018-2021_lonlat.nc'.format(path, param, dataset)
	elif dataset == 'WRF415':
		arq = '{0}/SAM-3km/postproc/cyclone/WRF415/{1}_SAM-3km_{2}_6hr_2018-2021_lonlat.nc'.format(path, param, dataset)	
	else:
		arq   = '{0}/SAM-3km/postproc/cyclone/ERA5/{1}_SAM-3km_{2}_6hr_2018-2021_lonlat.nc'.format(path, param, dataset)
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]	
	lon   = data.variables['lon'][:]

	if dataset == 'RegCM5':
		mean = var[:][:,0,:,:]
	else:
		mean = var[:][:,:,:]

	var_ = []
	for idx_i in indices:
		var_.append(np.squeeze(mean[idx_i,:,:]))

	return lat, lon, var_


def extract_timeseries(dataset, lat_array, lon_array, target_lat, target_lon):

	lat_idx = np.abs(lat_array - target_lat).argmin()
	lon_idx = np.abs(lon_array - target_lon).argmin()
	timeseries = dataset[:, lat_idx, lon_idx]

	return timeseries


def comp_stats(wind_dataset):

	wind_median = np.nanmedian(wind_dataset)
	wind_mean = np.nanmean(wind_dataset)
	wind_p90 = np.nanpercentile(wind_dataset, 90)
	wind_max = np.nanmax(wind_dataset)

	return wind_median, wind_mean, wind_p90, wind_max


# Generate list of dates from 2018 to 2021
hourly_dates = generate_hourly_dates(datetime(2018, 1, 1, 0), datetime(2021, 12, 31, 23))

# Import cyclone tracking info 
lat_era5, lon_era5, dt_era5 = open_file_dt('ERA5')
lat_regcm5, lon_regcm5, dt_regcm5 = open_file_dt('RegCM5')
lat_wrf415, lon_wrf415, dt_wrf415 = open_file_dt('WRF415')

era5_idx = find_indices_in_date_list(hourly_dates, dt_era5)
regcm5_idx = find_indices_in_date_list(hourly_dates, dt_regcm5)
wrf415_idx = find_indices_in_date_list(hourly_dates, dt_wrf415)

# Import model and obs dataset 
lat, lon, era5_u10 = import_data('u10', 'ERA5', era5_idx)
lat, lon, era5_v10 = import_data('v10', 'ERA5', era5_idx)
lat, lon, regcm5_uas = import_data('uas', 'RegCM5', regcm5_idx)
lat, lon, regcm5_vas = import_data('vas', 'RegCM5', regcm5_idx)
lat, lon, wrf415_u10 = import_data('U10', 'WRF415', wrf415_idx)
lat, lon, wrf415_v10 = import_data('V10', 'WRF415', wrf415_idx)

print(era5_u10.shape)
print(era5_v10.shape)
print(regcm5_uas.shape)
print(regcm5_vas.shape)
print(wrf415_u10.shape)
print(wrf415_v10.shape)

# Calculate wind speed
uv10_era5 = np.sqrt(era5_u10**2 + era5_v10**2)
uv10_regcm5 = np.sqrt(regcm5_uas**2 + regcm5_vas**2)
uv10_wrf415 = np.sqrt(wrf415_u10**2 + wrf415_v10**2)

print(uv10_era5.shape)
print(uv10_regcm5.shape)
print(uv10_wrf415.shape)

uv10_ts_era5 = extract_timeseries(uv10_era5, lat, lon, lat_era5, lon_era5)
uv10_ts_regcm5 = extract_timeseries(uv10_regcm5, lat, lon, lat_regcm5, lon_regcm5)
uv10_ts_wrf415 = extract_timeseries(uv10_wrf415, lat, lon, lat_wrf415, lon_wrf415)

print(uv10_ts_era5)
print(uv10_ts_regcm5)
print(uv10_ts_regcm5)

era5_wind_median, era5_wind_mean, era5_wind_p90, era5_wind_max = comp_stats(uv10_ts_era5)
regcm5_wind_median, regcm5_wind_mean, regcm5_wind_p90,regcm5_wind_max = comp_stats(uv10_ts_regcm5)
wrf415_wind_median, wrf415_wind_mean, wrf415_wind_p90, wrf415_wind_max = comp_stats(uv10_ts_regcm5)

print(era5_wind_median, era5_wind_mean, era5_wind_p90, era5_wind_max)
print(regcm5_wind_median, regcm5_wind_mean, regcm5_wind_p90,regcm5_wind_max)
print(wrf415_wind_median, wrf415_wind_mean, wrf415_wind_p90, wrf415_wind_max)

# Plot figure
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))
font_size = 10

ax1.hist(uv10_ts_era5, bins=40, color=color, edgecolor='black', label='ERA5')
ax1.axvline(era5_wind_median, color='red', label=f'Median: {median_val:.1f}')
ax1.axvline(era5_wind_mean, color='darkorange', label=f'Mean: {mean_val:.1f}')
ax1.axvline(era5_wind_p90, color='blue', label=f'90th percentile: {p90_val:.1f}')
ax1.axvline(era5_wind_max, color='purple', label=f'Max: {max_val:.1f}')
ax1.set_xlabel('Wind speed at 10 m (m s$^-$$^1$)', fontsize=font_size)
ax1.set_ylabel('Frequency', fontsize=font_size)
ax1.set_title('(a)', fontsize=font_size)
ax1.legend()
ax1.grid(True, linestyle='--', alpha=0.5)

ax2.hist(uv10_ts_regcm5, bins=40, color=color, edgecolor='black', label='RegCM5')
ax2.axvline(regcm5_wind_median, color='red', label=f'Median: {median_val:.1f}')
ax2.axvline(regcm5_wind_mean, color='darkorange', label=f'Mean: {mean_val:.1f}')
ax2.axvline(regcm5_wind_p90, color='blue', label=f'90th percentile: {p90_val:.1f}')
ax2.axvline(regcm5_wind_max, color='purple', label=f'Max: {max_val:.1f}')
ax1.set_xlabel('Wind speed at 10 m (m s$^-$$^1$)', fontsize=font_size)
ax1.set_ylabel('Frequency', fontsize=font_size)
ax1.set_title('(a)', fontsize=font_size)
ax1.legend()
ax1.grid(True, linestyle='--', alpha=0.5)

ax3.hist(uv10_ts_wrf415, bins=40, color=color, edgecolor='black', label='WRF415')
ax3.axvline(wrf415_wind_median, color='red', label=f'Median: {median_val:.1f}')
ax3.axvline(wrf415_wind_mean, color='darkorange', label=f'Mean: {mean_val:.1f}')
ax3.axvline(wrf415_wind_p90, color='blue', label=f'90th percentile: {p90_val:.1f}')
ax3.axvline(wrf415_wind_max, color='purple', label=f'Max: {max_val:.1f}')
ax1.set_xlabel('Wind speed at 10 m (m s$^-$$^1$)', fontsize=font_size)
ax1.set_ylabel('Frequency', fontsize=font_size)
ax1.set_title('(a)', fontsize=font_size)
ax1.legend()
ax1.grid(True, linestyle='--', alpha=0.5)

# Path out to save figure
path_out = '{0}/figs/cyclone'.format(path)
name_out = 'pyplt_graph_histogram_wind_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

