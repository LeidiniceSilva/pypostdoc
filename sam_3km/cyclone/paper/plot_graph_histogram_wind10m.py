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

font_size = 8
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
				lat_cyc.append(float(k[1][:]))
				lon_cyc.append(float(k[2][:]))
				dt_cyc.append(str(k[0][:]))

	return lat_cyc, lon_cyc, dt_cyc


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
	
	mean_ = np.ma.stack(var_, axis=0)

	return lat, lon, mean_


def extract_timeseries(dataset, lat_array, lon_array, target_lats, target_lons):

	values = []
	for t in range(len(target_lats)):
		lat_idx = np.abs(lat_array - target_lats[t]).argmin()
		lon_idx = np.abs(lon_array - target_lons[t]).argmin()
		val = dataset[t, lat_idx, lon_idx]
		values.append(val)

	return values


def comp_stats(wind_dataset):

	wind_median_ = np.nanmedian(wind_dataset)
	wind_mean_ = np.nanmean(wind_dataset)
	wind_p90_ = np.nanpercentile(wind_dataset, 90)
	wind_max_ = np.nanmax(wind_dataset)

	wind_median = round(wind_median_, 2)
	wind_mean = round(wind_mean_, 2)
	wind_p90 = round(wind_p90_, 2)
	wind_max = round(wind_max_, 2)

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
wrf415_wind_median, wrf415_wind_mean, wrf415_wind_p90, wrf415_wind_max = comp_stats(uv10_ts_wrf415)

print(era5_wind_median, era5_wind_mean, era5_wind_p90, era5_wind_max)
print(regcm5_wind_median, regcm5_wind_mean, regcm5_wind_p90,regcm5_wind_max)
print(wrf415_wind_median, wrf415_wind_mean, wrf415_wind_p90, wrf415_wind_max)

# Plot figure
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 4))

ax1.set_facecolor('lightgray')
ax1.hist(uv10_ts_era5, bins=10, color='black', alpha=0.75, edgecolor='black', label='ERA5')
ax1.axvline(era5_wind_median, color='gray', label='Median: {0}'.format(era5_wind_median))
ax1.axvline(era5_wind_mean, color='green', label='Mean: {0}'.format(era5_wind_mean))
ax1.axvline(era5_wind_p90, color='orange', label='90th percentile: {0}'.format(era5_wind_p90))
ax1.axvline(era5_wind_max, color='purple', label='Max: {0}'.format(era5_wind_max))
ax1.set_xlabel('Wind speed at 10 m (m s$^-$$^1$)', fontsize=font_size, fontweight='bold')
ax1.set_title('(a)', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_ylim(0, 80)
ax1.set_yticks(np.arange(0, 90, 10))
ax1.set_xlim(0, 20)
ax1.set_xticks(np.arange(0, 22, 2))
ax1.legend(fontsize=font_size, ncol=1, loc=2, shadow=True)
ax1.grid(True, color='k', linestyle='--', alpha=0.5)

ax2.set_facecolor('lightgray')
ax2.hist(uv10_ts_regcm5, bins=10, color='blue', alpha=0.75, edgecolor='black', label='RegCM5')
ax2.axvline(regcm5_wind_median, color='gray', label='Median: {0}'.format(regcm5_wind_median))
ax2.axvline(regcm5_wind_mean, color='green', label='Mean: {0}'.format(regcm5_wind_mean))
ax2.axvline(regcm5_wind_p90, color='orange', label='90th percentile: {0}'.format(regcm5_wind_p90))
ax2.axvline(regcm5_wind_max, color='purple', label='Max: {0}'.format(regcm5_wind_max))
ax2.set_xlabel('Wind speed at 10 m (m s$^-$$^1$)', fontsize=font_size, fontweight='bold')
ax2.set_title('(b)', loc='left', fontsize=font_size, fontweight='bold')
ax2.set_ylim(0, 80)
ax2.set_yticks(np.arange(0, 90, 10))
ax2.set_xlim(0, 20)
ax2.set_xticks(np.arange(0, 22, 2))
ax2.legend(fontsize=font_size, ncol=1, loc=2, shadow=True)
ax2.grid(True, color='k', linestyle='--', alpha=0.5)

ax3.set_facecolor('lightgray')
ax3.hist(uv10_ts_wrf415, bins=10, color='red', alpha=0.75, edgecolor='black', label='WRF415')
ax3.axvline(wrf415_wind_median, color='gray', label='Median: {0}'.format(wrf415_wind_median))
ax3.axvline(wrf415_wind_mean, color='green', label='Mean: {0}'.format(wrf415_wind_mean))
ax3.axvline(wrf415_wind_p90, color='orange', label='90th percentile: {0}'.format(wrf415_wind_p90))
ax3.axvline(wrf415_wind_max, color='purple', label='Max: {0}'.format(wrf415_wind_max))
ax3.set_xlabel('Wind speed at 10 m (m s$^-$$^1$)', fontsize=font_size, fontweight='bold')
ax3.set_title('(c)', loc='left', fontsize=font_size, fontweight='bold')
ax3.set_ylim(0, 80)
ax3.set_yticks(np.arange(0, 90, 10))
ax3.set_xlim(0, 20)
ax3.set_xticks(np.arange(0, 22, 2))
ax3.legend(fontsize=font_size, ncol=1, loc=2, shadow=True)
ax3.grid(True, color='k', linestyle='--', alpha=0.5)

# Path out to save figure
path_out = '{0}/SAM-3km/figs/cyclone'.format(path)
name_out = 'pyplt_graph_histogram_wind_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

