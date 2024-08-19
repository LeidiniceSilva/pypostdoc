# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot pdf of precipitation"

import os
import netCDF4
import datetime
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from datetime import datetime, timedelta

path='/marconi/home/userexternal/mdasilva'

skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,77,93,100,105,106,107,112,117,124,135,137,139,
149,152,155,158,168,174,177,183,186,199,204,210,212,224,226,239,240,248,249,253,254,276,277,280,293,298,
303,305,306,308,319,334,335,341,343,359,362,364,384,393,396,398,399,400,402,413,416,417,422,423,426,427,
443,444,446,451,453,457,458,467,474,479,483,488,489,490,495,505,509,513,514,516,529,534,544,559,566]			


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
				rows_list.append(rows)
		
		for j  in rows_list:
			for k in j:
				dt.append(str(k[0][:]))

	return dt


def import_data(param, dataset, indices):

	mean = []
	for station in range(1, 567):
		print(station, inmet[station][0])
		
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue
		
		if dataset == 'RegCM5':
			arq = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/'.format(path) + '{0}_SAM-3km_{1}_1hr_2018-2021_lonlat.nc'.format(param, dataset))
		elif dataset == 'WRF415':
			arq = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_cyclone/wrf/wrf/{1}/'.format(path, param) + '{0}_SAM-3km_{1}_1hr_2018-2021_lonlat.nc'.format(param, dataset))
		else:
			arq = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + '{0}_SAM-3km_{1}_1hr_2018-2021_lonlat.nc'.format(param, dataset))
		
		data   = arq[param]
		latlon = data.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		
		if dataset == 'WRF415':
			time = latlon.sel(XTIME=slice('2018-01-01','2021-12-31'))
		else:
			time = latlon.sel(time=slice('2018-01-01','2021-12-31'))
			
		var = time.values
		var_ = var[::6]

		for idx_i in indices:
			mean.append(var_[idx_i])
																
	return mean		


def import_ws(param, indices):

	mean = []
	for station in range(1, 567):
		print(station, inmet[station][0])
		
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		arq  = xr.open_dataset('{0}/user/mdasilva/WS-SA/INMET/automatic/nc/hourly/{1}/'.format(path, param) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(param, inmet[station][0]))
		data = arq[param]
		time = data.sel(time=slice('2018-01-01','2021-12-31'))
		var  = time.values
		var_ = var[::6]

		for idx_i in indices:
			mean.append(var_[idx_i])
																
	return mean


def comp_pdf(timeseries):

	ts_ = np.array(timeseries)
	ts_list = ts_.flatten()
	ts_round = np.round(ts_list,0)
	ts_filter = ts_round[ts_round > 0.]
	x_pdf_list, pdf_list = np.unique(ts_filter, return_counts=True)
	
	return x_pdf_list, pdf_list
	
		
# Generate list of daily dates from 2018 to 2021
hourly_dates = generate_hourly_dates(datetime(2018, 1, 1, 0), datetime(2021, 12, 31, 23))

# Import cyclone tracking date 
dt_era5 = open_dat_file('ERA5')
dt_regcm5 = open_dat_file('RegCM5')
dt_wrf415 = open_dat_file('WRF415')

era5_idx_i = find_indices_in_date_list(hourly_dates, dt_era5)
regcm5_idx_i = find_indices_in_date_list(hourly_dates, dt_regcm5)
wrf415_idx_i = find_indices_in_date_list(hourly_dates, dt_wrf415)

# Import model and obs dataset 
pr_inmet  = import_ws('pre', era5_idx_i)
pr_gpm    = import_data('precipitation', 'GPM', era5_idx_i)
pr_era5   = import_data('tp', 'ERA5', era5_idx_i)
pr_regcm5 = import_data('pr', 'RegCM5', regcm5_idx_i)
pr_wrf415 = import_data('PREC_ACC_NC', 'WRF415', wrf415_idx_i)

# Import pdf function
x_pdf_inmet, pdf_inmet = comp_pdf(pr_inmet)
x_pdf_gpm, pdf_gpm = comp_pdf(pr_gpm)
x_pdf_era5, pdf_era5 = comp_pdf(pr_era5)
x_pdf_regcm5, pdf_regcm5 = comp_pdf(pr_regcm5)
x_pdf_wrf415, pdf_wrf415 = comp_pdf(pr_wrf415)

# Compute 99.9th percentile
p99_inmet = np.percentile(pr_inmet, 99.9)
p99_gpm = np.percentile(pr_gpm, 99.9)
p99_era5 = np.percentile(pr_era5, 99.9)
p99_regcm5 = np.percentile(pr_regcm5, 99.9)
p99_wrf415 = np.percentile(pr_wrf415, 99.9)

# Plot figure
fig = plt.figure()
font_size = 8

ax = fig.add_subplot(1, 1, 1)  
plt.plot(x_pdf_inmet,  pdf_inmet,  marker='o', markersize=4, mfc='green',  mec='green',  alpha=0.75, linestyle='None', label='INMET')
plt.plot(x_pdf_gpm,    pdf_gpm,    marker='o', markersize=4, mfc='violet', mec='violet', alpha=0.75, linestyle='None', label='GPM')
plt.plot(x_pdf_era5,   pdf_era5,   marker='o', markersize=4, mfc='black',  mec='black',  alpha=0.75, linestyle='None', label='ERA5')
plt.plot(x_pdf_regcm5, pdf_regcm5, marker='o', markersize=4, mfc='blue',   mec='blue',   alpha=0.75, linestyle='None', label='RegCM5')
plt.plot(x_pdf_wrf415, pdf_wrf415, marker='o', markersize=4, mfc='red',    mec='red',    alpha=0.75, linestyle='None', label='WRF415')

plt.axvline(x=p99_inmet, color='green', linestyle='--')
plt.axvline(x=p99_gpm, color='black', linestyle='--')
plt.axvline(x=p99_era5, color='violet', linestyle='--')
plt.axvline(x=p99_regcm5, color='blue', linestyle='--')
plt.axvline(x=p99_wrf415, color='red', linestyle='--')

plt.title('(b)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('Precipitation (mm h$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.grid(axis='y', color='k', linestyle='--', alpha=0.3)
plt.legend(loc=1, ncol=1, fontsize=font_size, shadow=True)

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cyclone/paper'.format(path)
name_out = 'pyplt_graph_pdf_precipitation_1hr_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
