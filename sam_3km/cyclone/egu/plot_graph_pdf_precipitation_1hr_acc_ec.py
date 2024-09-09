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
from dict_inmet_stations import inmet
from datetime import datetime, timedelta
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

path='/marconi/home/userexternal/mdasilva'

skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,77, 93,100,105,106,107,112,117,124,135,137,139,
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
	
		data = read_dat_file('{0}/user/mdasilva/SAM-3km/post_cyclone/ECyclone/ECyclone_{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
	
		rows_list = []
		rows_list_i = []
		for i, (header, rows) in enumerate(data):
			rows_list.append(rows)
		
		for j  in rows_list:
			for k in j:
				dt.append(str(k[0][:]))
				
	return dt


def import_ws(param, indices):

	mean = []
	for station in range(1, 567):
		print(station, inmet[station][0])
		
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		arq  = xr.open_dataset('{0}/user/mdasilva/WS-SA/INMET/nc/hourly/{1}/'.format(path, param) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[station][0]))
		data = arq[param]
		time = data.sel(time=slice('2018-01-01','2021-12-31'))
		var  = time.values
		var_ = var[::6]

		for idx_i in indices:
			mean.append(var_[idx_i])
																
	return mean


def import_sat(param, indices):

	mean = []
	for station in range(1, 567):
		print(station, inmet[station][0])
		
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		arq    = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/'.format(path) + 'precipitation_SAM-3km_GPM_1hr_2018-2021_lonlat.nc')
		data   = arq[param]
		latlon = data.sel(lat=slice(inmet[station][2]-0.1,inmet[station][2]+0.1),lon=slice(inmet[station][3]-0.1,inmet[station][3]+0.1)).mean(('lat','lon'))
		time   = latlon.sel(time=slice('2018-01-01','2021-12-31'))
		var    = time.values
		var_   = var[::6]
		
		for idx_i in indices:
			mean.append(var_[idx_i])
																
	return mean	
	
	
def import_obs(param, indices):

	mean = []
	for station in range(1, 567):
		print(station, inmet[station][0])
		
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		arq    = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/'.format(path) + 'tp_SAM-3km_ERA5_1hr_2018-2021_lonlat.nc')
		data   = arq[param]
		latlon = data.sel(lat=slice(inmet[station][2]-0.25,inmet[station][2]+0.25),lon=slice(inmet[station][3]-0.25,inmet[station][3]+0.25)).mean(('lat','lon'))
		time   = latlon.sel(time=slice('2018-01-01','2021-12-31'))
		var    = time.values
		var_   = var[::6]

		for idx_i in indices:
			mean.append(var_[idx_i])
																
	return mean		


def import_rcm(param, indices):

	mean = []
	for station in range(1, 567):
		print(station, inmet[station][0])
		
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		arq    = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/'.format(path) + 'pr_SAM-3km_RegCM5_1hr_2018-2021_lonlat.nc')
		data   = arq[param]
		latlon = data.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		time   = latlon.sel(time=slice('2018-01-01','2021-12-31'))
		var    = time.values
		var_   = var[::6]

		for idx_i in indices:
			mean.append(var_[idx_i])
													
	return mean	
	
	

# List of hourly dates 
hourly_dates = generate_hourly_dates(datetime(2018, 1, 1, 0), datetime(2021, 12, 31, 23))

# Import cyclone tracking date 
dt_era5 = open_dat_file('ERA5')
dt_regcm5 = open_dat_file('RegCM5')

# List of indices
idx_era5 = find_indices_in_date_list(hourly_dates, dt_era5)
idx_regcm = find_indices_in_date_list(hourly_dates, dt_regcm5)

# Import model and obs dataset 
pr_inmet = import_ws('pre', idx_era5)
pr_gpm = import_sat('precipitation', idx_era5)
pr_era5 = import_obs('tp', idx_era5)
pr_regcm5 = import_rcm('pr', idx_regcm)

# Compute pdf 
pr_inmet_xr = np.array(pr_inmet)
pr_inmet_list = pr_inmet_xr.flatten()
pr_inmet_round = np.round(pr_inmet_list,0)
pr_inmet_filter = pr_inmet_round[pr_inmet_round > 0.]
x_pdf_inmet, pdf_inmet = np.unique(pr_inmet_filter, return_counts=True)

pr_era5_xr = np.array(pr_era5)
pr_era5_list = pr_era5_xr.flatten()
pr_era5_round = np.round(pr_era5_list,0)
pr_era5_filter = pr_era5_round[pr_era5_round > 0.]
x_pdf_era5, pdf_era5 = np.unique(pr_era5_filter, return_counts=True)

pr_gpm_xr = np.array(pr_gpm)
pr_gpm_list = pr_gpm_xr.flatten()
pr_gpm_round = np.round(pr_gpm_list,0)
pr_gpm_filter = pr_gpm_round[pr_gpm_round > 0.]
x_pdf_gpm, pdf_gpm = np.unique(pr_gpm_filter, return_counts=True)

pr_regcm5_xr = np.array(pr_regcm5)
pr_regcm5_list = pr_regcm5_xr.flatten()
pr_regcm5_round = np.round(pr_regcm5_list,0)
pr_regcm5_filter = pr_regcm5_round[pr_regcm5_round > 0.]
x_pdf_regcm5, pdf_regcm5 = np.unique(pr_regcm5_filter, return_counts=True)

# Plot figure
fig = plt.figure(figsize=(6, 9))
font_size = 8

ax = fig.add_subplot(3, 1, 1)  
plt.plot(x_pdf_inmet,  pdf_inmet,  marker='o', markersize=4, mfc='black', mec='black', alpha=0.75, linestyle='None', label='INMET')
plt.plot(x_pdf_era5,   pdf_era5,   marker='o', markersize=4, mfc='red',   mec='red',   alpha=0.75, linestyle='None', label='ERA5')
plt.plot(x_pdf_gpm,    pdf_gpm,    marker='o', markersize=4, mfc='green', mec='green', alpha=0.75, linestyle='None', label='GPM')
plt.plot(x_pdf_regcm5, pdf_regcm5, marker='o', markersize=4, mfc='blue',  mec='blue',  alpha=0.75, linestyle='None', label='RegCM5')
plt.xlabel('Precipitation (mm h$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.legend(loc=1, ncol=2, fontsize=font_size, shadow=True)

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cyclone/egu'.format(path)
name_out = 'pyplt_graph_pdf_1hr_precipitation_EC_ERA5_RegCM5_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
exit()
