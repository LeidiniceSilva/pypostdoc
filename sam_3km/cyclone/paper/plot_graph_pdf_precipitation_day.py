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

font_size = 10
path = '/leonardo/home/userexternal/mdasilva/leonardo_work'

skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,77,93,100,105,106,107,112,117,124,135,137,139,
149,152,155,158,168,174,177,183,186,199,204,210,212,224,226,239,240,248,249,253,254,276,277,280,293,298,
303,305,306,308,319,334,335,341,343,359,362,364,384,393,396,398,399,400,402,413,416,417,422,423,426,427,
443,444,446,451,453,457,458,467,474,479,483,488,489,490,495,505,509,513,514,516,529,534,544,559,566]			


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


def import_data(indices_i, indices_ii, indices_iii):

	mean_, mean_i, mean_ii, mean_iii, mean_iv = [], [], [], [], []
	
	for station in range(1, 567):
		print(station, inmet[station][0])
		
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		arq  = xr.open_dataset('{0}/FPS_SESA/database/obs/inmet/inmet_nc/hourly/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[station][0]))
		data = arq['pre']
		time = data.sel(time=slice('2018-01-01','2021-12-31'))
		var  = time.resample(time='1D').sum()
		var_ = var.values

		arq_i    = xr.open_dataset('{0}/SAM-3km/postproc/cyclone/CMORPH/'.format(path) + 'cmorph_SAM-3km_CMORPH_day_2018-2021_lonlat.nc')
		data_i   = arq_i['cmorph']
		latlon_i = data_i.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		time_i   = latlon_i.sel(time=slice('2018-01-01','2021-12-31'))
		var_i    = time_i.values

		arq_ii    = xr.open_dataset('{0}/SAM-3km/postproc/cyclone/ERA5/'.format(path) + 'tp_SAM-3km_ERA5_day_2018-2021_lonlat.nc')
		data_ii   = arq_ii['tp']
		latlon_ii = data_ii.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		time_ii   = latlon_ii.sel(time=slice('2018-01-01','2021-12-31'))
		var_ii    = time_ii.values

		arq_iii    = xr.open_dataset('{0}/SAM-3km/postproc/cyclone/RegCM5/'.format(path) + 'pr_SAM-3km_RegCM5_day_2018-2021_lonlat.nc')
		data_iii   = arq_iii['pr']
		latlon_iii = data_iii.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		time_iii   = latlon_iii.sel(time=slice('2018-01-01','2021-12-31'))
		var_iii    = time_iii.values

		arq_iv    = xr.open_dataset('{0}/SAM-3km/postproc/cyclone/WRF415/'.format(path) + 'PREC_ACC_NC_SAM-3km_WRF415_day_2018-2021_lonlat.nc')
		data_iv   = arq_iv['PREC_ACC_NC']
		latlon_iv = data_iv.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		time_iv   = latlon_iv.sel(XTIME=slice('2018-01-01','2021-12-31'))
		var_iv    = time_iv.values

		for idx_i in indices_i:
			mean_.append(var_[idx_i])
			mean_i.append(var_i[idx_i])
			mean_ii.append(var_ii[idx_i])

		for idx_ii in indices_ii:
			mean_iii.append(var_iii[idx_ii])

		for idx_iii in indices_iii:
			mean_iv.append(var_iv[idx_iii])
																
	return mean_, mean_i, mean_ii, mean_iii, mean_iv


def comp_pdf(timeseries):

	ts_ = np.array(timeseries)
	ts_list = ts_.flatten()
	ts_round = np.round(ts_list,0)
	ts_filter = ts_round[ts_round > 0.]
	x_pdf_list, pdf_list = np.unique(ts_filter, return_counts=True)
	
	return x_pdf_list, pdf_list
	
		
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
pr_inmet, pr_cmorph, pr_era5, pr_regcm5, pr_wrf415 = import_data(era5_idx_i, regcm5_idx_i, wrf415_idx_i)

# Import pdf function
x_pdf_inmet, pdf_inmet = comp_pdf(pr_inmet)
x_pdf_cmorph, pdf_cmorph = comp_pdf(pr_cmorph)
x_pdf_era5, pdf_era5 = comp_pdf(pr_era5)
x_pdf_regcm5, pdf_regcm5 = comp_pdf(pr_regcm5)
x_pdf_wrf415, pdf_wrf415 = comp_pdf(pr_wrf415)

# Compute 99.9th percentile
p99_inmet = np.nanpercentile(pr_inmet, 99)
p99_cmorph = np.nanpercentile(pr_cmorph, 99)
p99_era5 = np.nanpercentile(pr_era5, 99)
p99_regcm5 = np.nanpercentile(pr_regcm5, 99)
p99_wrf415 = np.nanpercentile(pr_wrf415, 99)

print(p99_inmet)
print(p99_cmorph)
print(p99_era5)
print(p99_regcm5)
print(p99_wrf415)

# Plot figure
fig = plt.figure()
font_size = 8

ax = fig.add_subplot(1, 1, 1)  
plt.plot(x_pdf_inmet,  pdf_inmet,  marker='o', markersize=4, mfc='green',  mec='green',  alpha=0.75, linestyle='None', label='INMET')
plt.plot(x_pdf_cmorph, pdf_cmorph, marker='o', markersize=4, mfc='violet', mec='violet', alpha=0.75, linestyle='None', label='CMORPH')
plt.plot(x_pdf_era5,   pdf_era5,   marker='o', markersize=4, mfc='black',  mec='black',  alpha=0.75, linestyle='None', label='ERA5')
plt.plot(x_pdf_regcm5, pdf_regcm5, marker='o', markersize=4, mfc='blue',   mec='blue',   alpha=0.75, linestyle='None', label='RegCM5')
plt.plot(x_pdf_wrf415, pdf_wrf415, marker='o', markersize=4, mfc='red',    mec='red',    alpha=0.75, linestyle='None', label='WRF415')

plt.axvline(x=p99_inmet, color='green', linestyle='--')
plt.axvline(x=p99_cmorph, color='violet', linestyle='--')
plt.axvline(x=p99_era5, color='black', linestyle='--')
plt.axvline(x=p99_regcm5, color='blue', linestyle='--')
plt.axvline(x=p99_wrf415, color='red', linestyle='--')

plt.title('(a)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('Precipitation (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.grid(axis='y', color='k', linestyle='--', alpha=0.3)
plt.legend(loc=1, ncol=1, fontsize=font_size, shadow=True)

# Path out to save figure
path_out = '{0}/SAM-3km/figs/cyclone'.format(path)
name_out = 'pyplt_graph_pdf_precipitation_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
exit()
