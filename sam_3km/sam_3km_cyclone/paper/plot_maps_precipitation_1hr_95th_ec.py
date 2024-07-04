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
import pandas as pd
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap

var = 'pr'
font_size = 10
path='/marconi/home/userexternal/mdasilva'

skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,7793,100,105,106,107,112,117,124,135,137,139,
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
	
	yy, xx, mean = [], [], [] 
	for station in range(1, 4): # 567	
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue
		print(station)

		yy.append(inmet[station][2])
		xx.append(inmet[station][3])

		arq  = xr.open_dataset('{0}/user/mdasilva/WS-SA/INMET/nc/hourly/{1}/'.format(path, param) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(param, inmet[station][0]))
		data = arq[param]
		time = data.sel(time=slice('2018-01-01','2021-12-31'))
		var  = time.resample(time='6H').sum()
		var_ = var.values
		
		var_i = []
		for idx_i in indices:
			var_i.append(var_[idx_i])
				
		mean.append(np.percentile(var_i, 95, axis=0))
		
	return yy, xx, mean


def import_obs(param, dataset, indices):

	arq = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + '{0}_SAM-3km_{1}_1hr_2018-2021_lonlat.nc'.format(param, dataset))
	data = arq[param]
	time = data.sel(time=slice('2018-01-01','2021-12-31'))
	var = time.resample(time='6H').sum()
	lat = var.lat
	lon = var.lon
	var_ = var.values

	var_i = []
	for idx_i in indices:
		var_i.append(np.squeeze(var_[idx_i,:,:]))
	
	mean = np.percentile(var_i, 95, axis=0)

	
	return lat, lon, mean


def import_rcm(param, dataset, indices):

	arq = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/'.format(path) + '{0}_SAM-3km_{1}_1hr_2018-2021_lonlat.nc'.format(param, dataset))
	data = arq[param]
	time = data.sel(time=slice('2018-01-01','2021-12-31'))
	var = time.resample(time='6H').sum()
	lat = var.lat
	lon = var.lon
	var_ = var.values

	var_i = []
	for idx_i in indices:
		var_i.append(np.squeeze(var_[idx_i,:,:]))
	
	mean = np.percentile(var_i, 95, axis=0)

	
	return lat, lon, mean
	
	
def basemap(lat, lon):
	
	map = Basemap(projection='cyl', llcrnrlon=-80., llcrnrlat=-38., urcrnrlon=-34.,urcrnrlat=-8., resolution='c')
	map.drawmeridians(np.arange(-80., -34., 12.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-38., -8., 6.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black')
	
	lons, lats = np.meshgrid(lon, lat)
	xx, yy = map(lons, lats)
	
	return map, xx, yy
		
	
# List of hourly dates 
hourly_dates = generate_hourly_dates(datetime(2018, 1, 1, 0), datetime(2021, 12, 31, 23))

# Import cyclone tracking date 
dt_era5 = open_dat_file('ERA5')
dt_regcm5 = open_dat_file('RegCM5')

# List of indices
idx_era5 = find_indices_in_date_list(hourly_dates, dt_era5)
idx_regcm = find_indices_in_date_list(hourly_dates, dt_regcm5)

# Import model and obs dataset 
lat_i, lon_i, pr_inmet = import_ws('pre', idx_era5)
lat_ii, lon_ii, pr_gpm = import_obs('precipitation', 'GPM', idx_era5)
lat_iii, lon_iii, pr_era5 = import_obs('tp', 'ERA5', idx_era5)
lat_iv, lon_iv, pr_regcm5 = import_rcm('pr', 'RegCM5', idx_era5)

# Plot figure
fig = plt.figure(figsize=(14, 8))

color = ['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']
dict_plot = {'pr': ['Precipitation (mm h$^-$$^1$)', np.arange(0, 20, 1), matplotlib.colors.ListedColormap(color)]}

ax = fig.add_subplot(2, 3, 1)  
map, xx, yy = basemap(lat_ii, lon_ii)
map.scatter(lon_i, lat_i, 12, pr_inmet, cmap=dict_plot[var][2], vmin=0, vmax=12, edgecolors='black', linewidth=0.5, marker='o') 
plt.title(u'(a) INMET', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 3, 2)  
map, xx, yy = basemap(lat_ii, lon_ii)
plt_map = map.contourf(xx, yy, pr_gpm, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(b) GPM', loc='left', fontsize=font_size)

ax = fig.add_subplot(2, 3, 3)  
map, xx, yy = basemap(lat_iii, lon_iii)
plt_map = map.contourf(xx, yy, pr_era5, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(c) ERA5', loc='left', fontsize=font_size)

ax = fig.add_subplot(2, 3, 4)  
map, xx, yy = basemap(lat_iv, lon_iv)
plt_map = map.contourf(xx, yy, pr_regcm5, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(d) RegCM5', loc='left', fontsize=font_size)

ax = fig.add_subplot(2, 3, 5)  
map, xx, yy = basemap(lat_ii, lon_ii)
plt_map = map.contourf(xx, yy, pr_gpm, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(e) WRF', loc='left', fontsize=font_size)

cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.94, 0.3, 0.015, 0.4]))
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cyclone/paper'.format(path)
name_out = 'pyplt_maps_95th_hourly_precipitation_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
