# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot composites"

import os
import netCDF4
import warnings
import argparse
import numpy as np
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeat

from datetime import datetime, timedelta
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

warnings.filterwarnings("ignore", category=RuntimeWarning)

parser = argparse.ArgumentParser(description='Plot composites for given variable')
parser.add_argument('--variable', type=str, required=True, choices=['cape', 'cin'], help='Plot: cape or cin')
args = parser.parse_args()
variable = args.variable

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


def open_file_dt(dataset):

	dt_hr = []
	for yr in range(2018, 2021+1):
	
		data = read_dat_file('{0}/SAM-3km/postproc/cyclone/{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
	
		for i, (header, rows) in enumerate(data):
			if (rows[0][2] < str(-56)):
				dt_hr.append(header[1][:])
				
	return dt_hr
	
	
def import_data(param, dataset, indices):

	if dataset == 'RegCM5':
		arq = '{0}/SAM-3km/postproc/cyclone/RegCM5/{1}_SAM-3km_{2}_6hr_2018-2021_lonlat.nc'.format(path, param, dataset)
	elif dataset == 'WRF415':
		arq = '{0}/SAM-3km/postproc/cyclone/WRF415/{1}_SAM-3km_{2}_6hr_2018-2021_lonlat.nc'.format(path, param, dataset)
	else:
		arq   = '{0}/SAM-3km/postproc/cyclone/ERA5/{1}_SAM-3km_{2}_6hr_2018-2021_lonlat.nc'.format(path, param, dataset)		
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:] 
	lat  = data.variables['lat'][:]	
	lon  = data.variables['lon'][:]
	avg  = var[:][:,:,:]
	
	if dataset == 'RegCM5':
		mean = np.where(avg < 0, np.nan, avg)
	elif dataset == 'WRF415':
		if param == 'CAPE':
			mean = np.where(avg < 0, 0, avg)
		else:
			mean = np.where(avg <= -99999., np.nan, avg)
	else:
		mean = np.where(avg < 0, np.nan, avg)
		
	var_i, var_ii, var_iii = [], [], [] 
	for idx in indices:
		var_i.append(np.squeeze(mean[idx-4,:,:]))
		var_ii.append(np.squeeze(mean[idx,:,:]))
		var_iii.append(np.squeeze(mean[idx+4,:,:]))

	mean_i = np.nanmean(var_i, axis=0)
	mean_ii = np.nanmean(var_ii, axis=0)
	mean_iii = np.nanmean(var_iii, axis=0)
		
	return lat, lon, mean_i, mean_ii, mean_iii


def configure_subplot(ax):

	states_provinces = cfeat.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')

	ax.set_extent([-76, -38.5, -34.5, -15], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(-76,-38.5,5), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(-34.5,-15,5), crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.add_feature(cfeat.BORDERS)
	ax.add_feature(states_provinces, edgecolor='0.25')
	ax.coastlines()	

	
# Generate list of dates from 2018 to 2021
hourly_dates = generate_hourly_dates(datetime(2018, 1, 1, 0), datetime(2021, 12, 31, 23))

# Import cyclone tracking date 
dt_era5 = open_file_dt('ERA5')
dt_regcm5 = open_file_dt('RegCM5')
dt_wrf415 = open_file_dt('WRF415')

era5_idx = find_indices_in_date_list(hourly_dates, dt_era5)
regcm5_idx = find_indices_in_date_list(hourly_dates, dt_regcm5)
wrf415_idx = find_indices_in_date_list(hourly_dates, dt_wrf415)

# Import model and obs dataset 
dict_var = {'cape': ['cape', 'cape', 'AFWA_CAPE_MU'],
'cin': ['cin', 'cin', 'AFWA_CIN_MU']}

lat, lon, var_era5_i, var_era5_ii, var_era5_iii = import_data(dict_var[variable][0], 'ERA5', era5_idx)
lat, lon, var_regcm5_i, var_regcm5_ii, var_regcm5_iii = import_data(dict_var[variable][1], 'RegCM5', regcm5_idx)
lat, lon, var_wrf415_i, var_wrf415_ii, var_wrf415_iii = import_data(dict_var[variable][2], 'WRF415', wrf415_idx)

print(np.nanmin(var_era5_i), np.nanmax(var_era5_i))
print(np.nanmin(var_regcm5_i), np.nanmax(var_regcm5_i))
print(np.nanmin(var_wrf415_i), np.nanmax(var_wrf415_i))

# Plot figure
fig, axes = plt.subplots(3,3, figsize=(14, 9), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9) = axes

states_provinces = cfeat.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')

if variable == 'cape':
	legend = 'CAPE (J Kg$^-$$^1$)'
	level = np.arange(0,1150,50)
	cmap_color = 'Greens'
else:
	legend = 'CIN (J Kg$^-$$^1$)'
	level = np.arange(0,1050,50)
	cmap_color = 'Reds'

ax1.set_title('(a) ERA5 (-24hr)', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
cf1 = ax1.contourf(lon, lat, var_era5_i, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
ct1 = ax1.contour(lon, lat, var_era5_i, levels=10, colors='black', linewidths=0.50)
ax1.clabel(ct1, inline=1, fontsize=8)
configure_subplot(ax1)

cf2 = ax2.contourf(lon, lat, var_regcm5_i, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
ct2 = ax2.contour(lon, lat, var_regcm5_i, levels=10, colors='black', linewidths=0.50)
ax2.clabel(ct2, inline=1, fontsize=8)
ax2.set_title('(b) RegCM5 (-24hr)', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

ax3.set_title('(c) WRF415 (-24hr)', loc='left', fontsize=font_size, fontweight='bold')
cf3 = ax3.contourf(lon, lat, var_wrf415_i, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
ct3 = ax3.contour(lon, lat, var_wrf415_i, levels=10, colors='black', linewidths=0.50)
ax3.clabel(ct3, inline=1, fontsize=8)
configure_subplot(ax3)

cf4 = ax4.contourf(lon, lat, var_era5_ii, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
ct4 = ax4.contour(lon, lat, var_era5_ii, levels=10, colors='black', linewidths=0.50)
ax4.clabel(ct4, inline=1, fontsize=8)
ax4.set_title('(d) ERA5 (cyclogenesis)', loc='left', fontsize=font_size, fontweight='bold')
ax4.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax4)

cf5 = ax5.contourf(lon, lat, var_regcm5_ii, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
ct5 = ax5.contour(lon, lat, var_regcm5_ii, levels=10, colors='black', linewidths=0.50)
ax5.clabel(ct5, inline=1, fontsize=8)
ax5.set_title('(e) RegCM5 (cyclogenesis)', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax5)

cf6 = ax6.contourf(lon, lat, var_wrf415_ii, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
ct6 = ax6.contour(lon, lat, var_wrf415_ii, levels=10, colors='black', linewidths=0.50)
ax6.clabel(ct6, inline=1, fontsize=8)
ax6.set_title('(f) WRF415 (cyclogenesis)', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax6)

cf7 = ax7.contourf(lon, lat, var_era5_iii, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
ct7 = ax7.contour(lon, lat, var_era5_iii, levels=10, colors='black', linewidths=0.50)
ax7.clabel(ct7, inline=1, fontsize=8)
ax7.set_title('(g) ERA5 (+24hr)', loc='left', fontsize=font_size, fontweight='bold')
ax7.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax7.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax7)

cf8 = ax8.contourf(lon, lat, var_regcm5_iii, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap='Greens')
ct8 = ax8.contour(lon, lat, var_regcm5_iii, levels=10, colors='black', linewidths=0.50)
ax8.clabel(ct8, inline=1, fontsize=8)
ax8.set_title('(h) RegCM5 (+24hr)', loc='left', fontsize=font_size, fontweight='bold')
ax8.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax8)

cf9 = ax9.contourf(lon, lat, var_wrf415_iii, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
ct9 = ax9.contour(lon, lat, var_wrf415_iii, levels=10, colors='black', linewidths=0.50)
ax9.clabel(ct9, inline=1, fontsize=8)
ax9.set_title('(i) WRF415 (+24hr)', loc='left', fontsize=font_size, fontweight='bold')
ax9.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax9)

cb = plt.colorbar(cf3, cax=fig.add_axes([0.91, 0.2, 0.015, 0.6]))
cb.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cb.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/SAM-3km/figs/cyclone'.format(path)
name_out = 'pyplt_maps_{0}_CP-RCM_SAM-3km_2018-2021.png'.format(variable)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
