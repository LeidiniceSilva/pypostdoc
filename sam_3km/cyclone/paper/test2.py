# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot composites"

import os
import netCDF4
import warnings
import numpy as np
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeat

from datetime import datetime, timedelta
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

warnings.filterwarnings("ignore", category=RuntimeWarning)

path='/marconi/home/userexternal/mdasilva'


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
	for yr in range(2018, 2018+1):
	
		data = read_dat_file('{0}/user/mdasilva/SAM-3km/post_cyclone/ECv2/{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
	
		for i, (header, rows) in enumerate(data):
			if (rows[0][2] < str(-56)):
				dt_hr.append(header[1][:])
				
	return dt_hr
	
	
def import_data(param, dataset, indices):

	if dataset == 'RegCM5':
		arq = '{0}/user/mdasilva/SAM-3km/post_cyclone/regcm5/regcm5/{1}/{1}_SAM-3km_{2}_6hr_2018_lonlat.nc'.format(path, param, dataset)
	elif dataset == 'WRF415':
		arq = '{0}/user/mdasilva/SAM-3km/post_cyclone/wrf/wrf/{1}/{1}_SAM-3km_{2}_6hr_2018_lonlat.nc'.format(path, param, dataset)
	else:
		arq   = '{0}/user/mdasilva/SAM-3km/post_cyclone/era5/era5/{1}_SAM-25km_{2}_6hr_2018_lonlat.nc'.format(path, param, dataset)		
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:] 
	lat  = data.variables['lat'][:]	
	lon  = data.variables['lon'][:]
	avg  = var[:][:,:,:]
	
	if dataset == 'RegCM5':
		mean = np.where(avg < 0, np.nan, avg)
	elif dataset == 'WRF415':
		if param == 'AFWA_CAPE_MU':
			mean = np.where(avg < 0, np.nan, avg)
		else:
			mean_ = np.where(avg <= -99999., np.nan, avg)
			mean = mean_ * -1.0
	else:
		mean = np.where(avg <= -32767.0, np.nan, avg)

	var_i, var_ii, var_iii = [], [], [] 
	for idx in indices:
		var_i.append(np.squeeze(mean[idx-4,:,:]))
		var_ii.append(np.squeeze(mean[idx,:,:]))
		var_iii.append(np.squeeze(mean[idx+4,:,:]))

	mean_i = np.nanmean(var_i, axis=0)
	mean_ii = np.nanmean(var_ii, axis=0)
	mean_iii = np.nanmean(var_iii, axis=0)
		
	return lat, lon, mean_i, mean_ii, mean_iii
	
	
# Generate list of dates from 2018 to 2021
hourly_dates = generate_hourly_dates(datetime(2018, 1, 1, 0), datetime(2018, 12, 31, 23))

# Import cyclone tracking date 
dt_era5 = open_file_dt('ERA5')
dt_regcm5 = open_file_dt('RegCM5')
dt_wrf415 = open_file_dt('WRF415')

era5_idx = find_indices_in_date_list(hourly_dates, dt_era5)
regcm5_idx = find_indices_in_date_list(hourly_dates, dt_regcm5)
wrf415_idx = find_indices_in_date_list(hourly_dates, dt_wrf415)

# Import model and obs dataset 
lat, lon, cape_era5_i, cape_era5_ii, cape_era5_iii = import_data('AFWA_CAPE_MU', 'WRF415', era5_idx)
lat, lon, cin_era5_i, cin_era5_ii, cin_era5_iii = import_data('AFWA_CIN_MU', 'WRF415', era5_idx)

print(np.nanmin(cape_era5_i), np.nanmax(cape_era5_i))
print(np.nanmin(cin_era5_i), np.nanmax(cin_era5_i))

# Plot figure
fig, axes = plt.subplots(3,3, figsize=(14, 9), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9) = axes
font_size = 10

states_provinces = cfeat.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')
level1 = np.arange(0,1020,20)
level2 = np.arange(0,1000,200)

ax1.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax1.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax1.xaxis.set_major_formatter(LongitudeFormatter())
ax1.yaxis.set_major_formatter(LatitudeFormatter())
ax1.grid(c='k', ls='--', alpha=0.3)
ax1.add_feature(cfeat.BORDERS)
ax1.add_feature(states_provinces, edgecolor='0.25')
ax1.coastlines()
ax1.set_title('(a) ERA5 (-24hr)', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
cf1 = ax1.contourf(lon, lat, cape_era5_i, levels=level1, transform=ccrs.PlateCarree(), extend='max', cmap='Blues')
ct1 = ax1.contour(lon, lat, cin_era5_i, levels=5, colors='red', linewidths=0.50)
ax1.clabel(ct1, inline=1, fontsize=8)

ax2.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax2.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax2.xaxis.set_major_formatter(LongitudeFormatter())
ax2.yaxis.set_major_formatter(LatitudeFormatter())
ax2.grid(c='k', ls='--', alpha=0.3)
ax2.add_feature(cfeat.BORDERS)
ax2.add_feature(states_provinces, edgecolor='0.25')
ax2.coastlines()
ax2.set_title('(b) RegCM5 (-24hr)', loc='left', fontsize=font_size, fontweight='bold')
cf2 = ax2.contourf(lon, lat, cape_era5_i, levels=level1, transform=ccrs.PlateCarree(), extend='max', cmap='Blues')
ct2 = ax2.contour(lon, lat, cin_era5_i, levels=level2, colors='red', linewidths=0.50)
ax2.clabel(ct2, inline=1, fontsize=8)

ax3.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax3.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax3.xaxis.set_major_formatter(LongitudeFormatter())
ax3.yaxis.set_major_formatter(LatitudeFormatter())
ax3.grid(c='k', ls='--', alpha=0.3)
ax3.add_feature(cfeat.BORDERS)
ax3.add_feature(states_provinces, edgecolor='0.25')
ax3.coastlines()
ax3.set_title('(c) WRF415 (-24hr)', loc='left', fontsize=font_size, fontweight='bold')
cf3 = ax3.contourf(lon, lat, cape_era5_i, levels=level1, transform=ccrs.PlateCarree(), extend='max', cmap='Blues')
ct3 = ax3.contour(lon, lat, cin_era5_i, levels=level2, colors='black', linewidths=0.50)
ax3.clabel(ct3, inline=1, fontsize=8)
cb = plt.colorbar(cf3, cax=fig.add_axes([0.91, 0.2, 0.015, 0.6]))

ax4.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax4.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax4.xaxis.set_major_formatter(LongitudeFormatter())
ax4.yaxis.set_major_formatter(LatitudeFormatter())
ax4.grid(c='k', ls='--', alpha=0.3)
ax4.add_feature(cfeat.BORDERS)
ax4.add_feature(states_provinces, edgecolor='0.25')
ax4.coastlines()
ax4.set_title('(d) ERA5 (cyclogenesis)', loc='left', fontsize=font_size, fontweight='bold')
ax4.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
cf4 = ax4.contourf(lon, lat, cape_era5_ii, levels=level1, transform=ccrs.PlateCarree(), extend='max', cmap='Blues')
ct4 = ax4.contour(lon, lat, cin_era5_ii, levels=level2, colors='black', linewidths=0.50)
ax4.clabel(ct4, inline=1, fontsize=8)

ax5.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax5.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax5.xaxis.set_major_formatter(LongitudeFormatter())
ax5.yaxis.set_major_formatter(LatitudeFormatter())
ax5.grid(c='k', ls='--', alpha=0.3)
ax5.add_feature(cfeat.BORDERS)
ax5.add_feature(states_provinces, edgecolor='0.25')
ax5.coastlines()
ax5.set_title('(e) RegCM5 (cyclogenesis)', loc='left', fontsize=font_size, fontweight='bold')
cf5 = ax5.contourf(lon, lat, cape_era5_ii, levels=level1, transform=ccrs.PlateCarree(), extend='max', cmap='Blues')
ct5 = ax5.contour(lon, lat, cin_era5_ii, levels=level2, colors='black', linewidths=0.50)
ax5.clabel(ct5, inline=1, fontsize=8)

ax6.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax6.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax6.xaxis.set_major_formatter(LongitudeFormatter())
ax6.yaxis.set_major_formatter(LatitudeFormatter())
ax6.grid(c='k', ls='--', alpha=0.3)
ax6.add_feature(cfeat.BORDERS)
ax6.add_feature(states_provinces, edgecolor='0.25')
ax6.coastlines()
ax6.set_title('(f) WRF415 (cyclogenesis)', loc='left', fontsize=font_size, fontweight='bold')
cf6 = ax6.contourf(lon, lat, cape_era5_ii, levels=level1, transform=ccrs.PlateCarree(), extend='max', cmap='Blues')
ct6 = ax6.contour(lon, lat, cin_era5_ii, levels=level2, colors='black', linewidths=0.50)
ax6.clabel(ct6, inline=1, fontsize=8)

ax7.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax7.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax7.xaxis.set_major_formatter(LongitudeFormatter())
ax7.yaxis.set_major_formatter(LatitudeFormatter())
ax7.grid(c='k', ls='--', alpha=0.3)
ax7.add_feature(cfeat.BORDERS)
ax7.add_feature(states_provinces, edgecolor='0.25')
ax7.coastlines()
ax7.set_title('(g) ERA5 (+24hr)', loc='left', fontsize=font_size, fontweight='bold')
ax7.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax7.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
cf7 = ax7.contourf(lon, lat, cape_era5_iii, levels=level1, transform=ccrs.PlateCarree(), extend='max', cmap='Blues')
ct7 = ax7.contour(lon, lat, cin_era5_iii, levels=level2, colors='black', linewidths=0.50)
ax7.clabel(ct7, inline=1, fontsize=8)

ax8.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax8.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax8.xaxis.set_major_formatter(LongitudeFormatter())
ax8.yaxis.set_major_formatter(LatitudeFormatter())
ax8.grid(c='k', ls='--', alpha=0.3)
ax8.add_feature(cfeat.BORDERS)
ax8.add_feature(states_provinces, edgecolor='0.25')
ax8.coastlines()
ax8.set_title('(h) RegCM5 (+24hr)', loc='left', fontsize=font_size, fontweight='bold')
ax8.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
cf8 = ax8.contourf(lon, lat, cape_era5_iii, levels=level1, transform=ccrs.PlateCarree(), extend='max', cmap='Blues')
ct8 = ax8.contour(lon, lat, cin_era5_iii, levels=level2, colors='black', linewidths=0.50)
ax8.clabel(ct8, inline=1, fontsize=8)

ax9.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax9.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax9.xaxis.set_major_formatter(LongitudeFormatter())
ax9.yaxis.set_major_formatter(LatitudeFormatter())
ax9.grid(c='k', ls='--', alpha=0.3)
ax9.add_feature(cfeat.BORDERS)
ax9.add_feature(states_provinces, edgecolor='0.25')
ax9.coastlines()
ax9.set_title('(i) WRF415 (+24hr)', loc='left', fontsize=font_size, fontweight='bold')
ax9.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
cf9 = ax9.contourf(lon, lat, cape_era5_iii, levels=level1, transform=ccrs.PlateCarree(), extend='max', cmap='Blues')
ct9 = ax9.contour(lon, lat, cin_era5_iii, levels=level2, colors='black', linewidths=0.50)
ax9.clabel(ct9, inline=1, fontsize=8)
cb = plt.colorbar(cf9, cax=fig.add_axes([0.91, 0.2, 0.015, 0.6]))

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cyclone/paper'.format(path)
name_out = 'pyplt_maps_cape_cin_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
