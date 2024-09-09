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

	if dataset == 'RegCM5':
		arq = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/'.format(path) + '{0}_SAM-3km_{1}_1hr_2018-2021_lonlat.nc'.format(param, dataset))
	elif dataset == 'WRF415':
		arq = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_cyclone/wrf/wrf/{1}/'.format(path, param) + '{0}_SAM-3km_{1}_1hr_2018-2021_lonlat.nc'.format(param, dataset))
	else:
		arq = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + '{0}_SAM-3km_{1}_1hr_2018-2021_lonlat.nc'.format(param, dataset))

	data = arq[param]
	
	if dataset == 'WRF415':
		time = data.isel(XTIME=slice(None,None, 6))
	else:
		time = data.isel(time=slice(None,None, 6))

	lat = time.lat
	lon = time.lon
	var = time.values

	var_i = []
	for idx_i in indices:
		var_i.append(np.squeeze(var[idx_i,:,:]))
	
	mean_95 = np.percentile(var_i, 95, axis=0)
	mean_99 = np.percentile(var_i, 99, axis=0)
		
	return lat, lon, mean_95, mean_99	


def import_ws(param, indices):
	
	yy, xx, mean_95, mean_99 = [], [], [], []
	for station in range(1, 567):
		print(inmet[station][0])
		
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		yy.append(inmet[station][2])
		xx.append(inmet[station][3])

		arq  = xr.open_dataset('{0}/user/mdasilva/WS-SA/INMET/automatic/nc/hourly/{1}/'.format(path, param) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(param, inmet[station][0]))
		data = arq[param]
		time = data.isel(time=slice(None,None, 6))
		var  = time.values

		var_i = []
		for idx_i in indices:
			var_i.append(var[idx_i])
				
		mean_95.append(np.percentile(var_i, 95, axis=0))
		mean_99.append(np.percentile(var_i, 99, axis=0))
		
	return yy, xx, mean_95, mean_99
	
	
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
lat_, lon_, inmet_idx_95, inmet_idx_99 = import_ws('pre', era5_idx_i)
lat, lon, gpm_idx_95, gpm_idx_99 = import_data('precipitation', 'GPM', era5_idx_i)
lat, lon, era5_idx_95, era5_idx_99 = import_data('tp', 'ERA5', era5_idx_i)
lat, lon, regcm5_idx_95, regcm5_idx_99 = import_data('pr', 'RegCM5', regcm5_idx_i)
lat, lon, wrf415_idx_95, wrf415_idx_99 = import_data('PREC_ACC_NC', 'WRF415', wrf415_idx_i)

# Plot figure
fig, axes = plt.subplots(5,2, figsize=(8, 13), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8), (ax9, ax10) = axes

states_provinces = cfeat.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')
color = ['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']
level = np.arange(0,17.5,0.5)

ax1.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax1.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax1.xaxis.set_major_formatter(LongitudeFormatter())
ax1.yaxis.set_major_formatter(LatitudeFormatter())
ax1.grid(c='k', ls='--', alpha=0.3)
ax1.add_feature(cfeat.BORDERS)
ax1.add_feature(states_provinces, edgecolor='0.25')
ax1.coastlines()
ax1.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
ax1.set_title('(a) INMET 95th', loc='left', fontsize=font_size, fontweight='bold')
cf = ax1.contourf(lon, lat, gpm_idx_95-gpm_idx_95, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))
sc = ax1.scatter(lon_, lat_, 12, inmet_idx_95, cmap=matplotlib.colors.ListedColormap(color), edgecolors='black', linewidth=0.5, marker='o', vmin=0, vmax=8) 

ax2.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax2.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax2.xaxis.set_major_formatter(LongitudeFormatter())
ax2.yaxis.set_major_formatter(LatitudeFormatter())
ax2.grid(c='k', ls='--', alpha=0.3)
ax2.add_feature(cfeat.BORDERS)
ax2.add_feature(states_provinces, edgecolor='0.25')
ax2.coastlines()
ax2.set_title('(b) INMET 99th', loc='left', fontsize=font_size, fontweight='bold')
cf = ax2.contourf(lon, lat, gpm_idx_95-gpm_idx_95, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))
sc = ax2.scatter(lon_, lat_, 12, inmet_idx_99, cmap=matplotlib.colors.ListedColormap(color), edgecolors='black', linewidth=0.5, marker='o', vmin=0, vmax=8) 

ax3.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax3.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax3.xaxis.set_major_formatter(LongitudeFormatter())
ax3.yaxis.set_major_formatter(LatitudeFormatter())
ax3.grid(c='k', ls='--', alpha=0.3)
ax3.add_feature(cfeat.BORDERS)
ax3.add_feature(states_provinces, edgecolor='0.25')
ax3.coastlines()
ax3.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
ax3.set_title('(c) GPM 95th', loc='left', fontsize=font_size, fontweight='bold')
cf = ax3.contourf(lon, lat, gpm_idx_95, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))

ax4.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax4.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax4.xaxis.set_major_formatter(LongitudeFormatter())
ax4.yaxis.set_major_formatter(LatitudeFormatter())
ax4.grid(c='k', ls='--', alpha=0.3)
ax4.add_feature(cfeat.BORDERS)
ax4.add_feature(states_provinces, edgecolor='0.25')
ax4.coastlines()
ax4.set_title('(d) GPM 99th', loc='left', fontsize=font_size, fontweight='bold')
cf = ax4.contourf(lon, lat, gpm_idx_99, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))

ax5.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax5.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax5.xaxis.set_major_formatter(LongitudeFormatter())
ax5.yaxis.set_major_formatter(LatitudeFormatter())
ax5.grid(c='k', ls='--', alpha=0.3)
ax5.add_feature(cfeat.BORDERS)
ax5.add_feature(states_provinces, edgecolor='0.25')
ax5.coastlines()
ax5.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
ax5.set_title('(e) ERA5 95th', loc='left', fontsize=font_size, fontweight='bold')
cf = ax5.contourf(lon, lat, era5_idx_95, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))

ax6.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax6.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax6.xaxis.set_major_formatter(LongitudeFormatter())
ax6.yaxis.set_major_formatter(LatitudeFormatter())
ax6.grid(c='k', ls='--', alpha=0.3)
ax6.add_feature(cfeat.BORDERS)
ax6.add_feature(states_provinces, edgecolor='0.25')
ax6.coastlines()
ax6.set_title('(f) ERA5 99th', loc='left', fontsize=font_size, fontweight='bold')
cf = ax6.contourf(lon, lat, era5_idx_99, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))

ax7.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax7.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax7.xaxis.set_major_formatter(LongitudeFormatter())
ax7.yaxis.set_major_formatter(LatitudeFormatter())
ax7.grid(c='k', ls='--', alpha=0.3)
ax7.add_feature(cfeat.BORDERS)
ax7.add_feature(states_provinces, edgecolor='0.25')
ax7.coastlines()
ax7.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
ax7.set_title('(g) RegCM5 95th', loc='left', fontsize=font_size, fontweight='bold')
cf = ax7.contourf(lon, lat, regcm5_idx_95, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))

ax8.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax8.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax8.xaxis.set_major_formatter(LongitudeFormatter())
ax8.yaxis.set_major_formatter(LatitudeFormatter())
ax8.grid(c='k', ls='--', alpha=0.3)
ax8.add_feature(cfeat.BORDERS)
ax8.add_feature(states_provinces, edgecolor='0.25')
ax8.coastlines()
ax8.set_title('(h) RegCM5 99th', loc='left', fontsize=font_size, fontweight='bold')
cf = ax8.contourf(lon, lat, regcm5_idx_99, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))

ax9.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax9.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax9.xaxis.set_major_formatter(LongitudeFormatter())
ax9.yaxis.set_major_formatter(LatitudeFormatter())
ax9.grid(c='k', ls='--', alpha=0.3)
ax9.add_feature(cfeat.BORDERS)
ax9.add_feature(states_provinces, edgecolor='0.25')
ax9.coastlines()
ax9.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
ax9.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax9.set_title('(i) WRF415 95th', loc='left', fontsize=font_size, fontweight='bold')
cf = ax9.contourf(lon, lat, wrf415_idx_95, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))

ax10.set_xticks(np.arange(-76,38.5,7), crs=ccrs.PlateCarree())
ax10.set_yticks(np.arange(-34.5,15,5), crs=ccrs.PlateCarree())
ax10.xaxis.set_major_formatter(LongitudeFormatter())
ax10.yaxis.set_major_formatter(LatitudeFormatter())
ax10.grid(c='k', ls='--', alpha=0.3)
ax10.add_feature(cfeat.BORDERS)
ax10.add_feature(states_provinces, edgecolor='0.25')
ax10.coastlines()
ax10.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax10.set_title('(j) WRF415 99th', loc='left', fontsize=font_size, fontweight='bold')
cf = ax10.contourf(lon, lat, wrf415_idx_99, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))
cb = plt.colorbar(cf, cax=fig.add_axes([0.92, 0.3, 0.015, 0.4]))

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cyclone/paper'.format(path)
name_out = 'pyplt_maps_precipitation_95-99th_1hr_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
