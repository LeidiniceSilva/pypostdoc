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
path = '/leonardo/home/userexternal/mdasilva/leonardo_work'

skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,7793,100,105,106,107,112,117,124,135,137,139,
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

	
def import_data(param, dataset, indices):

	if dataset == 'RegCM5':
		arq = '{0}/SAM-3km/postproc/cyclone/RegCM5/{1}_SAM-3km_{2}_day_2018-2021_lonlat.nc'.format(path, param, dataset)
	elif dataset == 'WRF415':
		arq = '{0}/SAM-3km/postproc/cyclone/WRF415/{1}_SAM-3km_{2}_day_2018-2021_lonlat.nc'.format(path, param, dataset)
	elif dataset == 'CMORPH':
		arq   = '{0}/SAM-3km/postproc/cyclone/CMORPH/{1}_SAM-3km_{2}_day_2018-2021_lonlat.nc'.format(path, param, dataset)		
	else:
		arq   = '{0}/SAM-3km/postproc/cyclone/ERA5/{1}_SAM-3km_{2}_day_2018-2021_lonlat.nc'.format(path, param, dataset)

	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]	
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]

	var_i = []
	for idx_i in indices:
		var_i.append(np.squeeze(mean[idx_i,:,:]))
	
	mean_99 = np.percentile(var_i, 99, axis=0)

	return lat, lon, mean_99	


def import_ws(param, indices):
	
	yy, xx, mean_99 = [], [], []
	for station in range(1, 567):
		print(inmet[station][0])
		
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		yy.append(inmet[station][2])
		xx.append(inmet[station][3])

		arq  = xr.open_dataset('{0}/FPS_SESA/database/obs/inmet/inmet_nc/hourly/{1}/'.format(path, param) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(param, inmet[station][0]))
		data = arq[param]
		time = data.sel(time=slice('2018-01-01','2021-12-31'))
		var  = time.resample(time='1D').sum()	
		var_ = var.values

		var_i = []
		for idx_i in indices:
			var_i.append(var_[idx_i])
				
		mean_99.append(np.nanpercentile(var_i, 99, axis=0))
		
	return yy, xx, mean_99
	

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
lat_, lon_, inmet_idx_99  = import_ws('pre', era5_idx_i)
lat, lon, cmorph_idx_99 = import_data('cmorph', 'CMORPH', era5_idx_i)
lat, lon, era5_idx_99 = import_data('tp', 'ERA5', era5_idx_i)
lat, lon, regcm5_idx_99 = import_data('pr', 'RegCM5', regcm5_idx_i)
lat, lon, wrf415_idx_99 = import_data('PREC_ACC_NC', 'WRF415', wrf415_idx_i)

# Plot figure
fig, axes = plt.subplots(2,3, figsize=(14, 6), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2, ax3), (ax4, ax5, ax6) = axes
fig.delaxes(ax6)

color = ['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']
level = np.arange(0,180,10)

sc1 = ax1.scatter(lon_, lat_, 12, inmet_idx_99, cmap=matplotlib.colors.ListedColormap(color), edgecolors='black', linewidth=0.5, marker='o', vmin=0, vmax=180) 
ax1.set_title('(a) INMET', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

cf2 = ax2.contourf(lon, lat, cmorph_idx_99, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))
ax2.set_title('(b) CMORPH', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

cf3 = ax3.contourf(lon, lat, era5_idx_99, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))
ax3.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax3.set_title('(c) ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

cf4 = ax4.contourf(lon, lat, regcm5_idx_99, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))
ax4.set_title('(d) RegCM5', loc='left', fontsize=font_size, fontweight='bold')
ax4.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax4.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax4)

cf5 = ax5.contourf(lon, lat, wrf415_idx_99, levels=level, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))
ax5.set_title('(e) WRF415', loc='left', fontsize=font_size, fontweight='bold')
ax5.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax5)

cb = plt.colorbar(cf2, cax=fig.add_axes([0.92, 0.3, 0.015, 0.4]))
cb.set_label('Precipitation 99th (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
cb.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/SAM-3km/figs/cyclone'.format(path)
name_out = 'pyplt_maps_precipitation_99th-day_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
exit()
