# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot map of genesis density"

import os
import netCDF4
import datetime
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeat

from scipy import signal, misc
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/SAM-3km'


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
			if rows:  # If we have rows, append them to data
				data.append((header, rows))
				rows = []  # Reset rows
			header = line
		else:
			rows.append(line)
	
	# Append the last header and rows to data
	if header and rows:
		data.append((header, rows))
	
	return data


def open_dat_file(dataset, yr_init, yr_end):

	lat, lon, dt = [], [], []
	for yr in range(yr_init, yr_end+1):
	
		data = read_dat_file('{0}/postproc/cyclone/{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
		for i, (header, rows) in enumerate(data):
			lat.append(rows[0][1])
			lon.append(rows[0][2])

			day, month, year = int(header[1][6:8]), int(header[1][4:6]), int(header[1][0:4])
			date = datetime.datetime(year, month, day)
			dt.append(date)

	return lat, lon, dt
	

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
	
	
# Import model and obs dataset
lat_era5, lon_era5, dt_era5 = open_dat_file('ERA5', 2018, 2021)
lat_regcm5, lon_regcm5, dt_regcm5 = open_dat_file('RegCM5', 2018, 2021)
lat_wrf415, lon_wrf415, dt_wrf415 = open_dat_file('WRF415', 2018, 2021)

lat_era5_ = [float(i) for i in lat_era5]
lon_era5_ = [float(i) for i in lon_era5]

lat_regcm5_ = [float(i) for i in lat_regcm5]
lon_regcm5_ = [float(i) for i in lon_regcm5]

lat_wrf415_ = [float(i) for i in lat_wrf415]
lon_wrf415_ = [float(i) for i in lon_wrf415]

# Convert dates to numeric days since a ref date
ref_date = datetime.datetime(2018, 1, 1)
date_era5 = [(d - ref_date).days for d in dt_era5]
date_regcm5 = [(d - ref_date).days for d in dt_regcm5]
date_wrf415 = [(d - ref_date).days for d in dt_wrf415]

# Plot figure
fig, axes = plt.subplots(2,2, figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2), (ax3, ax4) = axes
fig.delaxes(ax4)

font_size = 10
cmap = plt.cm.viridis_r
norm = matplotlib.colors.Normalize(vmin=0, vmax=(datetime.datetime(2021, 12, 31) - ref_date).days)

sc1 = ax1.scatter(lon_era5_, lat_era5_, c=date_era5, s=20, cmap=cmap, norm=norm, edgecolors='black', linewidth=0.5, marker='o') 
ax1.set_title('(a) ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax1.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

sc2 = ax2.scatter(lon_regcm5_, lat_regcm5_, c=date_regcm5, s=20, cmap=cmap, norm=norm, edgecolors='black', linewidth=0.5, marker='o') 
ax2.set_title('(b) RegCM5', loc='left', fontsize=font_size, fontweight='bold')
ax2.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

sc = ax3.scatter(lon_wrf415_, lat_wrf415_, c=date_wrf415, s=20, cmap=cmap, norm=norm, edgecolors='black', linewidth=0.5, marker='o') 
ax3.set_title('(c) WRF415', loc='left', fontsize=font_size, fontweight='bold')
ax3.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax3.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

cbar_ax = fig.add_axes([0.92, 0.25, 0.015, 0.5])
cbar = plt.colorbar(sc1, cax=cbar_ax)
cbar.set_label('Days since 2018-01-01')

# Path out to save figure
path_out = '{0}/figs/cyclone'.format(path)
name_out = 'pyplt_maps_track-v2_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
