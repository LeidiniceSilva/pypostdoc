# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot map of genesis"

import os
import netCDF4
import argparse
import datetime
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeat

from scipy import signal, misc
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

parser = argparse.ArgumentParser(description="Script to save cyclone genesis")
parser.add_argument("--domain", type=str, help="AUS, CAM, EUR, NAM, SAM, WAS")
args = parser.parse_args()

domain = args.domain

if domain == 'AUS':
    lon1, lon2, lat1, lat2 = 89, 177, -52, -1
elif domain == 'CAM':
    lon1, lon2, lat1, lat2 = 89, 177, -52, -1
elif domain == 'EUR':
    lon1, lon2, lat1, lat2 = 89, 177, -52, -1
elif domain == 'NAM':
    lon1, lon2, lat1, lat2 = 89, 177, -52, -1
elif domain == 'SAM':
    lon1, lon2, lat1, lat2 = -80, -30, -56, -16
else:
    lon1, lon2, lat1, lat2 = 89, 177, -52, -1

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/TRACK-CYCLONE/CORDEX-TF'


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
	
		data = read_dat_file('{0}/{1}/S2R-Vortrack/{2}/track/resultado_{3}.dat'.format(path, dataset, domain, yr))
		for i, (header, rows) in enumerate(data):
			lat.append(rows[0][1])
			lon.append(rows[0][2])

			day, month, year = int(header[1][6:8]), int(header[1][4:6]), int(header[1][0:4])
			date = datetime.datetime(year, month, day)
			dt.append(date)

	return lat, lon, dt
	

def configure_subplot(ax):

        states_provinces = cfeat.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')

        ax.set_extent([lon1, lon2, lat1, lat2], crs=ccrs.PlateCarree())
        ax.set_xticks(np.arange(lon1,lon2,20), crs=ccrs.PlateCarree())
        ax.set_yticks(np.arange(lat1,lat2,10), crs=ccrs.PlateCarree())
        ax.xaxis.set_major_formatter(LongitudeFormatter())
        ax.yaxis.set_major_formatter(LatitudeFormatter())
        ax.grid(c='k', ls='--', alpha=0.4)
        ax.add_feature(cfeat.BORDERS)
        ax.add_feature(states_provinces, edgecolor='0.25')
        ax.coastlines()
	
	
# Import model and obs dataset
lat_obs, lon_obs, dt_obs = open_dat_file('ERA5', 2000, 2009)
lat_gcm1, lon_gcm1, dt_gcm1 = open_dat_file('EC-Earth3-Veg', 2000, 2009)
lat_gcm2, lon_gcm2, dt_gcm2 = open_dat_file('MPI-ESM1-2-HR', 2000, 2009)
lat_gcm3, lon_gcm3, dt_gcm3 = open_dat_file('CNRM-ESM2-1', 2000, 2009)

lat_obs_ = [float(i) for i in lat_obs]
lon_obs_ = [float(i) for i in lon_obs]

lat_gcm1_ = [float(i) for i in lat_gcm1]
lon_gcm1_ = [float(i) for i in lon_gcm1]

lat_gcm2_ = [float(i) for i in lat_gcm2]
lon_gcm2_ = [float(i) for i in lon_gcm2]

lat_gcm3_ = [float(i) for i in lat_gcm3]
lon_gcm3_ = [float(i) for i in lon_gcm3]

# Convert dates to numeric days since a ref date
ref_date = datetime.datetime(2000, 1, 1)
date_obs = [(d - ref_date).days for d in dt_obs]
date_gcm1 = [(d - ref_date).days for d in dt_gcm1]
date_gcm2 = [(d - ref_date).days for d in dt_gcm2]
date_gcm3 = [(d - ref_date).days for d in dt_gcm3]

# Plot figure
fig, axes = plt.subplots(2,2, figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2), (ax3, ax4) = axes

font_size = 10
cmap = plt.cm.viridis_r
norm = matplotlib.colors.Normalize(vmin=0, vmax=(datetime.datetime(2009, 12, 31) - ref_date).days)

sc1 = ax1.scatter(lon_obs_, lat_obs_, c=date_obs, s=25, cmap=cmap, norm=norm, edgecolors='black', linewidth=0.5, marker='o') 
ax1.set_title('(a)', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax1.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

sc2 = ax2.scatter(lon_gcm1_, lat_gcm1_, c=date_gcm1, s=25, cmap=cmap, norm=norm, edgecolors='black', linewidth=0.5, marker='o') 
ax2.set_title('(b)', loc='left', fontsize=font_size, fontweight='bold')
ax2.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

sc = ax3.scatter(lon_gcm2_, lat_gcm2_, c=date_gcm2, s=25, cmap=cmap, norm=norm, edgecolors='black', linewidth=0.5, marker='o') 
ax3.set_title('(c)', loc='left', fontsize=font_size, fontweight='bold')
ax3.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax3.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

sc = ax4.scatter(lon_gcm3_, lat_gcm3_, c=date_gcm3, s=25, cmap=cmap, norm=norm, edgecolors='black', linewidth=0.5, marker='o')
ax4.set_title('(d)', loc='left', fontsize=font_size, fontweight='bold')
ax4.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax4.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax4)

cbar_ax = fig.add_axes([0.92, 0.25, 0.015, 0.5])
cbar = plt.colorbar(sc1, cax=cbar_ax)
cbar.set_label('Days since 2000-01-01')

# Path out to save figure
path_out = '{0}/figs/S2R-Vortrack'.format(path)
name_out = 'pyplt_maps_genesis_{0}_2000-2009.png'.format(domain)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
