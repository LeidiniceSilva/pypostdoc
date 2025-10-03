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
parser.add_argument("--domain", type=str, help="AFR, AUS, CAM, EAS, EUR, NAM, SAM, WAS")
args = parser.parse_args()

domain = args.domain

if domain == 'AFR':
    lon1, lon2, lat1, lat2 = -20, 60, -40, -1
elif domain == 'AUS':
    lon1, lon2, lat1, lat2 = 89, 177, -52, -1
elif domain == 'CAM':
    lon1, lon2, lat1, lat2 = -121, -28, 0, 33
elif domain == 'EAS':
    lon1, lon2, lat1, lat2 = 75, 165, -10, 60
elif domain == 'EUR':
    lon1, lon2, lat1, lat2 = -12, 36, 28, 67
elif domain == 'NAM':
    lon1, lon2, lat1, lat2 = -131.5, -53.5, 30.5, 75.5
elif domain == 'SAM':
    lon1, lon2, lat1, lat2 = -85, -35, -56, -16
else:
    lon1, lon2, lat1, lat2 = 40.25, 115.25, -15.75, 45.75

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

        ax.set_extent([lon1, lon2, lat1, lat2], crs=ccrs.PlateCarree())
        ax.set_xticks(np.arange(lon1,lon2,20), crs=ccrs.PlateCarree())
        ax.set_yticks(np.arange(lat1,lat2,10), crs=ccrs.PlateCarree())
        ax.xaxis.set_major_formatter(LongitudeFormatter())
        ax.yaxis.set_major_formatter(LatitudeFormatter())
        ax.grid(c='k', ls='--', alpha=0.4)
        ax.add_feature(cfeat.BORDERS)
        ax.coastlines()
	
	
# Import model and obs dataset
lat_obs, lon_obs, dt_obs = open_dat_file('ERA5', 2000, 2009)
lat_gcm1, lon_gcm1, dt_gcm1 = open_dat_file('EC-Earth3-Veg', 2000, 2009)
lat_gcm2, lon_gcm2, dt_gcm2 = open_dat_file('GFDL-ESM4', 2000, 2009)
lat_gcm3, lon_gcm3, dt_gcm3 = open_dat_file('HadGEM3-GC31-MM', 2000, 2009)
lat_gcm4, lon_gcm4, dt_gcm4 = open_dat_file('MPI-ESM1-2-HR', 2000, 2009)
lat_gcm5, lon_gcm5, dt_gcm5 = open_dat_file('MPI-ESM1-2-LR', 2000, 2009)
lat_gcm6, lon_gcm6, dt_gcm6 = open_dat_file('CNRM-ESM2-1', 2000, 2009)
lat_gcm7, lon_gcm7, dt_gcm7 = open_dat_file('UKESM1-0-LL', 2000, 2009)

lat_obs_ = [float(i) for i in lat_obs]
lon_obs_ = [float(i) for i in lon_obs]

lat_gcm1_ = [float(i) for i in lat_gcm1]
lon_gcm1_ = [float(i) for i in lon_gcm1]

lat_gcm2_ = [float(i) for i in lat_gcm2]
lon_gcm2_ = [float(i) for i in lon_gcm2]

lat_gcm3_ = [float(i) for i in lat_gcm3]
lon_gcm3_ = [float(i) for i in lon_gcm3]

lat_gcm4_ = [float(i) for i in lat_gcm4]
lon_gcm4_ = [float(i) for i in lon_gcm4]

lat_gcm5_ = [float(i) for i in lat_gcm5]
lon_gcm5_ = [float(i) for i in lon_gcm5]

lat_gcm6_ = [float(i) for i in lat_gcm6]
lon_gcm6_ = [float(i) for i in lon_gcm6]

lat_gcm7_ = [float(i) for i in lat_gcm7]
lon_gcm7_ = [float(i) for i in lon_gcm7]

# Convert dates to numeric days since a ref date
ref_date = datetime.datetime(2000, 1, 1)
date_obs = [(d - ref_date).days for d in dt_obs]
date_gcm1 = [(d - ref_date).days for d in dt_gcm1]
date_gcm2 = [(d - ref_date).days for d in dt_gcm2]
date_gcm3 = [(d - ref_date).days for d in dt_gcm3]
date_gcm4 = [(d - ref_date).days for d in dt_gcm4]
date_gcm5 = [(d - ref_date).days for d in dt_gcm5]
date_gcm6 = [(d - ref_date).days for d in dt_gcm6]
date_gcm7 = [(d - ref_date).days for d in dt_gcm7]

# AFR(14, 12); AUS(14, 12); CAM(14, 8); EAS(14, 8); EUR(14, 12); SAM(14, 12); NAM(16, 10); WAS(14, 8)
# Plot figure
fig, axes = plt.subplots(3,3, figsize=(16, 14), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9) = axes
fig.delaxes(ax9)

font_size = 10

sc1 = ax1.scatter(lon_obs_, lat_obs_, s=25, color='gray', edgecolors='black', linewidth=0.5, marker='o') 
ax1.set_title('(a) ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax1.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

sc2 = ax2.scatter(lon_gcm1_, lat_gcm1_, s=25, color='gray', edgecolors='black', linewidth=0.5, marker='o') 
ax2.set_title('(b) EC-Earth3-Veg', loc='left', fontsize=font_size, fontweight='bold')
ax2.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax3.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

sc = ax3.scatter(lon_gcm2_, lat_gcm2_, s=25, color='gray', edgecolors='black', linewidth=0.5, marker='o') 
ax3.set_title('(c) GFDL-ESM4', loc='left', fontsize=font_size, fontweight='bold')
ax3.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax3.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

sc = ax4.scatter(lon_gcm3_, lat_gcm3_, s=25, color='gray', edgecolors='black', linewidth=0.5, marker='o')
ax4.set_title('(d) HadGEM3-GC31-MM', loc='left', fontsize=font_size, fontweight='bold')
ax4.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax4.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax4)

sc = ax5.scatter(lon_gcm4_, lat_gcm4_, s=25, color='gray', edgecolors='black', linewidth=0.5, marker='o')
ax5.set_title('(e) MPI-ESM1-2-HR', loc='left', fontsize=font_size, fontweight='bold')
ax5.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax5.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax5)

sc = ax6.scatter(lon_gcm5_, lat_gcm5_, s=25, color='gray', edgecolors='black', linewidth=0.5, marker='o')
ax6.set_title('(f) MPI-ESM1-2-LR', loc='left', fontsize=font_size, fontweight='bold')
ax6.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax6.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax6)

sc = ax7.scatter(lon_gcm6_, lat_gcm6_, s=25, color='gray', edgecolors='black', linewidth=0.5, marker='o')
ax7.set_title('(g) NorESM-2MM', loc='left', fontsize=font_size, fontweight='bold')
ax7.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax7.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax7)

sc = ax8.scatter(lon_gcm7_, lat_gcm7_, s=25, color='gray', edgecolors='black', linewidth=0.5, marker='o')
ax8.set_title('(h) UKESM1-0-LL', loc='left', fontsize=font_size, fontweight='bold')
ax8.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax8.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax8)

# Path out to save figure
path_out = '{0}/figs/S2R-Vortrack'.format(path)
name_out = 'pyplt_maps_genesis_{0}_2000-2009.png'.format(domain)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
