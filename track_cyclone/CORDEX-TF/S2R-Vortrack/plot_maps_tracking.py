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
    lon1, lon2, lat1, lat2 = -85, -35, -56, -16
else:
    lon1, lon2, lat1, lat2 = 89, 177, -52, -1

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/TRACK-CYCLONE/CORDEX-TF'


def read_dat_file(filename):
	data = []
	with open(filename, 'r') as file:
		lines = file.readlines()
		header = []
		rows = []

	for line in lines:
		line = line.strip().split()
		if len(line) == 6:
			if rows:
				data.append((header, rows))
				rows = []
			header = line
		else:
			rows.append(line)
	if header and rows:
		data.append((header, rows))
	return data


def open_tracks(dataset, yr_init, yr_end):
	tracks = []
	ref_date = datetime.datetime(2000, 1, 1)

	for yr in range(yr_init, yr_end + 1):
		data = read_dat_file('{0}/{1}/S2R-Vortrack/{2}/track/resultado_{3}.dat'.format(path, dataset, domain, yr))
		for header, rows in data:
			track = {'lat': [], 'lon': [], 'time': []}
			for row in rows:
				track['lat'].append(float(row[1]))
				track['lon'].append(float(row[2]))
				print(row[0])

				# Format: YYYYMMDDHH
				date_str = row[0]
				year, month, day, hour = int(date_str[0:4]), int(date_str[4:6]), int(date_str[6:8]), int(date_str[8:10])
				dt = datetime.datetime(year, month, day, hour)
				track['time'].append((dt - ref_date).days)
			tracks.append(track)
	return tracks


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
tracks_obs  = open_tracks('ERA5', 2000, 2009)
tracks_gcm1 = open_tracks('EC-Earth3-Veg', 2000, 2009)
tracks_gcm2 = open_tracks('MPI-ESM1-2-HR', 2000, 2009)
tracks_gcm3 = open_tracks('CNRM-ESM2-1', 2000, 2009)

# Plot figure
fig, axes = plt.subplots(2, 2, figsize=(10, 7), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2), (ax3, ax4) = axes

font_size = 10

def plot_tracks(ax, tracks, title):
	for track in tracks:
		ax.plot(track['lon'], track['lat'], color='black', linewidth=0.75, alpha=0.25, transform=ccrs.PlateCarree())
	ax.set_title(title, loc='left', fontsize=font_size, fontweight='bold')
	ax.set_xlabel('Longitude', fontsize=font_size, fontweight='bold')
	ax.set_ylabel('Latitude', fontsize=font_size, fontweight='bold')
	configure_subplot(ax)

plot_tracks(ax1, tracks_obs, '(a)')
plot_tracks(ax2, tracks_gcm1, '(b)')
plot_tracks(ax3, tracks_gcm2, '(c)')
plot_tracks(ax4, tracks_gcm3, '(d)')

# Path out to save figure
path_out = '{0}/figs/S2R-Vortrack'.format(path)
name_out = 'pyplt_maps_tracking_{0}_2000-2009.png'.format(domain)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
