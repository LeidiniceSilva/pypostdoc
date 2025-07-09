# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plots cyclone tracks and genesis points on a map"

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
	ref_date = datetime.datetime(2018, 1, 1)

	for yr in range(yr_init, yr_end + 1):
		data = read_dat_file(f'{path}/postproc/cyclone/{dataset}/track/resultado_{yr}.dat')
		for header, rows in data:
			track = {'lat': [], 'lon': [], 'time': []}
			for row in rows:
				track['lat'].append(float(row[1]))
				track['lon'].append(float(row[2]))

				# Format: YYYYMMDDHH
				date_str = row[0]
				year, month, day, hour = int(date_str[0:4]), int(date_str[4:6]), int(date_str[6:8]), int(date_str[8:10])
				dt = datetime.datetime(year, month, day, hour)
				track['time'].append((dt - ref_date).days)
			tracks.append(track)
	return tracks


def configure_subplot(ax):
	states_provinces = cfeat.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')
	ax.set_extent([-76, -38.5, -34.5, -15], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(-76, -38.5, 5), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(-34.5, -15, 5), crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.add_feature(cfeat.BORDERS)
	ax.add_feature(states_provinces, edgecolor='0.25')
	ax.coastlines()


# Read cyclone tracks
tracks_era5   = open_tracks('ERA5', 2018, 2021)
tracks_regcm5 = open_tracks('RegCM5', 2018, 2021)
tracks_wrf415 = open_tracks('WRF415', 2018, 2021)

# Plot figure
fig, axes = plt.subplots(2, 2, figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2), (ax3, ax4) = axes
fig.delaxes(ax4)

font_size = 10
cmap = plt.cm.viridis_r
ref_date = datetime.datetime(2018, 1, 1)
norm = matplotlib.colors.Normalize(vmin=0, vmax=(datetime.datetime(2021, 12, 31) - ref_date).days)

def plot_tracks(ax, tracks, title):
	for track in tracks:
		ax.plot(track['lon'], track['lat'], marker='o', markersize=4, color='black',  markerfacecolor='gray', linewidth=0.75, alpha=0.5, transform=ccrs.PlateCarree())
	ax.set_title(title, loc='left', fontsize=font_size, fontweight='bold')
	ax.set_xlabel('Longitude', fontsize=font_size, fontweight='bold')
	ax.set_ylabel('Latitude', fontsize=font_size, fontweight='bold')
	configure_subplot(ax)

plot_tracks(ax1, tracks_era5, '(a) ERA5')
plot_tracks(ax2, tracks_regcm5, '(b) RegCM5')
plot_tracks(ax3, tracks_wrf415, '(c) WRF415')

# Save figure
path_out = f'{path}/figs/cyclone'
name_out = 'pyplt_maps_track_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

