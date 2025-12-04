# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "May 14, 2025"
__description__ = "This script plot study area"

import os
import sys
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib.pyplot as plt

from cartopy import config
from netCDF4 import Dataset as nc
from matplotlib.patches import Rectangle
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


def import_orog(dirnc, domain):

	number = 29
	
	if len(sys.argv) > 1:
		RCMf = nc(sys.argv[1], mode='r')
	else:
		RCMf = nc(os.path.join(dirnc,'{0}_DOMAIN000.nc'.format(domain)), mode='r')

	lat  = RCMf.variables['xlat'][:,:]
	lon  = RCMf.variables['xlon'][:,:]
	topo = RCMf.variables['topo'][:,:]
	RCMf.close()

	print(topo.shape)
	
	ny,nx = topo.shape
	border_mask = np.full((ny, nx), np.nan)
	border_mask[:number, :] = 1
	border_mask[-number:, :] = 1
	border_mask[:, :number] = 1
	border_mask[:, -number:] = 1

	return lat, lon, border_mask


# Import dataset
lat_i, lon_i, border_mask_i = import_orog('/leonardo/home/userexternal/mdasilva/leonardo_scratch/Otis_exp/domain_large/ctrl/input', 'Otis_exp')
lat_ii, lon_ii, border_mask_ii = import_orog('/leonardo/home/userexternal/mdasilva/leonardo_scratch/Otis_exp/domain_small/ctrl/input', 'Otis_exp')

# Plot figure
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
latlon = [-12, 42, -130, -65]
font_size = 10

ax.contourf(lon_i, lat_i, border_mask_i, colors=['white', 'blue'], levels=[0, 0.5, 1])
ax.contourf(lon_ii, lat_ii, border_mask_ii, colors=['white', 'red'], levels=[0, 0.5, 1])
ax.set_extent([latlon[2], latlon[3], latlon[0], latlon[1]], crs=ccrs.PlateCarree())
ax.set_xticks(np.arange(latlon[2],latlon[3],10), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(latlon[0], latlon[1],8), crs=ccrs.PlateCarree())

for label in ax.get_xticklabels() + ax.get_yticklabels():
	label.set_fontsize(font_size)

ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
ax.grid(c='k', ls='--')
ax.add_feature(cfeat.OCEAN, facecolor='lightblue')
ax.add_feature(cfeat.LAND, facecolor='green')
ax.add_feature(cfeat.BORDERS, linewidth=0.5)
ax.coastlines(linewidth=0.75)

plt.text(-125, -10, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
plt.text(-116, 2, u'd2', color='blue', fontsize=font_size, fontweight='bold')
plt.text(-105, 8, u'd1', color='red', fontsize=font_size, fontweight='bold')

# Path out to save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/figs'
name_out = 'pyplt_Hurricane_Otis_simulated_domains.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

