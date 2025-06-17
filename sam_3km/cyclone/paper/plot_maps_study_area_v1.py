# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "May 14, 2025"
__description__ = "This script plot study area"

import os
import sys
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

from cartopy import config
from netCDF4 import Dataset as nc
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

domain = 'SA'
latlon = [-60, 20, -100, -15]
number = 29


def import_regcm():

	dirnc = '/leonardo/home/userexternal/mdasilva/leonardo_work/SAM-3km/test/input'

	if len(sys.argv) > 1:
		RCMf = nc(sys.argv[1], mode='r')
	else:
		RCMf = nc(os.path.join(dirnc,'SAM-3km_DOMAIN000.nc'), mode='r')

	lat  = RCMf.variables['xlat'][:,:]
	lon  = RCMf.variables['xlon'][:,:]
	topo = RCMf.variables['topo'][:,:]
	lonc = RCMf.longitude_of_projection_origin
	latc = RCMf.latitude_of_projection_origin
	RCMf.close()

	ny,nx = topo.shape
	border_mask = np.full((ny, nx), np.nan)
	border_mask[:number, :] = 1
	border_mask[-number:, :] = 1
	border_mask[:, :number] = 1
	border_mask[:, -number:] = 1

	return lat, lon, border_mask


def import_wrf():

	dirnc = '/leonardo/home/userexternal/mdasilva/leonardo_work/WRF415/GRID'

	if len(sys.argv) > 1:
		RCMf = nc(sys.argv[1], mode='r')
	else:
		RCMf = nc(os.path.join(dirnc,'wrf_grid_4regrid.nc'), mode='r')

	lat  = RCMf.variables['latitude'][:,:]
	lon  = RCMf.variables['longitude'][:,:]
	topo = RCMf.variables['HGT'][:,:]
	RCMf.close()

	ny,nx = topo.shape
	border_mask = np.full((ny, nx), np.nan)
	border_mask[:number, :] = 1
	border_mask[-number:, :] = 1
	border_mask[:, :number] = 1
	border_mask[:, -number:] = 1

	return lat, lon, border_mask


# Import dataset
lat_i, lon_i, border_mask_i = import_regcm()
lat_ii, lon_ii, border_mask_ii = import_wrf()

# Plot study area
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
font_size = 10

ax.contourf(lon_i, lat_i, border_mask_i, cmap='gray', levels=[0, 1])
ax.contourf(lon_ii, lat_ii, border_mask_ii, cmap='gray', levels=[0, 1])
ax.set_extent([latlon[2], latlon[3], latlon[0], latlon[1]], crs=ccrs.PlateCarree())
ax.stock_img()
ax.coastlines()
ax.add_feature(cfeature.BORDERS, linestyle=':')
plt.text(-24, 11, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
plt.text(-56, -39, u'RegCM5', color='gray', fontsize=font_size, fontweight='bold')
plt.text(-56, -54, u'WRF415', color='gray', fontsize=font_size, fontweight='bold')

# Add gridlines 
gridlines = ax.gridlines(draw_labels=True, color='k', linestyle='--', alpha=0.4)
gridlines.right_labels = False
gridlines.bottom_labels = True
gridlines.left_labels = True
gridlines.top_labels = False
gridlines.xlabel_style = {'size': font_size}
gridlines.ylabel_style = {'size': font_size}

# Path out to save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/SAM-3km/figs/cyclone'
name_out = 'pyplt_maps_study_area_CPMs-{0}_2018-2021.png'.format(domain)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


