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
from matplotlib.patches import Rectangle
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Specify directories 
domain = 'EURR-3' # CSAM-3


def import_regcm(dirnc, number):

	if len(sys.argv) > 1:
		RCMf = nc(sys.argv[1], mode='r')
	else:
		RCMf = nc(os.path.join(dirnc,domain+'_DOMAIN000.nc'), mode='r')

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


def import_wrf(dirnc, number):

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


if domain == 'CSAM-3':
	dirnc_i = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/icbc'
	dirnc_ii = '/leonardo/home/userexternal/mdasilva/leonardo_work/WRF415/GRID'
	latlon = [-60, 20, -100, -15]
	number_i = 29

	# Import dataset
	lat_i, lon_i, border_mask_i = import_regcm(dirnc_i, number_i)
	lat_ii, lon_ii, border_mask_ii = import_wrf(dirnc_ii, number_i)

else:
	dirnc_i = '/leonardo_work/ICT25_ESP/jdeleeuw/EURR-3/ERA5/high_soil_moisture/ERA5/EURR-3/input'
	latlon = [30, 70, -30, 45]
	number_i = 39

	# Import dataset
	lat_i, lon_i, border_mask_i = import_regcm(dirnc_i, number_i)

# Plot study area
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
font_size = 10

if domain == 'CSAM-3':
	ax.contourf(lon_i, lat_i, border_mask_i, cmap='gray', levels=[0, 1])
	ax.contourf(lon_ii, lat_ii, border_mask_ii, cmap='gray', levels=[0, 1])
	ax.set_extent([latlon[2], latlon[3], latlon[0], latlon[1]], crs=ccrs.PlateCarree())
	ax.stock_img()
	ax.coastlines()
	ax.add_feature(cfeature.BORDERS, linestyle=':')
	plt.text(-24, 11, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
	plt.text(-56, -39.5, u'SESA', color='red', fontsize=font_size, fontweight='bold')
	plt.text(-56, -12, u'RegCM5', color='gray', fontsize=font_size, fontweight='bold')
	plt.text(-56, -55, u'WRF415', color='gray', fontsize=font_size, fontweight='bold')
	ax.add_patch(Rectangle((-60, -35), 10, 13, linewidth=1.5, edgecolor='red', linestyle='--', facecolor='none', transform=ccrs.PlateCarree()))
else:
	ax.contourf(lon_i, lat_i, border_mask_i, cmap='gray', levels=[0, 1])
	ax.set_extent([latlon[2], latlon[3], latlon[0], latlon[1]], crs=ccrs.PlateCarree())
	ax.stock_img()
	ax.coastlines()
	ax.add_feature(cfeature.BORDERS, linestyle=':')
	plt.text(-25, 31.25, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
	plt.text(8, 51, u'APL-3', color='red', fontsize=font_size, fontweight='bold')
	plt.text(1, 66, u'RegCM5', color='gray', fontsize=font_size, fontweight='bold')
	ax.add_patch(Rectangle((1, 40), 16, 10, linewidth=1.5, edgecolor='red', linestyle='--', facecolor='none', transform=ccrs.PlateCarree()))

# Add gridlines 
gridlines = ax.gridlines(draw_labels=True, color='k', linestyle='--', alpha=0.4)
gridlines.right_labels = False
gridlines.bottom_labels = True
gridlines.left_labels = True
gridlines.top_labels = False
gridlines.xlabel_style = {'size': font_size}
gridlines.ylabel_style = {'size': font_size}

# Path out to save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/figs/hpe'
name_out = 'pyplt_maps_study_area_RegCM5_{0}_2000-2009.png'.format(domain)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


