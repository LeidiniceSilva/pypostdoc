# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "April 14, 2026"
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

dirnc_i = '/leonardo_work/ICT26_ESP/jdeleeuw/CAR-4/ERA5/high_soil_moisture_OCN/ERA5/CAR-4/input'
dirnc_ii = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/ERA5/icbc'
dirnc_iii = '/leonardo_work/ICT26_ESP/jdeleeuw/EURR-3/ERA5/high_soil_moisture/ERA5/EURR-3/input'
number_i = 29


def import_rcm(dirnc, domain, number):

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


# Import dataset
lat_i, lon_i, border_mask_i = import_rcm(dirnc_i, 'CAR-4', number_i)
lat_ii, lon_ii, border_mask_ii = import_rcm(dirnc_ii, 'CSAM-3', number_i)
lat_iii, lon_iii, border_mask_iii = import_rcm(dirnc_iii, 'EURR-3', number_i)

# Plot study area
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
font_size = 10

ax.contourf(lon_i, lat_i, border_mask_i, cmap='grey', levels=[0, 1])
ax.contourf(lon_ii, lat_ii, border_mask_ii, cmap='grey', levels=[0, 1])
ax.contourf(lon_iii, lat_iii, border_mask_iii, cmap='grey', levels=[0, 1])
ax.stock_img()
ax.coastlines(linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)

plt.text(140, -80, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
plt.text(-160, 10, u'CAR-4', color='grey', fontsize=font_size, fontweight='bold')
plt.text(-56, -50, u'CSAM-3', color='grey', fontsize=font_size, fontweight='bold')
plt.text(-55, 35, u'EURR-3', color='grey', fontsize=font_size, fontweight='bold')

# Add gridlines 
gridlines = ax.gridlines(draw_labels=True, color='k', linestyle='--', alpha=0.5)
gridlines.right_labels = False
gridlines.bottom_labels = True
gridlines.left_labels = True
gridlines.top_labels = False
gridlines.xlabel_style = {'size': font_size}
gridlines.ylabel_style = {'size': font_size}

# Path out to save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/figs'
name_out = 'pyplt_maps_study_area_RegCM5_domains_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

