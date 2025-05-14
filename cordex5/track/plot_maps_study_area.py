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

# Specify directories 
#domain = 'CSAM-3'
domain = 'EURR-3'

if domain == 'CSAM-3':
	dirnc = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/ERA5/icbc'
	latlon = [-60, 20, -85, -30]
	number = 29
else:
	dirnc = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/ERA5-EURR-3/input'
	latlon = [15, 75, -45, 65]
	number = 39

# RegCM file
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

# Creating mask of the border
border_mask = np.full((*lat[0].shape, *lon[0].shape), np.nan)
border_mask[:number, :] = 1
border_mask[-number:, :] = 1
border_mask[:, :number] = 1
border_mask[:, -number:] = 1

# Plot study area
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
font_size = 10

ax.contourf(lon, lat, border_mask, cmap='gray', levels=[0, 1])
ax.set_extent([latlon[2], latlon[3], latlon[0], latlon[1]], crs=ccrs.PlateCarree())
ax.stock_img()
ax.coastlines()
ax.add_feature(cfeature.BORDERS, linestyle=':')

# Add gridlines with labels on right and bottom
gridlines = ax.gridlines(draw_labels=True)
gridlines.right_labels = False
gridlines.bottom_labels = True
gridlines.left_labels = True
gridlines.top_labels = False
gridlines.xlabel_style = {'size': 10}
gridlines.ylabel_style = {'size': 10}

# Path out to save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/figs/track'
name_out = 'pyplt_maps_study_area_RegCM5_{0}_2000-2009.png'.format(domain)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

