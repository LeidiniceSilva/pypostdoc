# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 16, 2023"
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
dirnc = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/icbc'
domname = 'CSAM-3'

# RegCM file
if len(sys.argv) > 1:
    RCMf = nc(sys.argv[1], mode='r')
else:
    RCMf = nc(os.path.join(dirnc,domname+'_DOMAIN000.nc'), mode='r')
    
lat  = RCMf.variables['xlat'][:,:]
lon  = RCMf.variables['xlon'][:,:]
topo = RCMf.variables['topo'][:,:]
lonc = RCMf.longitude_of_projection_origin
latc = RCMf.latitude_of_projection_origin
RCMf.close()

# Creating mask of the border
topo = topo[30:-30, 30:-30]
border_mask = np.full((783, 1231), np.nan)
border_mask[:29, :] = 1
border_mask[-29:, :] = 1
border_mask[:, :29] = 1
border_mask[:, -29:] = 1
x_mod = lon[30:-30,30:-30]
y_mod = lat[30:-30,30:-30]

lon_bounds = [-85, -30]
lat_bounds = [-60, 20]

# Plot study area
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
font_size = 10

ax.contourf(lon, lat, border_mask, cmap='gray', levels=[0, 1])
ax.text(-36, 11, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
ax.plot(-47.9000, -16.0000, 'o', ms=4, markeredgewidth=0.75, color='white', mec='black') # Brazilia
ax.plot(-46.6396, -23.5558, 'o', ms=4, markeredgewidth=0.75, color='white', mec='black') # Sao Paulo
ax.plot(-56.1674, -34.8335, 'o', ms=4, markeredgewidth=0.75, color='white', mec='black') # Montevideo
ax.plot(-57.5759, -25.2637, 'o', ms=4, markeredgewidth=0.75, color='white', mec='black') # Asuncion
ax.plot(-58.4004, -34.6051, 'o', ms=4, markeredgewidth=0.75, color='white', mec='black') # Buenos Aires
ax.plot(-68.1193, -16.4897, 'o', ms=4, markeredgewidth=0.75, color='white', mec='black') # La Paz
ax.plot(-70.6693, -33.4489, 'o', ms=4, markeredgewidth=0.75, color='white', mec='black') # Santiago

# Add features
ax.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]], crs=ccrs.PlateCarree())
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
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/figs/evaluate'
name_out = 'pyplt_maps_study_area_RegCM5_CSAM-3_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

