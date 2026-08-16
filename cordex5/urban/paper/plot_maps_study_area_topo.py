# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 28, 2026"
__description__ = "This script plot study area"

import os
import sys
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

from cartopy import config
from netCDF4 import Dataset as nc
from matplotlib.patches import Polygon
from matplotlib.patches import PathPatch
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Specify directories 
dirnc = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/ERA5/icbc'
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
border_mask = np.full((783, 1231), np.nan)
border_mask[:29, :] = 1
border_mask[-29:, :] = 1
border_mask[:, :29] = 1
border_mask[:, -29:] = 1
x_mod = lon[30:-30,30:-30]
y_mod = lat[30:-30,30:-30]

lon_bounds = [np.min(lon), np.max(lon)]
lat_bounds = [np.min(lat), np.max(lat)]

# Plot study area
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
font_size = 10

ct=ax.contourf(lon, lat, topo, np.arange(0, 3100, 100), cmap='terrain', extend='max')
ax.contourf(lon, lat, border_mask, cmap='gray', levels=[0, 1])
ax.text(-42, -33, u'\u25B2 \nN', color='white', fontsize=font_size, fontweight='bold')
ax.plot(-46.6396, -23.5558, 'o', ms=5, markeredgewidth=0.75, color='white', mec='black') # Sao Paulo
ax.plot(-58.4004, -34.6051, 'o', ms=5, markeredgewidth=0.75, color='white', mec='black') # Buenos Aires
ax.set_xlabel(u'Longitude', fontweight='bold')
ax.set_ylabel(u'Latitude', fontweight='bold')
ax.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]], crs=ccrs.PlateCarree())
ax.set_xticks(np.arange(lon_bounds[0], lon_bounds[1], 10), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(lat_bounds[0], lat_bounds[1], 5), crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
ax.grid(c='gray', ls='--', alpha=0.5)  

# Add features
ax.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linewidth=0.75)
ax.coastlines(linewidth=0.75)

c1,d1 = (-63,-35.25)
c2,d2 = (-63,-14.4)
c3,d3 = (-39.5,-14.4)
c4,d4 = (-39.5,-35.25)
poly2 = Polygon([(c1,d1),(c2,d2),(c3,d3),(c4,d4)], facecolor='none', edgecolor='red', linewidth=0.75)
plt.gca().add_patch(poly2)

# Bar topo
cbar = plt.colorbar(ct, ax=ax, orientation='vertical', shrink=0.7, pad=0.05)
cbar.set_label('Topography (m)', fontweight='bold')

# Path out to save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/figs/urban/paper'
name_out = 'pyplt_maps_study_area_RegCM5_CSAM-3_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

