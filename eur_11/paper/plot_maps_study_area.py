# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot study area"

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from matplotlib import gridspec
from matplotlib.path import Path
from netCDF4 import Dataset as nc
from matplotlib.patches import Polygon
from matplotlib.patches import PathPatch

# Specify directories 
dirnc = '/leonardo/home/userexternal/mdasilva/leonardo_scratch/EUR-11/icbc'

domname = 'EUR-11'
dt = '2000-2009'

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

lat_start = np.min(lat)
lat_end   = np.max(lat)
lon_start = np.min(lon)
lon_end   = np.max(lon)

# Plot study area
fig = plt.figure() 
font_size = 10

ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([lon_start, lon_end, lat_start, lat_end], crs=ccrs.PlateCarree())

gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='black', alpha=0.75, linestyle='--')
gl.top_labels = False
gl.bottom_labels = True
gl.left_labels = True
gl.right_labels = False
gl.xlabel_style = {'size': font_size}
gl.ylabel_style = {'size': font_size}
gl.xlocator = plt.FixedLocator(np.arange(lon_start, lon_end, 10.))
gl.ylocator = plt.FixedLocator(np.arange(lat_start, lat_end, 10.))

ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='black')
ax.add_feature(cfeature.COASTLINE, linewidth=0.5, edgecolor='black')

# Plot the topography
im = ax.contourf(lon, lat, topo, np.arange(0, 2900, 50), cmap=plt.cm.terrain, extend='neither', transform=ccrs.PlateCarree())
plt.text(-40, 25, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold', transform=ccrs.PlateCarree())

ax.plot([1, 17, 17, 1, 1],
        [40, 40, 50, 50, 40],
        color='red', linewidth=1, transform=ccrs.PlateCarree())

cbar = fig.colorbar(im, cax=fig.add_axes([0.92, 0.2, 0.026, 0.6]), drawedges=True, fraction=0.030, pad=0.04, aspect=20)
cbar.set_label('Topography (meters)', fontsize=font_size, fontweight='bold')

# Path out to save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11/figs/paper'
name_out = 'pyplt_maps_study_area_{0}_RegCM5_{1}.png'.format(domname, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
