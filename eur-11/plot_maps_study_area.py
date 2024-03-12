# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot study area"

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import gridspec
from matplotlib.path import Path
from netCDF4 import Dataset as nc
from matplotlib.patches import Polygon
from matplotlib.patches import PathPatch
from mpl_toolkits.basemap import Basemap, cm

# Specify directories 
dirnc = '/marconi/home/userexternal/mdasilva/user/mdasilva/EUR-11/icbc'

domname = 'EUR-11'
dt = '2000-2000'

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

lat_start = 15
lat_end   = 75
lon_start = -45
lon_end   = 65

# Plot study area
fig = plt.figure() 
font_size = 10

# Creating mask of the border
border_mask = np.full((415, 427), np.nan)
border_mask[:39, :] = 1
border_mask[-39:, :] = 1
border_mask[:, :39] = 1
border_mask[:, -39:] = 1

my_map = Basemap(llcrnrlon=lon_start, llcrnrlat=lat_start, urcrnrlon=lon_end, urcrnrlat=lat_end, resolution='c')	
my_map.drawparallels(np.arange(lat_start, lat_end, 10.), labels=[1,1,1,1], fontsize=font_size, linewidth=0.4, color='black')
my_map.drawmeridians(np.arange(lon_start, lon_end, 10.), labels=[1,1,1,1], fontsize=font_size, linewidth=0.4, color='black')                  
my_map.readshapefile('/marconi/home/userexternal/mdasilva/github_projects/shp/shp_world/world', 'world', drawbounds=True, color='black', linewidth=0.5)

# Plot the topography
xx, yy = my_map(lon,lat)

im = my_map.contourf(xx, yy, topo, np.arange(0, 2900, 50), cmap=plt.cm.terrain, extend='max')
plt.xlabel(u'Longitude', labelpad=20, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=30, fontsize=font_size, fontweight='bold')
plt.text(-36, 16, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')

# Mask the border points with gray color
cs_gray = my_map.contourf(xx, yy, border_mask, cmap='gray', levels=[0, 1])

cbar = fig.colorbar(im, cax=fig.add_axes([0.999, 0.2, 0.026, 0.6]), drawedges=True, fraction=0.030, pad=0.04, aspect=20)
cbar.set_label('Topography (meters)', fontsize=font_size, fontweight='bold')

# Path out to save figure
path_out = '/marconi/home/userexternal/mdasilva/user/mdasilva/EUR-11/figs'
name_out = 'pyplt_maps_study_area_{0}_RegCM5_{1}.png'.format(domname, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

