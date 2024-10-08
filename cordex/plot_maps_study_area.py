# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 16, 2023"
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
dirnc = '/marconi/home/userexternal/mdasilva/user/mdasilva/CORDEX/ERA5/icbc'

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

topo = topo[30:-30, 30:-30]

lat_start = -58
lat_end   = 18
lon_start = -90
lon_end   = -30

# Plot study area
fig = plt.figure() 
font_size = 10

# Creating mask of the border
border_mask = np.full((783, 1231), np.nan)
border_mask[:29, :] = 1
border_mask[-29:, :] = 1
border_mask[:, :29] = 1
border_mask[:, -29:] = 1

my_map = Basemap(llcrnrlon=lon_start, llcrnrlat=lat_start, urcrnrlon=lon_end, urcrnrlat=lat_end, resolution='c')	
my_map.drawparallels(np.arange(lat_start, lat_end, 10.), labels=[1,1,1,1], fontsize=font_size, linewidth=0.4, color='black')
my_map.drawmeridians(np.arange(lon_start, lon_end, 10.), labels=[1,1,1,1], fontsize=font_size, linewidth=0.4, color='black')                  
my_map.readshapefile('/marconi/home/userexternal/mdasilva/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='black', linewidth=0.5)
my_map.readshapefile('/marconi/home/userexternal/mdasilva/github_projects/shp/lim_unid_fed/lim_unid_fed', 'lim_unid_fed', drawbounds=True, color='black', linewidth=0.5)

# Plot the topography
xx, yy = my_map(lon,lat)

x_mod = xx[30:-30,30:-30]
y_mod = yy[30:-30,30:-30]

# Mask the border points with gray color
#masked_topo_array = np.ma.masked_array(topo.copy(), mask=border_mask)
cs_gray = my_map.contourf(xx, yy, border_mask, cmap='gray', levels=[0, 1])

im = my_map.contourf(x_mod, y_mod, topo, np.arange(0, 3050, 50), cmap=plt.cm.terrain, extend='max')
plt.xlabel(u'Longitude', labelpad=20, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=30, fontsize=font_size, fontweight='bold')
plt.text(-36, 11, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')

my_map.plot(-47.9000, -16.0000, 'o', ms=4, markeredgewidth=0.75, color='white', mec='black') # Brazilia
my_map.plot(-56.1674, -34.8335, 'o', ms=4, markeredgewidth=0.75, color='white', mec='black') # Montevideo
my_map.plot(-57.5759, -25.2637, 'o', ms=4, markeredgewidth=0.75, color='white', mec='black') # Asuncion
my_map.plot(-58.4004, -34.6051, 'o', ms=4, markeredgewidth=0.75, color='white', mec='black') # Buenos Aires
my_map.plot(-68.1193, -16.4897, 'o', ms=4, markeredgewidth=0.75, color='white', mec='black') # La Paz
my_map.plot(-70.6693, -33.4489, 'o', ms=4, markeredgewidth=0.75, color='white', mec='black') # Santiago

cbar = fig.colorbar(im, cax=fig.add_axes([0.82, 0.2, 0.026, 0.6]), drawedges=True, fraction=0.030, pad=0.04, aspect=20)
cbar.set_label('Topography (meters)', fontsize=font_size, fontweight='bold')

# Path out to save figure
path_out = '/marconi/home/userexternal/mdasilva/user/mdasilva/CORDEX/figs'
name_out = 'pyplt_maps_study_area_RegCM5_CSAM-3_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

