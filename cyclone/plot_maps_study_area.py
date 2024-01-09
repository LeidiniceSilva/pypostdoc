# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jan 02, 2024"
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
dirnc = '/marconi/home/userexternal/mdasilva/user/mdasilva/cyclone/input'
domname = 'SAM-3km-cyclone'

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

lat_start = -60
lat_end   = -10
lon_start = -85
lon_end   = 5

# Plot study area
fig = plt.figure() 
font_size = 10

lat_cyclone_i = (-24,-24,-23,-25,-25,-25,-25,-27,-29,-30,-31,-32,-33,-34,-36,-37,-38,-39,-40,-40,-41,-41,-39,-39,-39,-40,-40,-42,-43,-43,-46,-48,-52,-55.5,-58)
lon_cyclone_i = (-39,-41,-45,-45,-46,-46,-46,-47,-49,-49,-48,-47,-45,-44,-43,-40,-37,-36,-34,-33,-33,-32,-30,-27,-25,-24,-22,-21,-19,-16,-12,-6,-3,3,-1)

lat_cyclone_ii = (-27,-29,-30,-31,-32,-33,-33,-35,-38,-39,-40,-41,-42,-43,-44,-46)
lon_cyclone_ii = (-58,-54,-53,-48,-46,-43,-39,-34,-31,-26,-23,-21,-17,-13,-9,-5)

my_map = Basemap(llcrnrlon=lon_start, llcrnrlat=lat_start, urcrnrlon=lon_end, urcrnrlat=lat_end, resolution='c')	
my_map.drawparallels(np.arange(lat_start, lat_end, 5), labels=[1,0,0,0], fontsize=font_size, dashes=[4, 4], linewidth=0.4, color='black')
my_map.drawmeridians(np.arange(lon_start, lon_end, 10), labels=[0,0,0,1], fontsize=font_size, dashes=[4, 4], linewidth=0.4, color='black')                  
my_map.readshapefile('/marconi/home/userexternal/mdasilva/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='gray', linewidth=0.5)
my_map.readshapefile('/marconi/home/userexternal/mdasilva/github_projects/shp/lim_unid_fed/lim_unid_fed', 'lim_unid_fed', drawbounds=True, color='black', linewidth=0.5)
	
x, y = my_map(lon,lat)

llevels = (0, 25, 50, 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000)
im = my_map.contourf(x, y, topo, llevels, cmap=plt.cm.terrain, extend='max')

plt.plot(lon_cyclone_i, lat_cyclone_i, linewidth=1, color='red', marker='o', markerfacecolor='red', markersize=3, label='Cyclone I')
plt.plot(lon_cyclone_ii, lat_cyclone_ii, linewidth=1, color='black', marker='o', markerfacecolor='black', markersize=3, label='Cyclone II')

plt.xlabel(u'Longitude', labelpad=20, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=30, fontsize=font_size, fontweight='bold')
cbar = fig.colorbar(im, drawedges=True, fraction=0.030, pad=0.04, aspect=20)
cbar.set_label('Topography (meters)', fontsize=font_size, fontweight='bold')

plt.text(-81, -58, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')

plt.legend(loc=1, ncol=1, fontsize=font_size)

# Path out to save figure
path_out = '/marconi/home/userexternal/mdasilva/user/mdasilva/cyclone/figs'
name_out = 'pyplt_maps_study_area_SAM-3km_RegCM5_cyclone_study_case.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

