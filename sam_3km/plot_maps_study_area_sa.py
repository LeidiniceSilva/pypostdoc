# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Sep 16, 2023"
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
dirnc = '/home/mda_silv/scratch/sam_3km/input'
domname = 'sam-3km-pbl1'

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

lat_start = -58
lat_end   = 14
lon_start = -85
lon_end   = -30

# Plot study area
fig = plt.figure() 
font_size = 10

my_map = Basemap(llcrnrlon=lon_start, llcrnrlat=lat_start, urcrnrlon=lon_end, urcrnrlat=lat_end, resolution='c')	
my_map.drawparallels(np.arange(lat_start, lat_end, 8.), labels=[1,0,0,0], fontsize=font_size, linewidth=0.4, color='black')
my_map.drawmeridians(np.arange(lon_start, lon_end, 10.), labels=[0,0,0,1], fontsize=font_size, linewidth=0.4, color='black')                  
my_map.readshapefile('/home/mda_silv/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='black', linewidth=1.)
x, y = my_map(lon,lat)

llevels = (1, 25, 50, 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000)
im = my_map.contourf(x, y, topo, llevels, cmap=plt.cm.terrain, extend='max')
plt.xlabel(u'Longitude', labelpad=20, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=30, fontsize=font_size, fontweight='bold')
cbar = fig.colorbar(im, drawedges=True, fraction=0.030, pad=0.04, aspect=20)
cbar.set_label('Topography (meters)', fontsize=font_size, fontweight='bold')

plt.text(-36, 8, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')

a1,b1 = (-79.67558,-34.43281)
a2,b2 = (-79.67558,-11.05942)
a3,b3 = (-36.76505,-11.05942)
a4,b4 = (-36.76505,-34.43281)
poly1 = Polygon([(a1,b1),(a2,b2),(a3,b3),(a4,b4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly1)


# Path out to save figure
path_out = '/home/mda_silv/figs/sam_3km'
name_out = 'pyplt_maps_study_area_sam_3km.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
