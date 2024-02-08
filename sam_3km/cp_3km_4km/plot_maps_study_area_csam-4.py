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
dirnc = '/marconi_work/ICT23_ESP/fraffael/REGCM5/output-py/pbl1/output/CSAM-4/ICTP/ECMWF-ERA5/evaluation/r1i1p1f1/ICTP-RegCM5/v0/fixed/orog'
domname = 'orog_CSAM-4'

# RegCM file
if len(sys.argv) > 1:
    RCMf = nc(sys.argv[1], mode='r')
else:
    RCMf = nc(os.path.join(dirnc,domname+'_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5_v0_fixed_2018010100-2018020100.nc'), mode='r')
    
lat  = RCMf.variables['lat'][:,:]
lon  = RCMf.variables['lon'][:,:]
topo = RCMf.variables['orog'][:,:]
lonc = RCMf.longitude_of_projection_origin
latc = RCMf.latitude_of_projection_origin
RCMf.close()

lat_start = -58
lat_end   = 18
lon_start = -90
lon_end   = -30

# Plot study area
fig = plt.figure() 
font_size = 10

my_map = Basemap(llcrnrlon=lon_start, llcrnrlat=lat_start, urcrnrlon=lon_end, urcrnrlat=lat_end, resolution='c')	
my_map.drawparallels(np.arange(lat_start, lat_end, 10.), labels=[1,0,0,0], fontsize=font_size, linewidth=0.4, color='black')
my_map.drawmeridians(np.arange(lon_start, lon_end, 10.), labels=[0,0,0,1], fontsize=font_size, linewidth=0.4, color='black')                  
my_map.readshapefile('/marconi/home/userexternal/mdasilva/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='black', linewidth=0.5)
my_map.readshapefile('/marconi/home/userexternal/mdasilva/github_projects/shp/lim_unid_fed/lim_unid_fed', 'lim_unid_fed', drawbounds=True, color='black', linewidth=0.5)
x, y = my_map(lon,lat)
 
im = my_map.contourf(x, y, topo, np.arange(0, 3000, 50), cmap=plt.cm.terrain, extend='max')
plt.xlabel(u'Longitude', labelpad=20, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=30, fontsize=font_size, fontweight='bold')
plt.text(-56, -39, u'SESA', color='red', fontsize=font_size, fontweight='bold')
plt.text(-36, 11, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
cbar = fig.colorbar(im, drawedges=True, fraction=0.030, pad=0.04, aspect=20)
cbar.set_label('Topography (meters)', fontsize=font_size, fontweight='bold')

# CSAM
# my_map.drawgreatcircle(-78.4, -34.5, -75.0, -11.2, alpha=.6, linewidth=4, color='gray')
# my_map.drawgreatcircle(-75.0, -11.0, -38.0, -11.0, alpha=.6, linewidth=4, color='gray')
# my_map.drawgreatcircle(-35.0, -34.5, -38.0, -11.2, alpha=.6, linewidth=4, color='gray')
# my_map.drawgreatcircle(-78.4, -34.4, -36.0, -34.4, alpha=.6, linewidth=4, color='gray') 

# SESA
c1,d1 = (-65,-35)
c2,d2 = (-65,-24)
c3,d3 = (-52,-24)
c4,d4 = (-52,-35)
poly2 = Polygon([(c1,d1),(c2,d2),(c3,d3),(c4,d4)], facecolor='none', edgecolor='red', linewidth=1.)
plt.gca().add_patch(poly2)

# Path out to save figure
path_out = '/marconi/home/userexternal/mdasilva/user/mdasilva/SAM-3km/figs/cp_3km-4km'
name_out = 'pyplt_maps_study_area_CP-RegCM5_CSAM-4_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

