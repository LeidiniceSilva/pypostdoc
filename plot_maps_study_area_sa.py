# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 16, 2023"
__description__ = "This script plot study area"

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import gridspec
from matplotlib.path import Path
from netCDF4 import Dataset as nc
from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii
from matplotlib.patches import Polygon
from matplotlib.patches import PathPatch
from mpl_toolkits.basemap import Basemap, cm

# Select lat and lon 
ix = []		  
iy = []
jx = []
jy = []
kx = []
ky = []

for i in range(1, 100):

	ix.append(inmet[i][3])
	iy.append(inmet[i][2])

for j in range(1, 73):

	jx.append(smn_i[j][2])
	jy.append(smn_i[j][1])
	
for k in range(1, 87):

	kx.append(smn_ii[k][2])
	ky.append(smn_ii[k][1])
	
# Specify directories 
dirnc = '/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp'
domname = 'orog_CSAM-4i_ECMWF-ERA5'

# RegCM file
if len(sys.argv) > 1:
    RCMf = nc(sys.argv[1], mode='r')
else:
    RCMf = nc(os.path.join(dirnc,domname+'_evaluation_r1i1p1f1-USP-RegCM471_v0.nc'), mode='r')
    
lat  = RCMf.variables['xlat'][:,:]
lon  = RCMf.variables['xlon'][:,:]
topo = RCMf.variables['topo'][:,:]
lonc = RCMf.longitude_of_projection_origin
latc = RCMf.latitude_of_projection_origin
RCMf.close()

lat_start  = -35
lat_end    = -17
lon_start  = -75
lon_end    = -48

# Plot study area
fig = plt.figure(figsize=(10, 8)) 
gs = gridspec.GridSpec(1, 2, width_ratios=[1.5, 3]) 
font_size = 10

ax = plt.subplot(gs[0])
my_map = Basemap(ax=ax, llcrnrlon=-82., llcrnrlat=-56, urcrnrlon=-34., urcrnrlat=12, resolution='c', area_thresh=10000., projection='cyl', lon_0=lonc, lat_0=latc, lat_ts=0)	
my_map.drawparallels(np.arange(-56., 12,  10.), labels=[1,0,0,0], fontsize=font_size, linewidth=1., color='black')
my_map.drawmeridians(np.arange(-82, -34, 10.), labels=[0,0,0,1], fontsize=font_size, linewidth=1., color='black')                  
my_map.readshapefile('/home/nice/Documentos/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='black', linewidth=1.)

my_map.plot(ix, iy, 'o', color='blue', label='INMET', markersize=2)
my_map.plot(jx, jy, 'o', color='gray', label='SMN', markersize=2)	
my_map.plot(kx, ky, 'o', color='gray', markersize=2)	
plt.xlabel(u'Longitude', labelpad=20, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=30, fontsize=font_size, fontweight='bold')
plt.text(-40, 6, u'\u25B2 \nN', fontsize=font_size, fontweight='bold')
# ~ plt.text(-56, -39, u'SESA', color='red', fontsize=font_size, fontweight='bold')
plt.text(-56, -44, u'CSAM', color='black', fontsize=font_size, fontweight='bold')
plt.legend(loc=4, fontsize=font_size, frameon=False)

# CSAM
a1,b1 = (-70,-40)
a2,b2 = (-70,-15)
a3,b3 = (-45,-15)
a4,b4 = (-45,-40)
poly1 = Polygon([(a1,b1),(a2,b2),(a3,b3),(a4,b4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly1)

# ~ # SESA
# ~ a1,b1 = (-65,-35)
# ~ a2,b2 = (-65,-17)
# ~ a3,b3 = (-48,-17)
# ~ a4,b4 = (-48,-35)
# ~ poly1 = Polygon([(a1,b1),(a2,b2),(a3,b3),(a4,b4)], facecolor='none', edgecolor='red', linewidth=1.)
# ~ plt.gca().add_patch(poly1)

ax = plt.subplot(gs[1])
my_map = Basemap(ax=ax, llcrnrlon=lon_start, llcrnrlat=lat_start, urcrnrlon=lon_end, urcrnrlat=lat_end, resolution='c', area_thresh=10000., projection='cyl', lon_0=lonc, lat_0=latc, lat_ts=0)	
my_map.drawparallels(np.arange(lat_start, lat_end,  5.), labels=[1,0,0,0], fontsize=font_size, linewidth=1., color='black')
my_map.drawmeridians(np.arange(lon_start, lon_end, 5.), labels=[0,0,0,1], fontsize=font_size, linewidth=1., color='black')                  
my_map.readshapefile('/home/nice/Documentos/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='black', linewidth=1.)
x, y = my_map(lon,lat)

llevels = (1, 25, 50, 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000)
im = my_map.contourf(x, y, topo, llevels, cmap=plt.cm.terrain, extend='max')
plt.xlabel(u'Longitude', labelpad=20, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=30, fontsize=font_size, fontweight='bold')
cbar = fig.colorbar(im, drawedges=True, fraction=0.030, pad=0.04, aspect=20)
cbar.set_label('Topography (meters)', fontsize=font_size, fontweight='bold')

# ~ # SESA
# ~ a1,b1 = (-65,-34.9)
# ~ a2,b2 = (-65,-17.1)
# ~ a3,b3 = (-48.1,-17.1)
# ~ a4,b4 = (-48.1,-34.9)
# ~ poly1 = Polygon([(a1,b1),(a2,b2),(a3,b3),(a4,b4)], facecolor='none', edgecolor='red', linewidth=2.)
# ~ plt.gca().add_patch(poly1)

# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/paper_cp'
name_out = 'pyplt_maps_study_area_sesa.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
