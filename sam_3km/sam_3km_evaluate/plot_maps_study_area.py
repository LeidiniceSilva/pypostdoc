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
dirnc = '/marconi/home/userexternal/mdasilva/user/mdasilva/SAM-3km/input'
domname = 'SAM-3km'

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
lat_end   = 18
lon_start = -90
lon_end   = -30

# Plot study area
fig = plt.figure() 
font_size = 10

my_map = Basemap(llcrnrlon=lon_start, llcrnrlat=lat_start, urcrnrlon=lon_end, urcrnrlat=lat_end, resolution='c')	
my_map.drawparallels(np.arange(lat_start, lat_end, 10.), labels=[1,1,1,1], fontsize=font_size, linewidth=0.4, color='black')
my_map.drawmeridians(np.arange(lon_start, lon_end, 10.), labels=[1,1,1,1], fontsize=font_size, linewidth=0.4, color='black')                  
my_map.readshapefile('/marconi/home/userexternal/mdasilva/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='black', linewidth=0.5)
my_map.readshapefile('/marconi/home/userexternal/mdasilva/github_projects/shp/lim_unid_fed/lim_unid_fed', 'lim_unid_fed', drawbounds=True, color='black', linewidth=0.5)
	
x, y = my_map(lon,lat)
im = my_map.contourf(x, y, topo, np.arange(0, 3050, 50), cmap=plt.cm.terrain, extend='max')
plt.xlabel(u'Longitude', labelpad=20, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=30, fontsize=font_size, fontweight='bold')
plt.text(-56, -39, u'SESA', color='red', fontsize=font_size, fontweight='bold')
plt.text(-36, 11, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')

my_map.plot(-40.3385,-20.3197,'o',ms=3,markeredgewidth=0.75,color='white',mec='black') # Vitoria
my_map.plot(-43.9409,-19.9129,'o',ms=3,markeredgewidth=0.75,color='white',mec='black') # Belo Horizonte
my_map.plot(-54.6464,-20.4428,'o',ms=3,markeredgewidth=0.75,color='white',mec='black') # Campo Grande
my_map.plot(-56.0926,-15.5954,'o',ms=3,markeredgewidth=0.75,color='white',mec='black') # Goiania
my_map.plot(-46.6252,-23.5337,'o',ms=3,markeredgewidth=0.75,color='white',mec='black') # Sao Paulo
my_map.plot(-43.1729,-22.9068,'o',ms=3,markeredgewidth=0.75,color='white',mec='black') # Rio de Janeiro
my_map.plot(-49.2700,-25.4372,'o',ms=3,markeredgewidth=0.75,color='white',mec='black') # Curitiba
my_map.plot(-48.5569,-27.5948,'o',ms=3,markeredgewidth=0.75,color='white',mec='black') # Florianolopis
my_map.plot(-51.2090,-30.0368,'o',ms=3,markeredgewidth=0.75,color='white',mec='black') # Porto Alegre
my_map.plot(-56.1674,-34.8335,'o',ms=3,markeredgewidth=0.75,color='white',mec='black') # Montevideo
my_map.plot(-57.5759,-25.2637,'o',ms=3,markeredgewidth=0.75,color='white',mec='black') # Asuncion
my_map.plot(-58.4004,-34.6051,'o',ms=3,markeredgewidth=0.75,color='white',mec='black') # Buenos Aires
my_map.plot(-68.1193,-16.4897,'o',ms=3,markeredgewidth=0.75,color='white',mec='black') # La Paz
my_map.plot(-70.6693,-33.4489,'o',ms=3,markeredgewidth=0.75,color='white',mec='black') # Santiago

# SESA
c1,d1 = (-65,-35.5)
c2,d2 = (-65,-24)
c3,d3 = (-52,-24)
c4,d4 = (-52,-35.5)
poly2 = Polygon([(c1,d1),(c2,d2),(c3,d3),(c4,d4)], facecolor='none', edgecolor='red', linewidth=1.)
plt.gca().add_patch(poly2)

cbar = fig.colorbar(im, cax=fig.add_axes([0.82, 0.2, 0.026, 0.6]), drawedges=True, fraction=0.030, pad=0.04, aspect=20)
cbar.set_label('Topography (meters)', fontsize=font_size, fontweight='bold')

# Path out to save figure
path_out = '/marconi/home/userexternal/mdasilva/user/mdasilva/SAM-3km/figs/evaluate'
name_out = 'pyplt_maps_study_area_SAM-3km_RegCM5_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

