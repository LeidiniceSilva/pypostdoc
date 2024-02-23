# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jan 02, 2024"
__description__ = "This script plot altimetry of stations"

import os
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap

# Coordinates and altimetry
lat_tot = (-32.53,-31.80,-30.05,-30.01,-27.60,-24.96,-23.72,-23.22,-22.98,-32.08,-25.45,-28.60,-28.86)
lon_tot = (-53.38,-52.41,-51.17,-50.14,-48.62,-48.42,-46.68,-44.73,-42.02,-52.17,-49.23,-53.67,-52.54)
alt_tot = (31.48,13,41.18,4.56,4.87,659.89,771,3,5,4.92,922.91,426.69,660.44)

# Plot study area
fig = plt.figure() 
font_size = 10

lat_start = -35
lat_end   = -10
lon_start = -75
lon_end   = -35

my_map = Basemap(llcrnrlon=lon_start, llcrnrlat=lat_start, urcrnrlon=lon_end, urcrnrlat=lat_end, resolution='c')	
my_map.drawparallels(np.arange(lat_start, lat_end, 5), labels=[1,0,0,0], fontsize=font_size, dashes=[4, 4], linewidth=0.4, color='black')
my_map.drawmeridians(np.arange(lon_start, lon_end, 5), labels=[0,0,0,1], fontsize=font_size, dashes=[4, 4], linewidth=0.4, color='black')                  
my_map.readshapefile('/marconi/home/userexternal/mdasilva/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='gray', linewidth=0.5)
my_map.readshapefile('/marconi/home/userexternal/mdasilva/github_projects/shp/lim_unid_fed/lim_unid_fed', 'lim_unid_fed', drawbounds=True, color='black', linewidth=0.5)

sc=my_map.scatter(lon_tot, lat_tot, 40, alt_tot, cmap='jet', marker='o', vmin=0, vmax=1000)
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=30, fontsize=font_size, fontweight='bold')
plt.text(-38, -33, u'\u25B2 \nN', fontsize=font_size, fontweight='bold')
cbar=plt.colorbar(sc, cax=fig.add_axes([0.15, 0.27, 0.018, 0.45]), extend='max')
cbar.set_label('Altimetry (meters)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '/marconi/home/userexternal/mdasilva/user/mdasilva/SAM-3km-cyclone/figs'
name_out = 'pyplt_maps_stations_alt.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
