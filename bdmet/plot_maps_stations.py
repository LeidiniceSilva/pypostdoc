# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 16, 2023"
__description__ = "This script plot inmet stations dots"

import os
import numpy as np
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap

path = '/home/nice/Documentos'

# Select lat and lon 
ix = []		  
iy = []

for i in range(1, 567):

	iy.append(inmet[i][2])
	ix.append(inmet[i][3])
		
# Plot figure
fig = plt.figure()

my_map = Basemap(projection='cyl', llcrnrlon=-85., llcrnrlat=-60., urcrnrlon=-30.,urcrnrlat=15., resolution='c')
my_map.drawmeridians(np.arange(-85.,-30.,10.), labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-60.,15.,10.), labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.plot(ix, iy, 'o', color='blue', label='INMET', markersize=2)
plt.legend(loc=1, fontsize=10)

my_map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black', linewidth=.5)

plt.xlabel(u'Longitude', labelpad=20, fontsize=10, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=30, fontsize=10, fontweight='bold')
plt.text(-36, -57, u'\u25B2 \nN', fontsize=10, fontweight='bold')

# Path out to save figure
path_out = '{0}/FPS_SESA/figs/bdmet'.format(path)
name_out = 'pyplt_maps_stations_dots.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

