# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 16, 2023"
__description__ = "This script plot altimetry of inmet stations "

import os
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_inmet_conv_stations import inmet_conv

path = '/marconi/home/userexternal/mdasilva'

# Select lat and lon 
ix = []		  
iy = []
iz = []

for i in range(1, 141):
	iy.append(inmet_conv[i][2])
	ix.append(inmet_conv[i][3])
	iz.append(inmet_conv[i][4])
		
# Plot figure
fig = plt.figure()
font_size = 8

ax = fig.add_subplot(1, 1, 1)
my_map = Basemap(projection='cyl', llcrnrlon=-85., llcrnrlat=-60., urcrnrlon=-30.,urcrnrlat=15., resolution='c')
my_map.drawmeridians(np.arange(-85.,-30.,10.), labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-60.,15.,10.), labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black', linewidth=.5)
my_map.readshapefile('{0}/github_projects/shp/GEOFT_REGIAO_POLITICO_ADM_2013/GEOFT_REGIAO_POLITICO_ADM_2013'.format(path), 'GEOFT_REGIAO_POLITICO_ADM_2013', drawbounds=True, color='gray', linewidth=1.)

sc = my_map.scatter(ix, iy, 4, iz, cmap='jet', marker='o')
plt.title('(a)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=20, fontsize=font_size, fontweight='bold')
plt.text(-36, -57, u'\u25B2 \nN', fontsize=font_size, fontweight='bold')

# CSAM
a1,b1 = (-78,-35)
a2,b2 = (-78,-11)
a3,b3 = (-35,-11)
a4,b4 = (-35,-35)
poly1 = Polygon([(a1,b1),(a2,b2),(a3,b3),(a4,b4)], facecolor='none', edgecolor='red', linewidth=1.)
plt.gca().add_patch(poly1)

cbar=plt.colorbar(sc, extend='max')
cbar.set_label('Altimetry (meters)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=10)

# Path out to save figure
path_out = '{0}/user/mdasilva/OBS/WS-SA/figs/v2'.format(path)
name_out = 'pyplt_maps_stations_inmet_conv.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

