# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 16, 2023"
__description__ = "This script plot altimetry of inmet stations "

import os
import numpy as np
import matplotlib.pyplot as plt

from dict_inmet_auto_stations import inmet_auto
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap

path = '/marconi/home/userexternal/mdasilva'

skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,7793,100,105,106,107,112,117,124,135,137,139,
149,152,155,158,168,174,177,183,186,199,204,210,212,224,226,239,240,248,249,253,254,276,277,280,293,298,
303,305,306,308,319,334,335,341,343,359,362,364,384,393,396,398,399,400,402,413,416,417,422,423,426,427,
443,444,446,451,453,457,458,467,474,479,483,488,489,490,495,505,509,513,514,516,529,534,544,559,566]

# Select lat and lon 
ix, iy, iz = [], [], []
for i in range(1, 567):
	if i in skip_list:
		continue
	iy.append(inmet_auto[i][2])
	ix.append(inmet_auto[i][3])
	iz.append(inmet_auto[i][4])
					
# Plot figure
fig = plt.figure()
font_size = 8

ax = fig.add_subplot(1, 2, 1)
my_map = Basemap(projection='cyl', llcrnrlon=-85., llcrnrlat=-60., urcrnrlon=-30.,urcrnrlat=15., resolution='c')
my_map.drawmeridians(np.arange(-85.,-30.,10.), labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-60.,15.,10.), labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black', linewidth=.5)

my_map.plot(ix, iy, 'o', color='blue', label='INMET', markersize=2)
plt.title('(a)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=20, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=30, fontsize=font_size, fontweight='bold')
plt.text(-36, -57, u'\u25B2 \nN', fontsize=font_size, fontweight='bold')
plt.legend(loc=1, fontsize=font_size)

# CSAM
a1,b1 = (-78,-35)
a2,b2 = (-78,-11)
a3,b3 = (-35,-11)
a4,b4 = (-35,-35)
poly1 = Polygon([(a1,b1),(a2,b2),(a3,b3),(a4,b4)], facecolor='none', edgecolor='red', linewidth=1.)
plt.gca().add_patch(poly1)

ax = fig.add_subplot(1, 2, 2)
my_map = Basemap(projection='cyl', llcrnrlon=-85., llcrnrlat=-60., urcrnrlon=-30.,urcrnrlat=15., resolution='c')
my_map.drawmeridians(np.arange(-85.,-30.,10.), labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-60.,15.,10.), labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black', linewidth=.5)

sc = my_map.scatter(ix, iy, 4, iz, cmap='jet', label='INMET', marker='o')
plt.title('(b)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=20, fontsize=font_size, fontweight='bold')
plt.text(-36, -57, u'\u25B2 \nN', fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(sc, cax=fig.add_axes([0.91, 0.25, 0.015, 0.50]), extend='max')
cbar.set_label('Altimetry (meters)', fontsize=10, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# CSAM
a1,b1 = (-78,-35)
a2,b2 = (-78,-11)
a3,b3 = (-35,-11)
a4,b4 = (-35,-35)
poly1 = Polygon([(a1,b1),(a2,b2),(a3,b3),(a4,b4)], facecolor='none', edgecolor='red', linewidth=1.)
plt.gca().add_patch(poly1)

# Path out to save figure
path_out = '{0}/user/mdasilva/WS-SA/figs/figs_v2'.format(path)
name_out = 'pyplt_maps_stations_inmet_auto.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

