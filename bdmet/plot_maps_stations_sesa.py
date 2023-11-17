# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 16, 2023"
__description__ = "This script plot weather stations dots"

import os
import numpy as np
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap

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
		
# Plot figure
fig = plt.figure()

my_map = Basemap(projection='cyl', llcrnrlon=-90., llcrnrlat=-60., urcrnrlon=-30.,urcrnrlat=20., resolution='c')
my_map.drawmeridians(np.arange(-90.,-20.,10.), labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-60.,30.,10.), labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.plot(ix, iy, 'o', color='blue', label='INMET', markersize=2)
my_map.plot(jx, jy, 'o', color='gray', label='SMN', markersize=2)
my_map.plot(kx, ky, 'o', color='gray', markersize=2)
plt.legend(loc=1, fontsize=10)

my_map.readshapefile('/afs/ictp.it/home/m/mda_silv/Documents/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='black', linewidth=.5)

plt.xlabel(u'Longitude', labelpad=20, fontsize=10, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=30, fontsize=10, fontweight='bold')
plt.text(-56, -44, u'CSAM', fontsize=10, fontweight='bold')
plt.text(-36, -57, u'\u25B2 \nN', fontsize=10, fontweight='bold')

# CSAM
a1,b1 = (-70,-40)
a2,b2 = (-70,-15)
a3,b3 = (-45,-15)
a4,b4 = (-45,-40)
poly1 = Polygon([(a1,b1),(a2,b2),(a3,b3),(a4,b4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly1)

# Path out to save figure
path_out = '/afs/ictp.it/home/m/mda_silv/Documents/FPS_SESA/figs/sesa_v2'
name_out = 'pyplt_maps_stations_sesa.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

