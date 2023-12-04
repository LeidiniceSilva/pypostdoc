# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 20, 2023"
__description__ = "This script plot dendrogram of the stations"

import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import linkage, dendrogram

path='/marconi/home/userexternal/mdasilva'


def import_inmet():
		
	ix = []		  
	iy = []
	clim_i = []

	# Select lat and lon 
	for station in range(1, 567):

		if station == 15:
			continue
		if station == 23:
			continue
		if station == 47:
			continue
		if station == 105:
			continue
		if station == 112:
			continue
		if station == 117:
			continue
		if station == 124:
			continue
		if station == 137:
			continue
		if station == 149:
			continue
		if station == 158:
			continue
		if station == 174:
			continue
		if station == 183:
			continue
		if station == 335:
			continue
		if station == 343:
			continue
		if station == 359:
			continue
		if station == 398:
			continue
		if station == 399:
			continue
		if station == 413:
			continue
		if station == 417:
			continue
		if station == 422:
			continue
		if station == 426:
			continue
		if station == 444:
			continue
		if station == 453:
			continue
		if station == 457:
			continue
		if station == 458:
			continue
		if station == 479:
			continue
		if station == 490:
			continue
		if station == 495:
			continue
		if station == 505:
			continue
		if station == 529:
			continue
		if station == 566:
			continue
		
		iy.append(inmet[station][2])
		ix.append(inmet[station][3])
				
		print('Reading weather station:', station, inmet[station][0])
		# Reading inmet 
		d_i = xr.open_dataset('{0}/OBS/BDMET/nc/hourly/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[station][0]))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		clim_i.append(values_i*24)
			
	return iy, ix, clim_i


# Import latitude, longitude and database
lat_x, lon_x, clim_tot = import_inmet()
lon_xx = lon_x 
lat_yy = lat_x

df = pd.DataFrame(clim_tot)

print(df)
print()

print(lon_xx)
print(len(lon_xx))
print()
print(lat_yy)
print(len(lat_yy))
print()

# Linkage hierarchical 
z = linkage(df, method='ward', metric='euclidean')

# Agglomerative clustering
Agg_hc = AgglomerativeClustering(n_clusters=5, affinity='euclidean', linkage='ward')
y_hc = Agg_hc.fit_predict(df)
print(y_hc)
print(len(y_hc))

# Plot figure   
plt.figure(figsize=(30,10))
dendrogram(z, leaf_rotation=90., leaf_font_size=10., color_threshold=3800)
plt.title('Dendrogram', fontsize=20) 
plt.xlabel('Weather stations (INMET+SMN)', fontsize=20) 
plt.ylabel('Euclidean distances', fontsize=20) 
# ~ plt.yticks(np.arange(0, 10, 1), fontsize=20)

# Path out to save figure
path_out = '{0}/figs/bdmet'.format(path)
name_out = 'pyplt_dendrogram_stations.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
