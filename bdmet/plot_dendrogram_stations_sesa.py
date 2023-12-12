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

path='/afs/ictp.it/home/m/mda_silv/Documents/BDMET'

skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,7793,100,105,106,107,112,117,124,135,137,139,
149,152,155,158,168,174,177,183,186,199,204,210,212,224,226,239,240,248,249,253,254,276,277,280,293,298,
303,305,306,308,319,334,335,341,343,359,362,364,384,393,396,398,399,400,402,413,416,417,422,423,426,427,
443,444,446,451,453,457,458,467,474,479,483,488,489,490,495,505,509,513,514,516,529,534,544,559,566]


def import_inmet():
		
	ix = []		  
	iy = []
	clim_i = []

	# Select lat and lon 
	for station in range(1, 567):
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue
		
		iy.append(inmet[station][2])
		ix.append(inmet[station][3])
				
		print('Reading weather station:', station, inmet[station][0])
		# Reading inmet 
		d_i = xr.open_dataset('{0}/database/nc/hourly/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[station][0]))
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
plt.xlabel('Weather stations (INMET)', fontsize=20) 
plt.ylabel('Euclidean distances', fontsize=20) 

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_dendrogram_stations_sesa.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
