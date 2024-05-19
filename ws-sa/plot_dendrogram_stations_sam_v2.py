# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "May 17, 2024"
__description__ = "This script plot dendrogram of the stations"

import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import linkage, dendrogram

path='/marconi/home/userexternal/mdasilva/OBS/WS-SA'


def import_inmet():
		
	ix = []		  
	iy = []
	clim_i = []
	
	skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,7793,100,105,106,107,112,117,124,135,137,139,149,152,155,158,168,174,177,183,186,199,204,210,212,224,226,239,240,248,249,253,254,276,277,280,293,298,303,305,306,308,319,334,335,341,343,359,362,364,384,393,396,398,399,400,402,413,416,417,422,423,426,427,443,444,446,451,453,457,458,467,474,479,483,488,489,490,495,505,509,513,514,516,529,534,544,559,566]
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
		d_i = xr.open_dataset('{0}/INMET/nc/hourly/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[station][0]))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		clim_i.append(values_i*24)
			
	return iy, ix, clim_i


def import_smn_i():
		
	ix = []		  
	iy = []
	clim_i = []

	# Select lat and lon 
	for station in range(1, 73):

		iy.append(smn_i[station][1])
		ix.append(smn_i[station][2])
				
		print('Reading weather station:', station, smn_i[station][0])
		# Reading smn 
		d_i = xr.open_dataset('{0}/SMN/hourly/nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[station][0]))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		clim_i.append(values_i*24)
					
	return iy, ix, clim_i


def import_smn_ii():
		
	ix = []		  
	iy = []
	clim_i = []

	skip_list = [47,69,72,80,83,86,87,88,89,90,91,95,96,97,98,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128]
	# Select lat and lon 
	for station in range(1, 129):
		if station in skip_list:
			continue
			
		iy.append(smn_ii[station][1])
		ix.append(smn_ii[station][2])
				
		print('Reading weather station:', station, smn_ii[station][0])
		# Reading smn 
		d_i = xr.open_dataset('{0}/SMN/daily/nc/pre/'.format(path) + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[station][0]))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		clim_i.append(values_i)
		
	return iy, ix, clim_i
	
	
# Import latitude, longitude and database
lat_x, lon_x, clim_x = import_inmet()
lat_y, lon_y, clim_y = import_smn_i()
lat_z, lon_z, clim_z = import_smn_ii()

lat_tot = lat_x + lat_y + lat_z
lon_tot = lon_x + lon_y + lon_z 
clim_tot = clim_x + clim_y + clim_z

df = pd.DataFrame(clim_tot)

print(df)
print()

print(lat_tot)
print(len(lat_tot))
print()
print(lon_tot)
print(len(lon_tot))
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

# Path out to save figure
path_out = '{0}/figs/figs_v2'.format(path)
name_out = 'pyplt_dendrogram_stations_sam.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
