# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot pressure center of EC over SAM"

import os
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from netCDF4 import Dataset
from datetime import datetime

dataset = 'RegCM5'
path = '/marconi/home/userexternal/mdasilva/user/mdasilva/SAM-3km/post_cyclone/ECyclone'


def read_dat_file(filename):

	data = []
	with open(filename, 'r') as file:
		lines = file.readlines()
		header = []
		rows = []
	
	# Iterate over lines in the file
	for line in lines:
		line = line.strip().split()
		
		# If the line contains 6 elements, it's considered a header
		if len(line) == 6:
			
			# If we have rows, append them to data
			if rows:
				data.append((header, rows))
				rows = []
			header = line
		else:  
			rows.append(line)
	
	# Append the last header and rows to data
	if header and rows:
		data.append((header, rows))
	
	return data


def open_dat_file(dataset):

	vo = []
	for yr in range(2018, 2021+1):
	
		data = read_dat_file('{0}/ECyclone_{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
	
		rows_list = []
		rows_list_i = []
		for i, (header, rows) in enumerate(data):
			rows_list.append(rows[0])
		
		for j  in rows_list:
			vo.append(float(j[4]))
	
	return vo


def plot_dual_histogram(list1, list2, bins=10, xlabel='', ylabel=''):
    plt.hist([list1, list2], bins=bins, color=['red', 'blue'], alpha=0.75, edgecolor='black', label=['ERA5', 'RegCM5'])
    plt.xlabel(xlabel, fontweight='bold')
    plt.ylabel(ylabel, fontweight='bold')
    plt.legend()
    plt.show()
     

# Import model and obs dataset
list_era5 = open_dat_file('ERA5')
list_regcm5 = open_dat_file('RegCM5')

era5 = [1017.0, 1018.0, 1016.0, 1017.0, 1018.0, 1018.0, 1020.0, 1020.0, 1021.0, 1021.0, 1016.0, 1017.0, 1020.0, 1020.0, 1018.0, 1020.0, 1019.0, 1018.0, 1019.0, 1018.0, 1016.0, 1013.0, 1023.0, 1017.0, 1017.0, 1016.0, 1016.0, 1016.0, 1014.0, 1016.0, 1012.0, 1015.0, 1014.0, 1018.0, 1016.0, 1018.0, 1021.0, 1021.0, 1020.0, 1019.0, 1018.0, 1018.0, 1021.0, 1012.0, 1021.0, 1015.0, 1021.0, 1018.0, 1018.0, 1021.0, 1016.0, 1015.0, 1014.0, 1016.0, 1013.0, 1016.0, 1017.0, 1018.0, 1019.0, 1017.0, 1018.0, 1018.0, 1019.0, 1025.0, 1021.0, 1018.0, 1017.0, 1018.0, 1016.0, 1018.0, 1017.0, 1011.0, 1020.0, 1018.0, 1016.0, 1016.0, 1019.0, 1014.0, 1013.0, 1013.0, 1013.0, 1012.0, 1016.0, 1012.0, 1012.0, 1016.0, 1016.0, 1015.0, 1019.0, 1019.0, 1020.0, 1013.0, 1019.0, 1017.0, 1017.0, 1019.0, 1015.0, 1012.0]
regcm5 = [1013.0, 1010.0, 1012.0, 1015.0, 1015.0, 1017.0, 1019.0, 1019.0, 1020.0, 1019.0, 1019.0, 1016.0, 1021.0, 1020.0, 1016.0, 1018.0, 1017.0, 1022.0, 1014.0, 1015.0, 1017.0, 1016.0, 1021.0, 1019.0, 1013.0, 1016.0, 1010.0, 1015.0, 1015.0, 1010.0, 1013.0, 1014.0, 1016.0, 1013.0, 1012.0, 1010.0, 1014.0, 1016.0, 1015.0, 1015.0, 1015.0, 1017.0, 1019.0, 1017.0, 1015.0, 1019.0, 1017.0, 1016.0, 1019.0, 1016.0, 1012.0, 1018.0, 1018.0, 1018.0, 1016.0, 1019.0, 1020.0, 1019.0, 1018.0, 1017.0, 1017.0, 1017.0, 1015.0, 1017.0, 1015.0, 1020.0, 1017.0, 1017.0, 1012.0, 1016.0, 1013.0, 1013.0, 1016.0, 1015.0, 1017.0, 1017.0, 1017.0, 1017.0, 1017.0, 1020.0, 1018.0, 1020.0, 1022.0, 1017.0, 1022.0, 1018.0, 1019.0, 1018.0, 1016.0, 1018.0, 1018.0, 1018.0, 1016.0, 1016.0, 1015.0, 1017.0, 1012.0, 1014.0, 1015.0, 1012.0, 1010.0, 1012.0, 1014.0, 1014.0, 1017.0, 1016.0, 1017.0, 1020.0, 1022.0, 1018.0, 1020.0, 1021.0, 1018.0, 1016.0, 1020.0, 1022.0, 1021.0, 1017.0, 1018.0, 1016.0, 1015.0, 1020.0, 1022.0, 1020.0, 1016.0, 1018.0, 1012.0, 1014.0, 1015.0, 1014.0]

# Plot figure
fig, ax = plt.subplots()
font_size = 10

plt.hist([era5, regcm5], bins=10, color=['red', 'blue'], alpha=0.75, edgecolor='black', linewidth=1.5, label=['ERA5', 'RegCM5'])
ax.set_xlabel('Pressure center hPa', fontsize=font_size, fontweight='bold')
ax.set_ylabel('Frequency', fontsize=font_size, fontweight='bold')
ax.legend(fontsize=font_size, ncol=1, loc=2, shadow=True)
ax.spines[['right', 'top']].set_visible(False)

# Path out to save figure
path_out = '/marconi/home/userexternal/mdasilva/user/mdasilva/SAM-3km/figs/cyclone'
name_out = 'pyplt_graph_histogram_pcen_EC_ERA5_RegCM5_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()




			
