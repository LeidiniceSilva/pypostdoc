# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot vorticity of EC over SAM"

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
			vo.append(float(j[3]))
	
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

era5 = [-1.39, -1.1, -1.48, -1.31, -1.04, -1.52, -1.63, -1.69, -1.61, -1.61, -1.07, -1.13, -2.77, -1.8, -1.64, -1.31, -1.28, -1.3, -1.27, -1.35, -1.16, -1.03, -1.46, -1.38, -1.21, -1.22, -1.59, -1.04, -1.63, -1.01, -1.07, -1.16, -2.84, -1.9, -1.3, -1.73, -1.54, -2.89, -2.58, -2.04, -1.29, -1.6, -2.31, -1.37, -2.09, -1.19, -1.15, -1.21, -1.2, -1.88, -1.28, -2.13, -1.76, -2.4, -1.02, -1.37, -1.47, -1.08, -1.1, -1.35, -1.46, -1.39, -1.54, -1.26, -1.46, -2.01, -1.63, -1.71, -1.36, -1.57, -1.79, -2.44, -1.77, -2.16, -1.29, -1.4, -1.21, -1.49, -1.23, -1.05, -1.52, -1.6, -1.46, -1.5, -1.1, -1.71, -1.82, -1.56, -1.93, -1.14, -1.54, -1.28, -1.94, -1.73, -1.46, -1.52, -1.5, -1.67]
regcm5 = [-1.78, -2.04, -1.35, -1.44, -1.99, -1.92, -1.86, -2.29, -2.08, -1.68, -2.14, -1.95, -2.22, -1.73, -1.34, -1.56, -1.17, -1.37, -1.2, -1.27, -1.28, -1.25, -1.59, -2.31, -1.83, -1.35, -1.34, -1.43, -1.83, -1.72, -1.07, -1.6, -1.05, -1.69, -1.4, -1.34, -1.44, -1.29, -1.61, -2.31, -1.89, -2.18, -2.35, -2.09, -1.21, -1.09, -2.08, -1.05, -1.14, -1.97, -1.22, -1.15, -1.23, -2.02, -1.85, -1.79, -3.39, -1.82, -1.46, -1.01, -1.18, -1.82, -1.48, -1.93, -2.48, -2.58, -1.59, -1.69, -1.41, -1.54, -2.3, -2.36, -1.27, -1.54, -1.62, -1.55, -1.48, -1.46, -1.04, -1.1, -1.27, -1.3, -1.76, -1.12, -1.84, -2.15, -1.22, -1.6, -1.44, -2.42, -1.7, -1.8, -2.52, -1.73, -1.78, -1.3, -1.39, -1.32, -1.25, -1.29, -1.5, -1.79, -1.62, -1.52, -1.74, -2.23, -1.78, -1.02, -1.61, -1.22, -1.23, -1.46, -1.85, -1.45, -2.2, -1.81, -1.17, -2.21, -2.31, -1.42, -1.5, -1.06, -2.39, -2.36, -2.12, -2.11, -1.6, -1.2, -1.44, -1.26]

# Plot figure
fig, ax = plt.subplots()
font_size = 10

plt.hist([era5, regcm5], bins=10, color=['red', 'blue'], alpha=0.75, edgecolor='black', linewidth=1.5, label=['ERA5', 'RegCM5'])
ax.set_xlabel('Relative vorticity 10$^-$$^5$ s$^-$$^1$', fontsize=font_size, fontweight='bold')
ax.set_ylabel('Frequency', fontsize=font_size, fontweight='bold')
ax.legend(fontsize=font_size, ncol=1, loc=2, shadow=True)
ax.spines[['right', 'top']].set_visible(False)

# Path out to save figure
path_out = '/marconi/home/userexternal/mdasilva/user/mdasilva/SAM-3km/figs/cyclone'
name_out = 'pyplt_graph_histogram_vorticity_EC_ERA5_RegCM5_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()




			
