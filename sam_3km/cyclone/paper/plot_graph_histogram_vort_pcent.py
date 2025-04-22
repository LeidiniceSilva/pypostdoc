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

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/SAM-3km'


def read_dat_file(filename):

	data = []
	with open(filename, 'r') as file:
		lines = file.readlines()
		header = []
		rows = []
	
	# Iterate over lines in the file
	for line in lines:
		line = line.strip().split()
		
		if len(line) == 6:
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


def open_dat_file(dataset, yr_init, yr_end):

	vo = []
	pc = []
	for yr in range(yr_init, yr_end+1):
		
		data = read_dat_file('{0}/postproc/cyclone/{1}/track/resultado_{2}.dat'.format(path, dataset, yr))

		rows_list = []
		for i, (header, rows) in enumerate(data):
			if (rows[0][2] < str(-56)):
				rows_list.append(rows[0])
		
		for j  in rows_list:
			vo.append(float(j[3]))
			pc.append(float(j[4]))

	return vo, pc
	
    
# Import model and obs dataset
list_era5_i, list_era5_ii = open_dat_file('ERA5', 2018, 2021)
list_regcm5_i, list_regcm5_ii = open_dat_file('RegCM5', 2018, 2021)
list_wrf415_i, list_wrf415_ii = open_dat_file('WRF415', 2018, 2021)

list_era5_i_ = [x for x in list_era5_i if x != -99.0]
list_regcm5_i_ = [x for x in list_regcm5_i if x != -99.0]
list_wrf415_i_ = [x for x in list_wrf415_i if x != -99.0]

# Plot figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
font_size = 10

ax1.hist([list_era5_i_, list_regcm5_i_, list_wrf415_i_], bins=10, color=['black', 'blue', 'red'], alpha=0.75, edgecolor='black', linewidth=1., label=['ERA5', 'RegCM5', 'WRF415'])
ax1.set_title('(a)', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_xlabel('Relative vorticity (10$^-$$^5$ s$^-$$^1$)', fontsize=font_size, fontweight='bold')
ax1.set_ylabel('Frequency', fontsize=font_size, fontweight='bold')
ax1.set_ylim(0, 20)
ax1.axvspan(-1, -1.5, color='k', alpha=0.1)
ax1.grid(axis='y', c='k', ls='--', alpha=0.3)
ax1.legend(fontsize=font_size, ncol=3, loc=9, shadow=True)

ax2.hist([list_era5_ii, list_regcm5_ii, list_wrf415_ii], bins=10, color=['black', 'blue', 'red'], alpha=0.75, edgecolor='black', linewidth=1., label=['ERA5', 'RegCM5', 'WRF415'])
ax2.set_title('(b)', loc='left', fontsize=font_size, fontweight='bold')
ax2.set_xlabel('Pressure center (hPa)', fontsize=font_size, fontweight='bold')
ax2.set_ylabel('Frequency', fontsize=font_size, fontweight='bold')
ax2.set_ylim(0, 20)
ax2.axvspan(1010, 1018, color='k', alpha=0.1)
ax2.grid(axis='y', c='k', ls='--', alpha=0.3)

# Path out to save figure
path_out = '{0}/figs/cyclone'.format(path)
name_out = 'pyplt_graph_histogram_pcenter_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()




			
