# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot annual cycle of EC over SAM"

import os
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from netCDF4 import Dataset
from datetime import datetime
from scipy.stats.stats import pearsonr

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

	header_list_i, header_list_ii = [], []
	for yr in range(yr_init, yr_end+1):
	
		data = read_dat_file('{0}/postproc/cyclone/{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
		for i, (header, rows) in enumerate(data):
			if (rows[0][2] < str(-56)):			
				header_list_i.append(header[1])
				header_list_ii.append(header)
	
	return header_list_i, header_list_ii
	

def cyclone_number(list_dataset):

	jan_dt, feb_dt, mar_dt, apr_dt, may_dt, jun_dt, jul_dt, aug_dt, sep_dt, oct_dt, nov_dt, dec_dt = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
	for ii  in list_dataset:
		date_obj = datetime.strptime(ii, '%Y%m%d%H')
		if date_obj.month == 1:
               		jan_dt += 1
		if date_obj.month == 2:
             	   feb_dt += 1
		if date_obj.month == 3:
                	mar_dt += 1
		if date_obj.month == 4:
                	apr_dt += 1
		if date_obj.month == 5:
               		may_dt += 1
		if date_obj.month == 6:
                	jun_dt += 1
		if date_obj.month == 7:
                	jul_dt += 1
		if date_obj.month == 8:
                	aug_dt += 1
		if date_obj.month == 9:
                	sep_dt += 1
		if date_obj.month == 10:
                	oct_dt += 1
		if date_obj.month == 11:
                	nov_dt += 1
		if date_obj.month == 12:
                	dec_dt += 1
	
	list_numbers = jan_dt, feb_dt, mar_dt, apr_dt, may_dt, jun_dt, jul_dt, aug_dt, sep_dt, oct_dt, nov_dt, dec_dt
	
	return list_numbers
    

def cyclone_lifetime(list_dataset):
	
	jan_dt, feb_dt, mar_dt, apr_dt, may_dt, jun_dt, jul_dt, aug_dt, sep_dt, oct_dt, nov_dt, dec_dt = [], [], [], [], [], [], [], [], [], [], [], []
	for ii  in list_dataset:
		date_obj = datetime.strptime(ii[1], '%Y%m%d%H')
		if date_obj.month == 1:
			jan_dt.append(int(ii[4]))
		if date_obj.month == 2:
			feb_dt.append(int(ii[4]))
		if date_obj.month == 3:
                	mar_dt.append(int(ii[4]))
		if date_obj.month == 4:
                	apr_dt.append(int(ii[4]))
		if date_obj.month == 5:
               		may_dt.append(int(ii[4]))
		if date_obj.month == 6:
                	jun_dt.append(int(ii[4]))
		if date_obj.month == 7:
                	jul_dt.append(int(ii[4]))
		if date_obj.month == 8:
                	aug_dt.append(int(ii[4]))
		if date_obj.month == 9:
                	sep_dt.append(int(ii[4]))
		if date_obj.month == 10:
                	oct_dt.append(int(ii[4]))
		if date_obj.month == 11:
                	nov_dt.append(int(ii[4]))
		if date_obj.month == 12:
                	dec_dt.append(int(ii[4]))
	
	jan_dur = np.sum(jan_dt) / 18
	feb_dur = np.sum(feb_dt) / 18
	mar_dur = np.sum(mar_dt) / 18
	apr_dur = np.sum(apr_dt) / 18
	may_dur = np.sum(may_dt) / 18
	jun_dur = np.sum(jun_dt) / 18
	jul_dur = np.sum(jul_dt) / 18
	aug_dur = np.sum(aug_dt) / 18
	sep_dur = np.sum(sep_dt) / 18
	oct_dur = np.sum(oct_dt) / 18
	nov_dur = np.sum(nov_dt) / 18
	dec_dur = np.sum(dec_dt) / 18

	list_lifetime = jan_dur, feb_dur, mar_dur, apr_dur, may_dur, jun_dur, jul_dur, aug_dur, sep_dur, oct_dur, nov_dur, dec_dur
	
	return list_lifetime
	
	
# Import model and obs dataset
list_era5_i, list_era5_ii = open_dat_file('ERA5', 2018, 2021)
list_regcm5_i, list_regcm5_ii = open_dat_file('RegCM5', 2018, 2021)
list_wrf415_i, list_wrf415_ii = open_dat_file('WRF415', 2018, 2021)

list_number_era5 = cyclone_number(list_era5_i)
list_number_regcm5 = cyclone_number(list_regcm5_i)
list_number_wrf415 = cyclone_number(list_wrf415_i)

print(pearsonr(list_number_era5, list_number_regcm5))
print(pearsonr(list_number_era5, list_number_wrf415))
exit()

list_lifetime_era5 = cyclone_lifetime(list_era5_ii)
list_lifetime_regcm5 = cyclone_lifetime(list_regcm5_ii)
list_lifetime_wrf415 = cyclone_lifetime(list_wrf415_ii)

# Plot figure
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
ax1, ax2 = axes

labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
x = np.arange(len(labels))  
font_size = 10
width = 0.25 

ax1.set_facecolor('lightgray')
ax1.bar(x-width, list_number_era5, width, color='black', alpha=0.75, edgecolor='black', linewidth=1.25, label='ERA5')
ax1.bar(x, list_number_regcm5, width, color='blue', alpha=0.75, edgecolor='black', linewidth=1.25, label='RegCM5')
ax1.bar(x+width, list_number_wrf415, width, color='red', alpha=0.75, edgecolor='black', linewidth=1.25, label='WRF415')
ax1.set_title('(a)', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_ylabel('Number of cyclones', fontsize=font_size, fontweight='bold')
ax1.set_xlabel('Months', fontsize=font_size, fontweight='bold')
ax1.set_xticks(x)
ax1.set_xticklabels(labels, fontsize=font_size)
ax1.set_ylim(0, 20)
ax1.set_yticks(np.arange(0, 22, 2))
ax1.grid(axis='y', c='k', ls='--', alpha=0.5)
ax1.legend(fontsize=font_size, ncol=3, loc=9, shadow=True)

ax2.set_facecolor('lightgray')
ax2.bar(x-width, list_lifetime_era5, width, color='black', alpha=0.75, edgecolor='black', linewidth=1.25, label='ERA5')
ax2.bar(x, list_lifetime_regcm5, width, color='blue', alpha=0.75, edgecolor='black', linewidth=1.25, label='RegCM5')
ax2.bar(x+width, list_lifetime_wrf415, width, color='red', alpha=0.75, edgecolor='black', linewidth=1.25, label='WRF415')
ax2.set_title('(b)', loc='left', fontsize=font_size, fontweight='bold')
ax2.set_ylabel('Lifetime of cyclones (days)', fontsize=font_size, fontweight='bold')
ax2.set_xlabel('Months', fontsize=font_size, fontweight='bold')
ax2.set_xticks(x)
ax2.set_xticklabels(labels, fontsize=font_size)
ax2.set_ylim(0, 10)
ax2.set_yticks(np.arange(0, 11, 1))
ax2.grid(axis='y', c='k', ls='--', alpha=0.5)

# Path out to save figure
path_out = '{0}/figs/cyclone'.format(path)
name_out = 'pyplt_graph_annual_cycle_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()




			
