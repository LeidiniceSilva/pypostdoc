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

path = '/marconi/home/userexternal/mdasilva/user/mdasilva/SAM-3km'


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
			if rows:  # If we have rows, append them to data
				data.append((header, rows))
				rows = []  # Reset rows
			header = line
		else:
			rows.append(line)
	
	# Append the last header and rows to data
	if header and rows:
		data.append((header, rows))
	
	return data


def open_dat_file(dataset, yr):

	data = read_dat_file('{0}/post_cyclone/ECyclone_v2/{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
	
	header_list_i = []
	header_list_ii = []
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
	
	jan_dur = np.nanmean(jan_dt)
	feb_dur = np.nanmean(feb_dt)
	mar_dur = np.nanmean(mar_dt)
	apr_dur = np.nanmean(apr_dt)
	may_dur = np.nanmean(may_dt)
	jun_dur = np.nanmean(jun_dt)
	jul_dur = np.nanmean(jul_dt)
	aug_dur = np.nanmean(aug_dt)
	sep_dur = np.nanmean(sep_dt)
	oct_dur = np.nanmean(oct_dt)
	nov_dur = np.nanmean(nov_dt)
	dec_dur = np.nanmean(dec_dt)

	list_lifetime = jan_dur, feb_dur, mar_dur, apr_dur, may_dur, jun_dur, jul_dur, aug_dur, sep_dur, oct_dur, nov_dur, dec_dur
	
	return list_lifetime
	
	
# Import model and obs dataset
era5_2018_i, era5_2018_ii = open_dat_file('ERA5', '2018')
era5_2019_i, era5_2019_ii = open_dat_file('ERA5', '2019')
era5_2020_i, era5_2020_ii = open_dat_file('ERA5', '2020')
era5_2021_i, era5_2021_ii = open_dat_file('ERA5', '2021')
list_era5_i = era5_2018_i + era5_2019_i + era5_2020_i + era5_2021_i
list_era5_ii = era5_2018_ii + era5_2019_ii + era5_2020_ii + era5_2021_ii

regcm5_2018_i, regcm5_2018_ii = open_dat_file('RegCM5', '2018')
regcm5_2019_i, regcm5_2019_ii = open_dat_file('RegCM5', '2019')
regcm5_2020_i, regcm5_2020_ii = open_dat_file('RegCM5', '2020')
regcm5_2021_i, regcm5_2021_ii = open_dat_file('RegCM5', '2021')
list_regcm5_i = regcm5_2018_i + regcm5_2019_i + regcm5_2020_i + regcm5_2021_i
list_regcm5_ii = regcm5_2018_ii + regcm5_2019_ii + regcm5_2020_ii + regcm5_2021_ii

wrf415_2018_i, wrf415_2018_ii = open_dat_file('WRF415', '2018')
wrf415_2019_i, wrf415_2019_ii = open_dat_file('WRF415', '2019')
wrf415_2020_i, wrf415_2020_ii = open_dat_file('WRF415', '2020')
wrf415_2021_i, wrf415_2021_ii = open_dat_file('WRF415', '2021')
list_wrf415_i = wrf415_2018_i + wrf415_2019_i + wrf415_2020_i + wrf415_2021_i
list_wrf415_ii = wrf415_2018_ii + wrf415_2019_ii + wrf415_2020_ii + wrf415_2021_ii

list_number_era5_ = cyclone_number(list_era5_i)
list_number_regcm5_ = cyclone_number(list_regcm5_i)
list_number_wrf415_ = cyclone_number(list_wrf415_i)

list_lifetime_era5_ = cyclone_lifetime(list_era5_ii)
list_lifetime_regcm5_ = cyclone_lifetime(list_regcm5_ii)
list_lifetime_wrf415_ = cyclone_lifetime(list_wrf415_ii)

# Plot figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
x = np.arange(len(labels))  
font_size = 10
width = 0.25 

ax1.bar(x-width, list_number_era5_,   width, color='black', alpha=0.75, edgecolor='black', linewidth=1., label='ERA5')
ax1.bar(x,       list_number_regcm5_, width, color='blue',  alpha=0.75, edgecolor='black', linewidth=1., label='RegCM5')
ax1.bar(x+width, list_number_wrf415_, width, color='red',   alpha=0.75, edgecolor='black', linewidth=1., label='WRF415')
ax1.set_title('(a)', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_xlabel('Months', fontsize=font_size, fontweight='bold')
ax1.set_ylabel('Number of cyclones', fontsize=font_size, fontweight='bold')
ax1.set_xticks(x)
ax1.set_xticklabels(labels, fontsize=font_size)
ax1.set_ylim(0, 10)
ax1.grid(axis='y', c='k', ls='--', alpha=0.3)
ax1.axvspan(0, 3, color='k', alpha=0.1)
ax1.axvspan(8, 11, color='k', alpha=0.1)
ax1.legend(fontsize=font_size, ncol=1, loc=9, shadow=True)

ax2.bar(x-width, list_lifetime_era5_,   width, color='black', alpha=0.75, edgecolor='black', linewidth=1., label='ERA5')
ax2.bar(x,       list_lifetime_regcm5_, width, color='blue',  alpha=0.75, edgecolor='black', linewidth=1., label='RegCM5')
ax2.bar(x+width, list_lifetime_wrf415_, width, color='red',   alpha=0.75, edgecolor='black', linewidth=1., label='WRF415')
ax2.set_title('(b)', loc='left', fontsize=font_size, fontweight='bold')
ax2.set_xlabel('Months', fontsize=font_size, fontweight='bold')
ax2.set_ylabel('LIfetime of cyclones', fontsize=font_size, fontweight='bold')
ax2.set_xticks(x)
ax2.set_xticklabels(labels, fontsize=font_size)
ax2.set_ylim(0, 10)
ax2.grid(axis='y', c='k', ls='--', alpha=0.3)
ax2.axvspan(0, 3, color='k', alpha=0.1)
ax2.axvspan(8, 11, color='k', alpha=0.1)

# Path out to save figure
path_out = '{0}/figs/cyclone/paper'.format(path)
name_out = 'pyplt_graph_annual_cycle_count_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()




			
