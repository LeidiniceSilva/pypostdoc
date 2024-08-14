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
                if rows:  # If we have rows, append them to data
                    data.append((header, rows))
                    rows = []  # Reset rows
                header = line
            else:  # Otherwise, it's a data row
                rows.append(line)

        # Append the last header and rows to data
        if header and rows:
            data.append((header, rows))
	
    return data


def open_dat_file(dataset, yr):

	data = read_dat_file('{0}/{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
	
	header_list = []
	for i, (header, rows) in enumerate(data):
		header_list.append(header[1])
	
	return header_list


def list_cyclone(header_list):
	
	jan_dt, feb_dt, mar_dt, apr_dt, may_dt, jun_dt, jul_dt, aug_dt, sep_dt, oct_dt, nov_dt, dec_dt = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
	for ii  in header_list:
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
	
	cyclone_numbers = jan_dt, feb_dt, mar_dt, apr_dt, may_dt, jun_dt, jul_dt, aug_dt, sep_dt, oct_dt, nov_dt, dec_dt
	
	return cyclone_numbers
    

# Import model and obs dataset
ec_era5_2018 = open_dat_file('ECyclone_ERA5', '2018')
ec_era5_2019 = open_dat_file('ECyclone_ERA5', '2019')
ec_era5_2020 = open_dat_file('ECyclone_ERA5', '2020')
ec_era5_2021 = open_dat_file('ECyclone_ERA5', '2021')
header_list_era5 = ec_era5_2018 + ec_era5_2019 + ec_era5_2020 + ec_era5_2021
list_cyc_era5 = list_cyclone(header_list_era5)

ec_regcm5_2018 = open_dat_file('ECyclone_RegCM5', '2018')
ec_regcm5_2019 = open_dat_file('ECyclone_RegCM5', '2019')
ec_regcm5_2020 = open_dat_file('ECyclone_RegCM5', '2020')
ec_regcm5_2021 = open_dat_file('ECyclone_RegCM5', '2021')
header_list_regcm5 = ec_regcm5_2018 + ec_regcm5_2019 + ec_regcm5_2020 + ec_regcm5_2021
list_cyc_regcm5 = list_cyclone(header_list_regcm5)

# Plot figure
fig, ax = plt.subplots()

labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars
font_size = 10

rects1 = ax.bar(x - width/2, list_cyc_era5, width, color='red', alpha=0.75, edgecolor='black', linewidth=1.5, label='ERA5')
rects2 = ax.bar(x + width/2, list_cyc_regcm5, width, color='blue', alpha=0.75, edgecolor='black', linewidth=1.5, label='RegCM5')
ax.set_xlabel('Months', fontsize=font_size, fontweight='bold')
ax.set_ylabel('Number of cyclones (2018-2021)', fontsize=font_size, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=font_size)
ax.legend(fontsize=font_size, ncol=1, loc=2, shadow=True)
ax.spines[['right', 'top']].set_visible(False)

# Path out to save figure
path_out = '/marconi/home/userexternal/mdasilva/user/mdasilva/SAM-3km/figs/cyclone'
name_out = 'pyplt_annual_cycle_number_cyclones_ERA5_RegCM5_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()




			
