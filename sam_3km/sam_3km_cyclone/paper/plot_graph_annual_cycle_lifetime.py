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
		header_list.append(header)

	return header_list


def list_cyclone(header_list):
	
	jan_dt, feb_dt, mar_dt, apr_dt, may_dt, jun_dt, jul_dt, aug_dt, sep_dt, oct_dt, nov_dt, dec_dt = [], [], [], [], [], [], [], [], [], [], [], []
	for ii  in header_list:
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

	cyclone_duration = jan_dur, feb_dur, mar_dur, apr_dur, may_dur, jun_dur, jul_dur, aug_dur, sep_dur, oct_dur, nov_dur, dec_dur
	
	return cyclone_duration
    

# Import model and obs dataset
ec_era5_2018 = open_dat_file('ERA5', '2018')
ec_era5_2019 = open_dat_file('ERA5', '2019')
ec_era5_2020 = open_dat_file('ERA5', '2020')
ec_era5_2021 = open_dat_file('ERA5', '2021')
header_list_era5 = ec_era5_2018 + ec_era5_2019 + ec_era5_2020 + ec_era5_2021
list_cyc_era5 = list_cyclone(header_list_era5)

ec_regcm5_2018 = open_dat_file('RegCM5', '2018')
ec_regcm5_2019 = open_dat_file('RegCM5', '2019')
ec_regcm5_2020 = open_dat_file('RegCM5', '2020')
ec_regcm5_2021 = open_dat_file('RegCM5', '2021')
header_list_regcm5 = ec_regcm5_2018 + ec_regcm5_2019 + ec_regcm5_2020 + ec_regcm5_2021
list_cyc_regcm5 = list_cyclone(header_list_regcm5)

ec_wrf415_2018 = open_dat_file('WRF415', '2018')
ec_wrf415_2019 = open_dat_file('WRF415', '2019')
ec_wrf415_2020 = open_dat_file('WRF415', '2020')
ec_wrf415_2021 = open_dat_file('WRF415', '2021')
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
ax.set_ylabel('Duration of cyclones (2018-2021)', fontsize=font_size, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=font_size)
ax.legend(fontsize=font_size, ncol=2, loc=2, shadow=True)
ax.spines[['right', 'top']].set_visible(False)

# Path out to save figure
path_out = '/marconi/home/userexternal/mdasilva/user/mdasilva/SAM-3km/figs/cyclone'
name_out = 'pyplt_annual_cycle_duration_cyclones_ERA5_RegCM5_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()




			
