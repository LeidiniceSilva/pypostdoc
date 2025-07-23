# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot wind histogram"

import os
import netCDF4
import datetime
import numpy as np
import xarray as xr
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from scipy import signal, misc
from datetime import datetime

font_size = 10

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


def process_cyclone_csv_files(dataset, yr_init, yr_end):

	vort_values = []
	delta_days_vort_pres = []
	for yr in range(yr_init, yr_end+1):

		data = read_dat_file('{0}/postproc/cyclone/{1}/track/resultado_{2}.dat'.format(path, dataset, yr))

		rows_list = []
		current_event = []
		for i, (header, rows) in enumerate(data):
			if (rows[0][2] < str(-56)):
				rows_list.append(rows[0])
		
		for j  in rows_list:
			dt = datetime.strptime(j[0], '%Y%m%d%H')
			vo = float(j[3])
			pc = float(j[4])
			current_event.append({'date': dt, 'vort': vo, 'pres': pc})
		print(current_event)

		cyclones = []
		if current_event:
			cyclones.append(current_event)


		for event in cyclones:
			if len(event) < 2:
				continue  

			vorts = [d['vort'] for d in event]
			press = [d['pres'] for d in event]
			dates = [d['date'] for d in event]

			try:
				max_vort_idx = vorts.index(min(vorts))  
				min_pres_idx = press.index(min(press)) 

				delta_days = (dates[max_vort_idx] - dates[min_pres_idx]).days
				delta_days_vort_pres.append(delta_days)
				vort_values.append(vorts[max_vort_idx])
			except:
				continue

	return delta_days_vort_pres, vort_values


def clean_nan_values(vort_dataset, delta_days_dataset):

	vort_dataset_clean, delta_days_dataset_clean = zip(
	*[(v, d) for v, d in zip(vort_dataset, delta_days_dataset) if v != -99])

	vort_dataset_clean = list(vort_dataset_clean)
	delta_days_dataset_clean = list(delta_days_dataset_clean)
	
	return vort_dataset_clean, delta_days_dataset_clean


def compute_normalized_frequency(vort_values):

	vorts = np.array(vort_values)
	vorts_abs = np.abs(vorts)  
	freq_normalized = (vorts_abs - np.min(vorts_abs)) / (np.max(vorts_abs) - np.min(vorts_abs))

	return freq_normalized


# Import cyclone tracking values 
delta_days_era5, vort_era5 = process_cyclone_csv_files('ERA5', 2018, 2021)
delta_days_regcm5, vort_regcm5 = process_cyclone_csv_files('RegCM5', 2018, 2021)
delta_days_wrf415, vort_wrf415 = process_cyclone_csv_files('WRF415', 2018, 2021)

vort_era5_clean, delta_days_era5_clean = clean_nan_values(vort_era5, delta_days_era5)
vort_regcm5_clean, delta_days_regcm5_clean = clean_nan_values(vort_regcm5, delta_days_regcm5)
vort_wrf415_clean, delta_days_wrf415_clean = clean_nan_values(vort_wrf415, delta_days_wrf415)

freq_norm_era5 = compute_normalized_frequency(vort_era5_clean)
freq_norm_regcm5 = compute_normalized_frequency(vort_regcm5_clean)
freq_norm_wrf415 = compute_normalized_frequency(vort_wrf415_clean)

# Plot figure
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))

ax1.set_facecolor('lightgray')
ax1.scatter(delta_days_era5_clean, vort_era5_clean, c=freq_norm_era5, cmap=cm.viridis, marker='s', s=100, vmin=0, vmax=0.9)
ax1.set_xticks(range(-3, 4))  
ax1.set_yticks(range(-7, -1))  
ax1.set_xlabel('Days', fontsize=font_size, fontweight='bold')
ax1.set_ylabel('Relative vorticity (10$^-$$^5$ s$^-$$^1$)', fontsize=font_size, fontweight='bold')
ax1.set_title('(a)', loc='left', fontsize=font_size, fontweight='bold')
ax1.grid(True, c='k', ls='--', alpha=0.5)

ax2.set_facecolor('lightgray')
ax2.scatter(delta_days_regcm5_clean, vort_regcm5_clean, c=freq_norm_regcm5, cmap=cm.viridis, marker='s', s=100, vmin=0, vmax=0.9)
ax2.set_xticks(range(-3, 4))  
ax2.set_yticks(range(-7, -1))  
ax2.set_xlabel('Days', fontsize=font_size, fontweight='bold')
ax2.set_title('(b)', loc='left', fontsize=font_size, fontweight='bold')
ax2.grid(True, c='k', ls='--', alpha=0.5)

ax3.set_facecolor('lightgray')
sc3 = ax3.scatter(delta_days_wrf415_clean, vort_wrf415_clean, c=freq_norm_wrf415, cmap=cm.viridis, marker='s', s=100, vmin=0, vmax=0.9)
ax3.set_xticks(range(-3, 4)) 
ax3.set_yticks(range(-7, -1))  
ax3.set_xlabel('Days', fontsize=font_size, fontweight='bold')
ax3.set_title('(c)', loc='left', fontsize=font_size, fontweight='bold')
ax3.grid(True, c='k', ls='--', alpha=0.5)

cbar = fig.colorbar(sc3, ax=[ax1, ax2, ax3], pad=0.02)
cbar.set_label('Normalized frequency', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/SAM-3km/figs/cyclone'.format(path)
name_out = 'pyplt_graph_scatter_vort-pcent_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

