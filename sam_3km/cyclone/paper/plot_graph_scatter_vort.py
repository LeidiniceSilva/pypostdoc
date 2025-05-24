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
path = '/leonardo/home/userexternal/mdasilva/leonardo_work'


def process_cyclone_csv_files(dataset):

	vort_values = []
	delta_days_vort_pres = []
	for yr in range(2018, 2021+1):
		file_path = os.path.join('{0}/SAM-3km/postproc/cyclone/{1}/track'.format(path, dataset), 'resultado_{0}.dat'.format(yr))

		if not os.path.exists(file_path):
			print(f'Warming: file not found: {file_path}')
			continue

		cyclones = []
		current_event = []
		with open(file_path, 'r') as f:
			lines = f.readlines()

		for line in lines:
			parts = line.strip().split()

			if len(parts) == 6:
				if current_event:
					cyclones.append(current_event)
				current_event = []
			elif len(parts) == 5:
				try:
					date = datetime.strptime(parts[0], '%Y%m%d%H')
					lat = float(parts[1])
					lon = float(parts[2])
					vort = float(parts[3])
					pres = float(parts[4])
					if lon > -56:
						current_event.append({'date': date, 'vort': vort, 'pres': pres})
				except:
					continue

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
delta_days_era5, vort_era5 = process_cyclone_csv_files('ERA5')
delta_days_regcm5, vort_regcm5 = process_cyclone_csv_files('RegCM5')
delta_days_wrf415, vort_wrf415 = process_cyclone_csv_files('WRF415')

vort_era5_clean, delta_days_era5_clean = clean_nan_values(vort_era5, delta_days_era5)
vort_regcm5_clean, delta_days_regcm5_clean = clean_nan_values(vort_regcm5, delta_days_regcm5)
vort_wrf415_clean, delta_days_wrf415_clean = clean_nan_values(vort_wrf415, delta_days_wrf415)

freq_norm_era5 = compute_normalized_frequency(vort_era5_clean)
freq_norm_regcm5 = compute_normalized_frequency(vort_regcm5_clean)
freq_norm_wrf415 = compute_normalized_frequency(vort_wrf415_clean)

# Plot figure
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))

ax1.scatter(delta_days_era5_clean, vort_era5_clean, c=freq_norm_era5, cmap=cm.viridis, marker='s', s=100, vmin=0, vmax=0.9)
ax1.set_xticks(range(-5, 6))  
ax1.set_xlabel('Days', fontsize=font_size, fontweight='bold')
ax1.set_ylabel('Vorticity', fontsize=font_size, fontweight='bold')
ax1.set_title('(a)', loc='left', fontsize=font_size, fontweight='bold')
ax1.grid(True, linestyle='--', alpha=0.5)
ax1.legend()

ax2.scatter(delta_days_regcm5_clean, vort_regcm5_clean, c=freq_norm_regcm5, cmap=cm.viridis, marker='s', s=100, vmin=0, vmax=0.9)
ax2.set_xticks(range(-5, 6))  
ax2.set_xlabel('Days', fontsize=font_size, fontweight='bold')
ax2.set_title('(b)', loc='left', fontsize=font_size, fontweight='bold')
ax2.grid(True, linestyle='--', alpha=0.5)
ax2.legend()

sc3 = ax3.scatter(delta_days_wrf415_clean, vort_wrf415_clean, c=freq_norm_wrf415, cmap=cm.viridis, marker='s', s=100, vmin=0, vmax=0.9)
ax3.set_xticks(range(-5, 6)) 
ax3.set_xlabel('Days', fontsize=font_size)
ax3.set_title('(c)', loc='left', fontsize=font_size, fontweight='bold')
ax3.grid(True, linestyle='--', alpha=0.5)
ax3.legend()

cbar = fig.colorbar(sc3, ax=[ax1, ax2, ax3])

# Path out to save figure
path_out = '{0}/SAM-3km/figs/cyclone'.format(path)
name_out = 'pyplt_graph_scatter_vort_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

