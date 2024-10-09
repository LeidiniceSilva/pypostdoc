# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot diurnal cycle"

import os
import netCDF4
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.pyplot as plt

var = 'pr'
domain = 'CSAM-3'
idt, fdt = '2000', '2009'
dt = '{0}-{1}'.format(idt, fdt)

dict_ws = {
	   0: [-68.1193, -16.4897, 'La Paz'],
	   1: [-47.9000, -16.0000, 'Brazilia'],
	   2: [-70.6693, -33.4489, 'Santiago'],
	   3: [-57.5759, -25.2637, 'Asuncion'],
	   4: [-58.4004, -34.6051, 'Buenos Aires'],
	   5: [-56.0000, -34.0000, 'Montevideo']
	   }

path = '/marconi/home/userexternal/mdasilva'

		
def import_obs(param, dataset):
	
	ts = []
	for i in range(0, 6):
		yy=dict_ws[i][1]
		xx=dict_ws[i][0]

		arq  = xr.open_dataset('{0}/user/mdasilva/CORDEX/post_evaluate/obs/{1}_{2}_diurnal_cycle_2000-2009_lonlat.nc'.format(path, param, dataset))
		data = arq[param]
		var  = data.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		ts.append(var.values)

	return ts


def import_rcm(param, dataset):

	ts = []
	for i in range(0, 6):
		yy=dict_ws[i][1]
		xx=dict_ws[i][0]
		
		arq  = xr.open_dataset('{0}/user/mdasilva/CORDEX/post_evaluate/rcm/{1}_{2}_diurnal_cycle_2000-2005_lonlat.nc'.format(path, param, dataset))
		data = arq[param]
		var  = data.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		ts.append(var.values)
		
	return ts


# Import model and obs dataset
dict_var = {'pr': ['cmorph', 'pr']}

ts_cmorph = import_obs(dict_var[var][0], 'CSAM-3_CMORPH')
ts_era5   = import_obs(dict_var[var][1], 'CSAM-3_ERA5')
ts_rcm3   = import_rcm(var, 'CSAM-3_RegCM5')

# Plot figure
fig = plt.figure(figsize=(14, 10))
time = np.arange(0.5, 24 + 0.5)
font_size = 8

ax = fig.add_subplot(3, 3, 1)
plt.plot(time, ts_rcm3[0],   linewidth=1., color='black', markersize=2, marker='o', markerfacecolor='white', label='CPM3')
plt.plot(time, ts_cmorph[0], linewidth=1., color='red',   markersize=2, marker='o', markerfacecolor='white', label='CMORPH')
plt.plot(time, ts_era5[0],   linewidth=1., color='blue',  markersize=2, marker='o', markerfacecolor='white', label='ERA5')
plt.title('(a) {0}'.format(dict_ws[0][2]), loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(0.0, 0.50)
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.yticks(np.arange(0.0, 0.53, 0.03), fontsize=font_size)
plt.legend(loc=9, ncol=1, fontsize=font_size, shadow=True)
plt.grid(linestyle='--')

ax = fig.add_subplot(3, 3, 3)
plt.plot(time, ts_rcm3[1],   linewidth=1., color='black', markersize=2, marker='o', markerfacecolor='white', label='CPM3')
plt.plot(time, ts_cmorph[1], linewidth=1., color='red',   markersize=2, marker='o', markerfacecolor='white', label='CMORPH')
plt.plot(time, ts_era5[1],   linewidth=1., color='blue',  markersize=2, marker='o', markerfacecolor='white', label='ERA5')
plt.title('(b) {0}'.format(dict_ws[1][2]), loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(0.0, 0.50)
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.yticks(np.arange(0.0, 0.53, 0.03), fontsize=font_size)
plt.grid(linestyle='--')

ax = fig.add_subplot(3, 3, 4)
plt.plot(time, ts_rcm3[2],   linewidth=1., color='black', markersize=2, marker='o', markerfacecolor='white', label='CPM3')
plt.plot(time, ts_cmorph[2], linewidth=1., color='red',   markersize=2, marker='o', markerfacecolor='white', label='CMORPH')
plt.plot(time, ts_era5[2],   linewidth=1., color='blue',  markersize=2, marker='o', markerfacecolor='white', label='ERA5')
plt.title('(c) {0}'.format(dict_ws[2][2]), loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation (mm h$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.ylim(0.0, 0.50)
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.yticks(np.arange(0.0, 0.53, 0.03), fontsize=font_size)
plt.grid(linestyle='--')

ax = fig.add_subplot(3, 3, 6)
plt.plot(time, ts_rcm3[3],   linewidth=1., color='black', markersize=2, marker='o', markerfacecolor='white', label='CPM3')
plt.plot(time, ts_cmorph[3], linewidth=1., color='red',   markersize=2, marker='o', markerfacecolor='white', label='CMORPH')
plt.plot(time, ts_era5[3],   linewidth=1., color='blue',  markersize=2, marker='o', markerfacecolor='white', label='ERA5')
plt.title('(d) {0}'.format(dict_ws[3][2]), loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation (mm h$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.ylim(0.0, 0.50)
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.yticks(np.arange(0.0, 0.53, 0.03), fontsize=font_size)
plt.grid(linestyle='--')

ax = fig.add_subplot(3, 3, 7)
plt.plot(time, ts_rcm3[4],   linewidth=1., color='black', markersize=2, marker='o', markerfacecolor='white', label='CPM3')
plt.plot(time, ts_cmorph[4], linewidth=1., color='red',   markersize=2, marker='o', markerfacecolor='white', label='CMORPH')
plt.plot(time, ts_era5[4],   linewidth=1., color='blue',  markersize=2, marker='o', markerfacecolor='white', label='ERA5')
plt.title('(e) {0}'.format(dict_ws[4][2]), loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('Hours', fontsize=font_size, fontweight='bold')
plt.ylim(0.0, 0.50)
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.yticks(np.arange(0.0, 0.53, 0.03), fontsize=font_size)
plt.grid(linestyle='--')

ax = fig.add_subplot(3, 3, 9)
plt.plot(time, ts_rcm3[5],   linewidth=1., color='black', markersize=2, marker='o', markerfacecolor='white', label='CPM3')
plt.plot(time, ts_cmorph[5], linewidth=1., color='red',   markersize=2, marker='o', markerfacecolor='white', label='CMORPH')
plt.plot(time, ts_era5[5],   linewidth=1., color='blue',  markersize=2, marker='o', markerfacecolor='white', label='ERA5')
plt.title('(f) {0}'.format(dict_ws[5][2]), loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('Hours', fontsize=font_size, fontweight='bold')
plt.ylim(0.0, 0.50)
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.yticks(np.arange(0.0, 0.53, 0.03), fontsize=font_size)
plt.grid(linestyle='--')

# Path out to save figure
path_out = '{0}/user/mdasilva/CORDEX/figs'.format(path)
name_out = 'pyplt_graph_diurnal_cycle_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
