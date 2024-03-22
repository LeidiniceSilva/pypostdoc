# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jan 02, 2024"
__description__ = "This script plot diurnal cycle"

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap

path = '/marconi/home/userexternal/mdasilva/user/mdasilva/SAM-3km-cyclone'


def import_inmet(param):
	
	mean = []
	for i in range(1, 14):

		arq  = xr.open_dataset('{0}/post/obs/inmet/'.format(path, param) + '{0}_{1}_H_2023-01-01_2023-12-31.nc'.format(param, inmet[i][0]))
		data = arq[param]
		time = data.sel(time=slice('2023-06-01','2023-07-31'))
		var  = time.groupby('time.hour').mean('time')
		var_ = np.bincount(np.arange(len(var))//3, var)
		mean.append(var_)
										
	return mean


def import_era5(param):
	
	mean = []
	for i in range(1, 14):

		arq    = xr.open_dataset('{0}/post/obs/era5/'.format(path, param) + '{0}_SAM-3km-cyclone_ERA5_1hr_dc_2023060100-2023073100_lonlat.nc'.format(param))
		data   = arq[param]
		latlon = data.sel(lat=slice(inmet[i][2]-0.03,inmet[i][2]+0.03),lon=slice(inmet[i][3]-0.03,inmet[i][3]+0.03)).mean(('lat','lon'))
		var    = np.bincount(np.arange(len(latlon))//3, latlon)
		mean.append(var*1000)
										
	return mean
		

def import_persiann(param):
	
	mean = []
	for i in range(1, 14):

		arq  = xr.open_dataset('{0}/post/obs/persiann/'.format(path, param) + '{0}_persiann_1hr_dc_2023060100-2023073100_lonlat.nc'.format(param))
		data   = arq[param]
		latlon = data.sel(lat=slice(inmet[i][2]-0.03,inmet[i][2]+0.03),lon=slice(inmet[i][3]-0.03,inmet[i][3]+0.03)).mean(('lat','lon'))
		var    = np.bincount(np.arange(len(latlon))//3, latlon)
		mean.append(var)
										
	return mean


def import_regcm5(param):
	
	mean = []
	for i in range(1, 14):

		arq  = xr.open_dataset('{0}/post/rcm/'.format(path, param) + '{0}_SAM-3km-cyclone_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5_3hr_dc_2023060100-2023073100_lonlat.nc'.format(param))
		data   = arq[param]
		latlon = data.sel(lat=slice(inmet[i][2]-0.03,inmet[i][2]+0.03),lon=slice(inmet[i][3]-0.03,inmet[i][3]+0.03)).mean(('lat','lon'))
		var    = latlon.values
		mean.append(var*10800)
										
	return mean
	
	
# Import model and obs dataset
var = 'pr'

dict_var = {
'pr': ['pre', 'tp', 'prec'],
'srfWind': ['uv', 'srfWind'],
'psl': ['psl', 'msl']
}

inmet_    = import_inmet(dict_var[var][0])
era5_     = import_era5(dict_var[var][1])
persiann_ = import_persiann(dict_var[var][2])
regcm5_   = import_regcm5(var)

# Plot study area
fig = plt.figure(figsize=(14, 8)) 
time = np.arange(0.5, 8 + 0.5)
font_size = 8

ax = fig.add_subplot(4, 4, 1)
plt.plot(time, inmet_[0],    linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='black', label='INMET')
plt.plot(time, era5_[0],     linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='red',   label='ERA5')
plt.plot(time, persiann_[0], linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='green', label='PERSIANN')
plt.plot(time, regcm5_[0],   linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='blue',  label='RegCM5')
plt.title('(a) {0}'.format(inmet[1][1]), loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation (mm h⁻¹)', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '03', '06', '09', '12', '15', '18', '21'), fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(linestyle='--')
plt.legend(loc=1, ncol=2, fontsize=font_size, shadow=True)

ax = fig.add_subplot(4, 4, 2)
plt.plot(time, inmet_[1],    linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='black', label='INMET')
plt.plot(time, era5_[1],     linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='red',   label='ERA5')
plt.plot(time, persiann_[1], linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='green', label='PERSIANN')
plt.plot(time, regcm5_[1],   linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='blue',  label='RegCM5')
plt.title('(b) {0}'.format(inmet[2][1]), loc='left', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '03', '06', '09', '12', '15', '18', '21'), fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(linestyle='--')

ax = fig.add_subplot(4, 4, 3)
plt.plot(time, inmet_[2],    linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='black', label='INMET')
plt.plot(time, era5_[2],     linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='red',   label='ERA5')
plt.plot(time, persiann_[2], linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='green', label='PERSIANN')
plt.plot(time, regcm5_[2],   linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='blue',  label='RegCM5')
plt.title('(c) {0}'.format(inmet[3][1]), loc='left', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '03', '06', '09', '12', '15', '18', '21'), fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(linestyle='--')

ax = fig.add_subplot(4, 4, 4)
plt.plot(time, inmet_[3],    linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='black', label='INMET')
plt.plot(time, era5_[3],     linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='red',   label='ERA5')
plt.plot(time, persiann_[3], linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='green', label='PERSIANN')
plt.plot(time, regcm5_[3],   linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='blue',  label='RegCM5')
plt.title('(d) {0}'.format(inmet[4][1]), loc='left', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '03', '06', '09', '12', '15', '18', '21'), fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(linestyle='--')

ax = fig.add_subplot(4, 4, 5)
plt.plot(time, inmet_[4],    linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='black', label='INMET')
plt.plot(time, era5_[4],     linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='red',   label='ERA5')
plt.plot(time, persiann_[4], linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='green', label='PERSIANN')
plt.plot(time, regcm5_[4],   linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='blue',  label='RegCM5')
plt.title('(e) {0}'.format(inmet[5][1]), loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation (mm h⁻¹)', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '03', '06', '09', '12', '15', '18', '21'), fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(linestyle='--')

ax = fig.add_subplot(4, 4, 6)
plt.plot(time, inmet_[5],    linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='black', label='INMET')
plt.plot(time, era5_[5],     linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='red',   label='ERA5')
plt.plot(time, persiann_[5], linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='green', label='PERSIANN')
plt.plot(time, regcm5_[5],   linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='blue',  label='RegCM5')
plt.title('(f) {0}'.format(inmet[6][1]), loc='left', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '03', '06', '09', '12', '15', '18', '21'), fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(linestyle='--')

ax = fig.add_subplot(4, 4, 7)
plt.plot(time, inmet_[6],    linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='black', label='INMET')
plt.plot(time, era5_[6],     linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='red',   label='ERA5')
plt.plot(time, persiann_[6], linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='green', label='PERSIANN')
plt.plot(time, regcm5_[6],   linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='blue',  label='RegCM5')
plt.title('(g) {0}'.format(inmet[7][1]), loc='left', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '03', '06', '09', '12', '15', '18', '21'), fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(linestyle='--')

ax = fig.add_subplot(4, 4, 8)
plt.plot(time, inmet_[7],    linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='black', label='INMET')
plt.plot(time, era5_[7],     linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='red',   label='ERA5')
plt.plot(time, persiann_[7], linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='green', label='PERSIANN')
plt.plot(time, regcm5_[7],   linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='blue',  label='RegCM5')
plt.title('(h) {0}'.format(inmet[8][1]), loc='left', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '03', '06', '09', '12', '15', '18', '21'), fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(linestyle='--')

ax = fig.add_subplot(4, 4, 9)
plt.plot(time, inmet_[8],    linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='black', label='INMET')
plt.plot(time, era5_[8],     linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='red',   label='ERA5')
plt.plot(time, persiann_[8], linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='green', label='PERSIANN')
plt.plot(time, regcm5_[8],   linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='blue',  label='RegCM5')
plt.title('(i) {0}'.format(inmet[9][1]), loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation (mm h⁻¹)', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '03', '06', '09', '12', '15', '18', '21'), fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(linestyle='--')

ax = fig.add_subplot(4, 4, 10)
plt.plot(time, inmet_[9],    linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='black', label='INMET')
plt.plot(time, era5_[9],     linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='red',   label='ERA5')
plt.plot(time, persiann_[9], linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='green', label='PERSIANN')
plt.plot(time, regcm5_[9],   linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='blue',  label='RegCM5')
plt.title('(j) {0}'.format(inmet[10][1]), loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('Hours', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '03', '06', '09', '12', '15', '18', '21'), fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.grid(linestyle='--')

ax = fig.add_subplot(4, 4, 11)
plt.plot(time, inmet_[10],    linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='black', label='INMET')
plt.plot(time, era5_[10],     linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='red',   label='ERA5')
plt.plot(time, persiann_[10], linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='green', label='PERSIANN')
plt.plot(time, regcm5_[10],   linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='blue',  label='RegCM5')
plt.title('(k) {0}'.format(inmet[11][1]), loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('Hours', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '03', '06', '09', '12', '15', '18', '21'), fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.grid(linestyle='--')

ax = fig.add_subplot(4, 4, 12)
plt.plot(time, inmet_[11],    linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='black', label='INMET')
plt.plot(time, era5_[11],     linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='red',   label='ERA5')
plt.plot(time, persiann_[11], linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='green', label='PERSIANN')
plt.plot(time, regcm5_[11],   linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='blue',  label='RegCM5')
plt.title('(l) {0}'.format(inmet[12][1]), loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('Hours', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '03', '06', '09', '12', '15', '18', '21'), fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.grid(linestyle='--')

ax = fig.add_subplot(4, 4, 13)
plt.plot(time, inmet_[12],    linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='black', label='INMET')
plt.plot(time, era5_[12],     linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='red',   label='ERA5')
plt.plot(time, persiann_[12], linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='green', label='PERSIANN')
plt.plot(time, regcm5_[12],   linewidth=1., markersize=2, marker='o', markerfacecolor='white', color='blue',  label='RegCM5')
plt.title('(m) {0}'.format(inmet[13][1]), loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation (mm h⁻¹)', fontsize=font_size, fontweight='bold')
plt.xlabel('Hours', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '03', '06', '09', '12', '15', '18', '21'), fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.grid(linestyle='--')

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_graph_diurnal_cycle_{0}.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
