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

from dict_inmet_stations import inmet

var = 'pr'
dt = '2018-2021'
path = '/marconi/home/userexternal/mdasilva'

skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,7793,100,105,106,107,112,117,124,135,137,139,
149,152,155,158,168,174,177,183,186,199,204,210,212,224,226,239,240,248,249,253,254,276,277,280,293,298,
303,305,306,308,319,334,335,341,343,359,362,364,384,393,396,398,399,400,402,413,416,417,422,423,426,427,
443,444,446,451,453,457,458,467,474,479,483,488,489,490,495,505,509,513,514,516,529,534,544,559,566]
	
	
def import_ws(param):
	
	mean = []
	for station in range(1, 567):
		print(station, inmet[station][0])
		
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		arq  = xr.open_dataset('{0}/OBS/BDMET/database/nc/hourly/{1}/'.format(path, param) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(param, inmet[station][0]))
		data = arq[param]
		time = data.sel(time=slice('2018-06-01','2021-05-31'))
		var  = time.groupby('time.hour').mean('time')
		mean.append(var.values)
										
	return mean
	

def import_obs(param, domain, dataset):
	
	mean = []
	for station in range(1, 567):
		print(station, inmet[station][0])
		
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		arq    = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + '{0}_{1}_{2}_diurnal_cycle_{3}_lonlat.nc'.format(param, domain, dataset, dt))
		data   = arq[param]
		latlon = data.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		time   = latlon.sel(time=slice('2018-06-01','2021-05-31'))
		var    = time.values
		mean.append(var)
										
	return mean


def import_sam_3km(param, domain, dataset):
	
	mean = []
	for station in range(1, 567):
		print(station, inmet[station][0])
		
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		arq    = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/'.format(path) + '{0}_{1}_{2}_diurnal_cycle_{3}_lonlat.nc'.format(param, domain, dataset, dt))
		data   = arq[param]
		latlon = data.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		time   = latlon.sel(time=slice('2018-06-01','2021-05-31'))
		var    = time.values
		mean.append(var)
										
	return mean
	
	
# Import model and obs dataset
dict_var = {'pr': ['pre', 'precip', 'sat_gauge_precip', 'tp']}

inmet_ = import_ws(dict_var[var][0])
era5   = import_obs(dict_var[var][3], 'SAM-3km', 'ERA5')
regcm  = import_sam_3km(var, 'SAM-3km', 'RegCM5')
		
list_hc = [2, 3, 2, 3, 2, 2, 3, 0, 1, 3, 2, 3, 3, 4, 1, 2, 3, 1, 0, 3, 0, 0, 3, 3, 2, 2, 2, 2, 3, 1, 1, 0, 
1, 2, 3, 3, 1, 1, 2, 2, 0, 0, 3, 1, 3, 3, 0, 0, 2, 3, 0, 2, 3, 2, 2, 0, 0, 2, 3, 4, 2, 3, 2, 1, 3, 0, 0, 1, 
3, 4, 3, 2, 2, 3, 1, 0, 4, 3, 0, 0, 3, 3, 1, 0, 3, 0, 0, 3, 2, 1, 2, 3, 0, 4, 1, 0, 3, 3, 3, 0, 3, 0, 3, 1, 
1, 2, 2, 2, 3, 3, 3, 1, 0, 2, 1, 4, 0, 0, 4, 1, 2, 4, 0, 4, 2, 3, 2, 2, 1, 0, 4, 3, 3, 0, 3, 2, 2, 0, 2, 2, 
2, 2, 3, 0, 3, 2, 1, 1, 0, 0, 0, 4, 2, 2, 1, 3, 0, 2, 3, 0, 3, 2, 4, 3, 2, 3, 0, 3, 2, 2, 3, 3, 2, 0, 3, 0, 
3, 2, 1, 1, 2, 1, 0, 3, 2, 3, 0, 3, 2, 1, 1, 1, 2, 0, 3, 3, 1, 3, 2, 1, 3, 3, 0, 0, 3, 0, 4, 3, 2, 2, 0, 1, 
0, 0, 0, 3, 2, 3, 2, 2, 2, 2, 0, 0, 0, 0, 2, 3, 3, 2, 4, 0, 0, 3, 2, 0, 0, 0, 0, 3, 1, 0, 1, 0, 0, 0, 2, 2, 
2, 3, 2, 3, 3, 0, 1, 2, 0, 2, 2, 4, 1, 2, 1, 1, 0, 2, 1, 0, 2, 3, 2, 1, 0, 3, 0, 0, 3, 3, 3, 0, 1, 3, 3, 0, 
0, 0, 2, 0, 2, 3, 3, 2, 3, 2, 1, 2, 2, 0]

count_i, count_ii, count_iii, count_iv, count_v = [], [], [], [], []
for count, idx in enumerate(list_hc):
	if idx == 0:
		count_i.append(count)
	if idx == 1:
		count_ii.append(count)
	if idx == 2:
		count_iii.append(count)
	if idx == 3:
		count_iv.append(count)
	if idx == 4:
		count_v.append(count)
		
inmet_i, inmet_ii,  inmet_iii,  inmet_iv,  inmet_v  = [], [], [], [], []
era5_i,  era5_ii,   era5_iii,   era5_iv,   era5_v   = [], [], [], [], []
regcm_i, regcm_ii,  regcm_iii,  regcm_iv,  regcm_v  = [], [], [], [], []

for c_i in count_i:
	inmet_i.append(inmet_[c_i])
	era5_i.append(era5[c_i])
	regcm_i.append(regcm[c_i])

for c_ii in count_ii:
	inmet_ii.append(inmet_[c_ii])
	era5_ii.append(era5[c_ii])
	regcm_ii.append(regcm[c_ii])
	
for c_iii in count_iii:
	inmet_iii.append(inmet_[c_iii])
	era5_iii.append(era5[c_iii])
	regcm_iii.append(regcm[c_iii])
	
for c_iv in count_iv:
	inmet_iv.append(inmet_[c_iv])
	era5_iv.append(era5[c_iv])
	regcm_iv.append(regcm[c_iv])
	
for c_v in count_v:
	inmet_v.append(inmet_[c_v])
	era5_v.append(era5[c_v])
	regcm_v.append(regcm[c_v])
	
inmet_c_i = np.nanmean(inmet_i, axis=0)
era5_c_i  = np.nanmean(era5_i, axis=0)
regcm_c_i = np.nanmean(regcm_i, axis=0)

inmet_c_ii = np.nanmean(inmet_ii, axis=0)
era5_c_ii  = np.nanmean(era5_ii, axis=0)
regcm_c_ii = np.nanmean(regcm_ii, axis=0)

inmet_c_iii = np.nanmean(inmet_iii, axis=0)
era5_c_iii  = np.nanmean(era5_iii, axis=0)
regcm_c_iii = np.nanmean(regcm_iii, axis=0)

inmet_c_iv = np.nanmean(inmet_iv, axis=0)
era5_c_iv  = np.nanmean(era5_iv, axis=0)
regcm_c_iv = np.nanmean(regcm_iv, axis=0)

inmet_c_v = np.nanmean(inmet_v, axis=0)
era5_c_v  = np.nanmean(era5_v, axis=0)
regcm_c_v = np.nanmean(regcm_v, axis=0)

inmet_c_ii_iv = np.nanmean([inmet_c_ii, inmet_c_iii, inmet_c_iv], axis=0)
era5_c_ii_iv  = np.nanmean([era5_c_ii,  era5_c_iii,  era5_c_iv], axis=0)
regcm_c_ii_iv = np.nanmean([regcm_c_ii, regcm_c_iii, regcm_c_iv], axis=0)

# Plot figure
fig = plt.figure(figsize=(6, 9))
time = np.arange(0.5, 24 + 0.5)
font_size = 8

ax = fig.add_subplot(3, 1, 1)
plt.plot(time, inmet_c_i, linewidth=1.5, color='red',   label = 'INMET')
plt.plot(time, era5_c_i,  linewidth=1.5, color='black', label = 'ERA5')
plt.plot(time, regcm_c_i, linewidth=1.5, linestyle='--', markersize=2, marker='o', markerfacecolor='white', color='red', label = 'RegCM5-3km')
plt.title('(a) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(0, 0.3)
plt.yticks(np.arange(0, 0.33, 0.03), fontsize=font_size)
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.legend(loc=1, ncol=3, fontsize=font_size, shadow=True)

ax = fig.add_subplot(3, 1, 2)
plt.plot(time, inmet_c_ii, linewidth=1.5, color='red',   label = 'INMET')
plt.plot(time, era5_c_ii,  linewidth=1.5, color='black', label = 'ERA5')
plt.plot(time, regcm_c_ii, linewidth=1.5, linestyle='--', markersize=2, marker='o', markerfacecolor='white', color='red', label='RegCM5-3km')
plt.title('(b) Cluster II-III-IV', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation (mm h$^-$$^1$)', fontsize=8, fontweight='bold')
plt.ylim(0, 0.4)
plt.yticks(np.arange(0, 0.44, 0.04), fontsize=font_size)
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)

ax = fig.add_subplot(3, 1, 3)
plt.plot(time, inmet_c_v, linewidth=1.5, color='red',   label = 'INMET')
plt.plot(time, era5_c_v,  linewidth=1.5, color='black', label = 'ERA5')
plt.plot(time, regcm_c_v, linewidth=1.5, linestyle='--', markersize=2, marker='o', markerfacecolor='white', color='red', label='RegCM5-3km')
plt.xlabel('Time (hrs)', fontsize=font_size, fontweight='bold')
plt.title('(c) Cluster V', loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(0, 0.2)
plt.yticks(np.arange(0, 0.22, 0.02), fontsize=font_size)
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cp_3km-4km'.format(path)
name_out = 'pyplt_graph_diurnal_cycle_{0}_SAM-3km_RegCM5_{1}.png'.format(var, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
