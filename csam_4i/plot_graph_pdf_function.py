# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot pdf function"

import os
import netCDF4
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet

var = 'pr'
freq = 'hourly'
path = '/marconi/home/userexternal/mdasilva'

if freq == 'hourly':
	dt = '1hr_2018-2021'
	legend = 'Precipitation (mm h$^-$$^1$)'
else:
	dt = 'day_2018-2021'
	legend = 'Precipitation (mm d$^-$$^1$)'
	
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
		
		if freq == 'hourly':
			var = time.values
			mean.append(var)
		else:
			var = time.resample(time='1D').sum()
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

		arq    = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + '{0}_{1}_{2}_{3}_lonlat_new.nc'.format(param, domain, dataset, dt))
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

		arq    = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/'.format(path) + '{0}_{1}_{2}_{3}_lonlat_new.nc'.format(param, domain, dataset, dt))
		data   = arq[param]
		latlon = data.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		time   = latlon.sel(time=slice('2018-06-01','2021-05-31'))
		var    = time.values
		mean.append(var)
										
	return mean
	
	
# Import model and obs dataset
dict_var = {'pr': ['pre', 'precip', 'sat_gauge_precip', 'tp']}

inmet_ = import_ws(dict_var[var][0])
era5 = import_obs(dict_var[var][3], 'SAM-3km', 'ERA5')
regcm = import_sam_3km(var, 'SAM-3km', 'RegCM5')
		
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

inmet_i = np.array(inmet_i)
inmet_c_i = inmet_i.flatten()

era5_i = np.array(era5_i)
era5_c_i = era5_i.flatten()

regcm_i = np.array(regcm_i)
regcm_c_i = regcm_i.flatten()

inmet_ii = np.array(inmet_ii)
inmet_c_ii = inmet_ii.flatten()

era5_ii = np.array(era5_ii)
era5_c_ii = era5_ii.flatten()

regcm_ii = np.array(regcm_ii)
regcm_c_ii = regcm_ii.flatten()

inmet_iii = np.array(inmet_iii)
inmet_c_iii = inmet_iii.flatten()

era5_iii = np.array(era5_iii)
era5_c_iii = era5_iii.flatten()

regcm_iii = np.array(regcm_iii)
regcm_c_iii = regcm_iii.flatten()

inmet_iv = np.array(inmet_iv)
inmet_c_iv = inmet_iv.flatten()

era5_iv = np.array(era5_iv)
era5_c_iv = era5_iv.flatten()

regcm_iv = np.array(regcm_iv)
regcm_c_iv = regcm_iv.flatten()

inmet_v = np.array(inmet_v)
inmet_c_v = inmet_v.flatten()

era5_v = np.array(era5_v)
era5_c_v = era5_v.flatten()

regcm_v = np.array(regcm_v)
regcm_c_v = regcm_v.flatten()

# Round values to each cluster
round_inmet_c_i = np.round(inmet_c_i,0)
round_era5_c_i  = np.round(era5_c_i,0)
round_regcm_c_i = np.round(regcm_c_i,0)

round_inmet_c_ii_iv = np.round(inmet_c_ii,0)
round_era5_c_ii_iv  = np.round(era5_c_ii,0)
round_regcm_c_ii_iv = np.round(regcm_c_ii,0)

round_inmet_c_v = np.round(inmet_c_v,0)
round_era5_c_v  = np.round(era5_c_v,0)
round_regcm_c_v = np.round(regcm_c_v,0)

# Filter 0 mm/day
filter_inmet_c_i = round_inmet_c_i[round_inmet_c_i > 0.]
filter_era5_c_i  = round_era5_c_i[round_era5_c_i > 0.]
filter_regcm_c_i = round_regcm_c_i[round_regcm_c_i > 0.]

filter_inmet_c_ii_iv = round_inmet_c_ii_iv[round_inmet_c_ii_iv > 0.]
filter_era5_c_ii_iv  = round_era5_c_ii_iv[round_era5_c_ii_iv > 0.]
filter_regcm_c_ii_iv = round_regcm_c_ii_iv[round_regcm_c_ii_iv > 0.]

filter_inmet_c_v = round_inmet_c_v[round_inmet_c_v > 0.]
filter_era5_c_v  = round_era5_c_v[round_era5_c_v > 0.]
filter_regcm_c_v = round_regcm_c_v[round_regcm_c_v > 0.]

# Compute frequency
x_pdf_inmet_c_i, pdf_inmet_c_i = np.unique(filter_inmet_c_i, return_counts=True) 
x_pdf_era5_c_i,  pdf_era5_c_i  = np.unique(filter_era5_c_i, return_counts=True) 
x_pdf_regcm_c_i, pdf_regcm_c_i = np.unique(filter_regcm_c_i, return_counts=True) 

x_pdf_inmet_c_ii_iv, pdf_inmet_c_ii_iv = np.unique(filter_inmet_c_ii_iv, return_counts=True) 
x_pdf_era5_c_ii_iv,  pdf_era5_c_ii_iv  = np.unique(filter_era5_c_ii_iv, return_counts=True) 
x_pdf_regcm_c_ii_iv, pdf_regcm_c_ii_iv = np.unique(filter_regcm_c_ii_iv, return_counts=True) 

x_pdf_inmet_c_v, pdf_inmet_c_v = np.unique(filter_inmet_c_v, return_counts=True) 
x_pdf_era5_c_v,  pdf_era5_c_v  = np.unique(filter_era5_c_v, return_counts=True) 
x_pdf_regcm_c_v, pdf_regcm_c_v = np.unique(filter_regcm_c_v, return_counts=True) 

# Plot figure
fig = plt.figure(figsize=(6, 9))
time = np.arange(0.5, 12 + 0.5)
font_size = 8

ax = fig.add_subplot(3, 1, 1)  
plt.plot(x_pdf_inmet_c_i, pdf_inmet_c_i, marker='.', markersize=4, mfc='red',   mec='red',   linestyle='None', label='INMET')
plt.plot(x_pdf_era5_c_i,  pdf_era5_c_i,  marker='.', markersize=4, mfc='black', mec='black', linestyle='None', label='ERA5')
plt.plot(x_pdf_regcm_c_i, pdf_regcm_c_i, marker='.', markersize=4, mfc='blue', mec='blue',   linestyle='None', label='RegCM5-3km')
plt.title('(a) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.legend(loc=1, ncol=3, fontsize=font_size, shadow=True)

ax = fig.add_subplot(3, 1, 2)  
plt.plot(x_pdf_inmet_c_ii_iv, pdf_inmet_c_ii_iv, marker='.', markersize=4, mfc='red',   mec='red',   linestyle='None', label='INMET')
plt.plot(x_pdf_era5_c_ii_iv,  pdf_era5_c_ii_iv,  marker='.', markersize=4, mfc='black', mec='black', linestyle='None', label='ERA5')
plt.plot(x_pdf_regcm_c_ii_iv, pdf_regcm_c_ii_iv, marker='.', markersize=4, mfc='blue',  mec='blue',  linestyle='None', label='RegCM5-3km')
plt.title('(b) Cluster II-III-IV', loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 1, 3)  
plt.plot(x_pdf_inmet_c_v, pdf_inmet_c_v, marker='.', markersize=4, mfc='red',   mec='red',   linestyle='None', label='INMET')
plt.plot(x_pdf_era5_c_v,  pdf_era5_c_v,  marker='.', markersize=4, mfc='black', mec='black', linestyle='None', label='ERA5')
plt.plot(x_pdf_regcm_c_v, pdf_regcm_c_v, marker='.', markersize=4, mfc='blue',  mec='blue',  linestyle='None', label='RegCM5-3km')
plt.title('(c) Cluster V', loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.xlabel('{0}'.format(legend), fontsize=font_size, fontweight='bold')

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cp_3km-4km'.format(path)
name_out = 'pyplt_graph_pdf_{0}_SAM-3km_RegCM5_{1}.png'.format(var, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
exit()
