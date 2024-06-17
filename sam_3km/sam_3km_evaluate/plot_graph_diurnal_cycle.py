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
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii

var = 'pr'
domain = 'SAM-3km'
idt, fdt = '2018', '2021'
dt = '{0}-{1}'.format(idt, fdt)

path = '/marconi/home/userexternal/mdasilva'
	
	
def import_situ_i():
	
	mean_i, mean_ii, mean_iii = [], [], []

	skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,7793,100,105,106,107,112,117,124,135,137,139,
	149,152,155,158,168,174,177,183,186,199,204,210,212,224,226,239,240,248,249,253,254,276,277,280,293,298,
	303,305,306,308,319,334,335,341,343,359,362,364,384,393,396,398,399,400,402,413,416,417,422,423,426,427,
	443,444,446,451,453,457,458,467,474,479,483,488,489,490,495,505,509,513,514,516,529,534,544,559,566]
	
	for station in range(1, 567):
		print(station, inmet[station][1])
		
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		yy=inmet[station][2]
		xx=inmet[station][3]

		arq_i  = xr.open_dataset('{0}/OBS/WS-SA/INMET/nc/hourly/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[station][0]))
		data_i = arq_i['pre']
		time_i = data_i.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_i  = time_i.groupby('time.hour').mean('time')
		mean_i.append(var_i.values)
		
		arq_ii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'tp_{0}_ERA5_diurnal_cycle_{1}_lonlat.nc'.format(domain, dt))
		data_ii = arq_ii['tp']
		data_ii = data_ii.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		time_ii = data_ii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_ii  = time_ii.values
		mean_ii.append(var_ii)
		
		arq_iii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'tp_{0}_ERA5_diurnal_cycle_{1}_lonlat.nc'.format(domain, dt))
		data_iii = arq_iii['tp']
		data_iii = data_iii.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		time_iii = data_iii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_iii  = time_iii.values
		mean_iii.append(var_iii)
			
	return mean_i, mean_ii, mean_iii
	

def import_situ_ii():
	
	mean_i, mean_ii, mean_iii = [], [], []
	
	for station in range(1, 73):
		print(station, smn_i[station][0])

		yy=smn_i[station][1]
		xx=smn_i[station][2]

		arq_i  = xr.open_dataset('{0}/OBS/WS-SA/SMN/hourly/nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[station][0]))
		data_i = arq_i['pre']
		time_i = data_i.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_i  = time_i.groupby('time.hour').mean('time')
		mean_i.append(var_i.values)
		
		arq_ii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'tp_{0}_ERA5_diurnal_cycle_{1}_lonlat.nc'.format(domain, dt))
		data_ii = arq_ii['tp']
		data_ii = data_ii.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		time_ii = data_ii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_ii  = time_ii.values
		mean_ii.append(var_ii)
		
		arq_iii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'tp_{0}_ERA5_diurnal_cycle_{1}_lonlat.nc'.format(domain, dt))
		data_iii = arq_iii['tp']
		data_iii = data_iii.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		time_iii = data_iii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_iii  = time_iii.values
		mean_iii.append(var_iii)
			
	return mean_i, mean_ii, mean_iii
	

# Import model and obs dataset
clim_i_x, clim_ii_x, clim_iii_x = import_situ_i()			
clim_i_y, clim_ii_y, clim_iii_y = import_situ_ii()			

inmet_smn = clim_i_x 
era5 = clim_ii_x 
regcm5 = clim_iii_x 

inmet_smn_ = clim_i_y
era5_ = clim_ii_y 
regcm5_ = clim_iii_y 

list_hc = [2, 3, 2, 3, 2, 2, 3, 0, 1, 3, 2, 3, 3, 4, 1, 2, 3, 1, 0, 3, 0, 0, 3, 3, 2, 2, 2, 2, 3, 1, 1, 0, 1, 2, 3, 3, 
1, 1, 2, 2, 0, 0, 3, 1, 3, 3, 0, 0, 2, 3, 0, 2, 3, 2, 2, 0, 0, 2, 3, 4, 2, 3, 2, 1, 3, 0, 0, 1, 3, 4, 3, 2, 2, 3, 1, 0, 
4, 3, 0, 0, 3, 3, 1, 0, 3, 0, 0, 3, 2, 1, 2, 3, 0, 4, 1, 0, 3, 3, 3, 0, 3, 0, 3, 1, 1, 2, 2, 2, 3, 3, 3, 1, 0, 2, 1, 4, 
0, 0, 4, 1, 2, 4, 0, 4, 2, 3, 2, 2, 1, 0, 4, 3, 3, 0, 3, 2, 2, 0, 2, 2, 2, 2, 3, 0, 3, 2, 1, 1, 0, 0, 0, 4, 2, 2, 1, 3, 
0, 2, 3, 0, 3, 2, 4, 3, 2, 3, 0, 3, 2, 2, 3, 3, 2, 0, 3, 0, 3, 2, 1, 1, 2, 1, 0, 3, 2, 3, 0, 3, 2, 1, 1, 1, 2, 0, 3, 3, 
1, 3, 2, 1, 3, 3, 0, 0, 3, 0, 4, 3, 2, 2, 0, 1, 0, 0, 0, 3, 2, 3, 2, 2, 2, 2, 0, 0, 0, 0, 2, 3, 3, 2, 4, 0, 0, 3, 2, 0, 
0, 0, 0, 3, 1, 0, 1, 0, 0, 0, 2, 2, 2, 3, 2, 3, 3, 0, 1, 2, 0, 2, 2, 4, 1, 2, 1, 1, 0, 2, 1, 0, 2, 3, 2, 1, 0, 3, 0, 0, 
3, 3, 3, 0, 1, 3, 3, 0, 0, 0, 2, 0, 2, 3, 3, 2, 3, 2, 1, 2, 2, 0]

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
	inmet_i.append(inmet_smn[c_i])
	era5_i.append(era5[c_i])
	regcm_i.append(regcm5[c_i])

for c_ii in count_ii:
	inmet_ii.append(inmet_smn[c_ii])
	era5_ii.append(era5[c_ii])
	regcm_ii.append(regcm5[c_ii])
	
for c_iii in count_iii:
	inmet_iii.append(inmet_smn[c_iii])
	era5_iii.append(era5[c_iii])
	regcm_iii.append(regcm5[c_iii])
	
for c_iv in count_iv:
	inmet_iv.append(inmet_smn[c_iv])
	era5_iv.append(era5[c_iv])
	regcm_iv.append(regcm5[c_iv])
	
for c_v in count_v:
	inmet_v.append(inmet_smn[c_v])
	era5_v.append(era5[c_v])
	regcm_v.append(regcm5[c_v])
	
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

inmet_c_ii_iii_iv = np.nanmean([inmet_c_ii, inmet_c_iii, inmet_c_iv], axis=0)
era5_c_ii_iii_iv  = np.nanmean([era5_c_ii,  era5_c_iii,  era5_c_iv], axis=0)
regcm_c_ii_iii_iv = np.nanmean([regcm_c_ii, regcm_c_iii, regcm_c_iv], axis=0)

inmet_c_v = np.nanmean(inmet_v, axis=0)
era5_c_v  = np.nanmean(era5_v, axis=0)
regcm_c_v = np.nanmean(regcm_v, axis=0)

inmet_c_vi = np.nanmean(inmet_smn_, axis=0)
era5_c_vi  = np.nanmean(era5_, axis=0)
regcm_c_vi = np.nanmean(regcm5_, axis=0)

# Plot figure
fig = plt.figure(figsize=(8, 6))
time = np.arange(0.5, 24 + 0.5)
font_size = 8

ax = fig.add_subplot(2, 2, 1)
plt.plot(time, inmet_c_i, linewidth=1., color='blue',  markersize=2, markerfacecolor='white', marker='^', label='INMET+SMN')
plt.plot(time, era5_c_i,  linewidth=1., color='green', markersize=2, markerfacecolor='white', marker='s', label='ERA5')
plt.plot(time, regcm_c_i, linewidth=1., color='black', markersize=2, markerfacecolor='white', marker='o', label='RegCM5')
plt.ylabel('Precipitation (mm h$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.title('(a) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(0.06, 0.26)
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.yticks(np.arange(0.06, 0.27, 0.01), fontsize=font_size)
plt.grid(linestyle='--')
plt.legend(loc=4, ncol=2, fontsize=font_size, shadow=True)

ax = fig.add_subplot(2, 2, 2)
plt.plot(time, inmet_c_ii_iii_iv, linewidth=1., color='blue',  markersize=2, markerfacecolor='white', marker='^', label='INMET+SMN')
plt.plot(time, era5_c_ii_iii_iv,  linewidth=1., color='green', markersize=2, markerfacecolor='white', marker='s', label='ERA5')
plt.plot(time, regcm_c_ii_iii_iv, linewidth=1., color='black', markersize=2, markerfacecolor='white', marker='o', label='RegCM5')
plt.title('(b) Cluster II-III-IV', loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(0.06, 0.26)
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.yticks(np.arange(0.06, 0.27, 0.01), fontsize=font_size)
plt.grid(linestyle='--')

ax = fig.add_subplot(2, 2, 3)
plt.plot(time, inmet_c_v, linewidth=1., color='blue',  markersize=2, markerfacecolor='white', marker='^', label='INMET+SMN')
plt.plot(time, era5_c_v,  linewidth=1., color='green', markersize=2, markerfacecolor='white', marker='s', label='ERA5')
plt.plot(time, regcm_c_v, linewidth=1., color='black', markersize=2, markerfacecolor='white', marker='o', label='RegCM5')
plt.title('(c) Cluster V', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('Hours', fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation (mm h$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.ylim(0.06, 0.26)
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.yticks(np.arange(0.06, 0.27, 0.01), fontsize=font_size)
plt.grid(linestyle='--')

ax = fig.add_subplot(2, 2, 4)
plt.plot(time, inmet_c_vi, linewidth=1., color='blue',  markersize=2, markerfacecolor='white', marker='^', label='INMET+SMN')
plt.plot(time, era5_c_vi,  linewidth=1., color='green', markersize=2, markerfacecolor='white', marker='s', label='ERA5')
plt.plot(time, regcm_c_vi, linewidth=1., color='black', markersize=2, markerfacecolor='white', marker='o', label='RegCM5')
plt.title('(d) Cluster VI', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('Hours', fontsize=font_size, fontweight='bold')
plt.ylim(0.06, 0.26)
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.yticks(np.arange(0.06, 0.27, 0.01), fontsize=font_size)
plt.grid(linestyle='--')

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_diurnal_cycle_{0}_{1}_RegCM5_2018-2021.png'.format(var, domain)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
