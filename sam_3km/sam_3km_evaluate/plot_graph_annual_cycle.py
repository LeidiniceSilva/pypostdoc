# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot annual cycle"

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
	
	mean_i, mean_ii, mean_iii, mean_iv, mean_v = [], [], [], [], []

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
		var_i  = time_i.groupby('time.month').mean('time')
		mean_i.append(var_i.values*24)

		arq_ii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'pre_{0}_CRU_mon_{1}_lonlat.nc'.format(domain, dt))
		data_ii = arq_ii['pre']
		data_ii = data_ii.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		time_ii = data_ii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_ii  = time_ii.groupby('time.month').mean('time')
		mean_ii.append(var_ii.values)

		arq_iii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'precip_{0}_CPC_mon_{1}_lonlat.nc'.format(domain, dt))
		data_iii = arq_iii['precip']
		data_iii = data_iii.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		time_iii = data_iii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_iii  = time_iii.groupby('time.month').mean('time')
		mean_iii.append(var_iii.values)
		
		arq_iv  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'tp_{0}_ERA5_mon_{1}_lonlat.nc'.format(domain, dt))
		data_iv = arq_iv['tp']
		data_iv = data_iv.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		time_iv = data_iv.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_iv  = time_iv.groupby('time.month').mean('time')
		mean_iv.append(var_iv.values)
		
		arq_v  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/'.format(path) + 'pr_{0}_RegCM5_mon_{1}_lonlat.nc'.format(domain, dt))
		data_v = arq_v['pr']
		data_v = data_v.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		time_v = data_v.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_v  = time_v.groupby('time.month').mean('time')
		mean_v.append(var_v.values)
			
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v
	

def import_situ_ii():
	
	mean_i, mean_ii, mean_iii, mean_iv, mean_v = [], [], [], [], []
	
	for station in range(1, 73):
		print(station, smn_i[station][0])

		yy=smn_i[station][1]
		xx=smn_i[station][2]

		arq_i  = xr.open_dataset('{0}/OBS/WS-SA/SMN/hourly/nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[station][0]))
		data_i = arq_i['pre']
		time_i = data_i.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_i  = time_i.groupby('time.month').mean('time')
		mean_i.append(var_i.values*24)

		arq_ii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'pre_{0}_CRU_mon_{1}_lonlat.nc'.format(domain, dt))
		data_ii = arq_ii['pre']
		data_ii = data_ii.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		time_ii = data_ii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_ii  = time_ii.groupby('time.month').mean('time')
		mean_ii.append(var_ii.values)

		arq_iii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'precip_{0}_CPC_mon_{1}_lonlat.nc'.format(domain, dt))
		data_iii = arq_iii['precip']
		data_iii = data_iii.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		time_iii = data_iii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_iii  = time_iii.groupby('time.month').mean('time')
		mean_iii.append(var_iii.values)
		
		arq_iv  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'tp_{0}_ERA5_mon_{1}_lonlat.nc'.format(domain, dt))
		data_iv = arq_iv['tp']
		data_iv = data_iv.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		time_iv = data_iv.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_iv  = time_iv.groupby('time.month').mean('time')
		mean_iv.append(var_iv.values)
		
		arq_v  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/'.format(path) + 'pr_{0}_RegCM5_mon_{1}_lonlat.nc'.format(domain, dt))
		data_v = arq_v['pr']
		data_v = data_v.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		time_v = data_v.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_v  = time_v.groupby('time.month').mean('time')
		mean_v.append(var_v.values)
			
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v
	

def import_situ_iii():
	
	mean_i, mean_ii, mean_iii, mean_iv, mean_v = [], [], [], [], []
	
	skip_list = [86,87,88,89,90,91,95,96,97,98,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,
	117,118,119,120,121,122,123,124,125,126,127,128]	
		
	for station in range(1, 129):
		print(station, smn_ii[station][0])

		if station in skip_list:
			continue

		yy=smn_ii[station][1]
		xx=smn_ii[station][2]
		
		arq_i  = xr.open_dataset('{0}/OBS/WS-SA/SMN/daily/nc/pre/'.format(path) + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[station][0]))
		data_i = arq_i['pre']
		time_i = data_i.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_i  = time_i.groupby('time.month').mean('time')
		mean_i.append(var_i.values)

		arq_ii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'pre_{0}_CRU_mon_{1}_lonlat.nc'.format(domain, dt))
		data_ii = arq_ii['pre']
		data_ii = data_ii.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		time_ii = data_ii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_ii  = time_ii.groupby('time.month').mean('time')
		mean_ii.append(var_ii.values)

		arq_iii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'precip_{0}_CPC_mon_{1}_lonlat.nc'.format(domain, dt))
		data_iii = arq_iii['precip']
		data_iii = data_iii.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		time_iii = data_iii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_iii  = time_iii.groupby('time.month').mean('time')
		mean_iii.append(var_iii.values)
		
		arq_iv  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'tp_{0}_ERA5_mon_{1}_lonlat.nc'.format(domain, dt))
		data_iv = arq_iv['tp']
		data_iv = data_iv.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		time_iv = data_iv.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_iv  = time_iv.groupby('time.month').mean('time')
		mean_iv.append(var_iv.values)
		
		arq_v  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/'.format(path) + 'pr_{0}_RegCM5_mon_{1}_lonlat.nc'.format(domain, dt))
		data_v = arq_v['pr']
		data_v = data_v.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		time_v = data_v.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_v  = time_v.groupby('time.month').mean('time')
		mean_v.append(var_v.values)
			
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v



# Import model and obs dataset
clim_i_x, clim_ii_x, clim_iii_x, clim_iv_x, clim_v_x = import_situ_i()			
clim_i_y, clim_ii_y, clim_iii_y, clim_iv_y, clim_v_y = import_situ_ii()			
clim_i_z, clim_ii_z, clim_iii_z, clim_iv_z, clim_v_z = import_situ_iii()			

inmet_smn = clim_i_x + clim_i_y + clim_i_z
cru = clim_ii_x + clim_ii_y + clim_ii_z
cpc = clim_iii_x + clim_iii_y + clim_iii_z
era5 = clim_iv_x + clim_iv_y + clim_iv_z
regcm5 = clim_v_x + clim_v_y + clim_v_z

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
cru_i,   cru_ii,    cru_iii,    cru_iv,    cru_v    = [], [], [], [], []
cpc_i,   cpc_ii,    cpc_iii,    cpc_iv,    cpc_v    = [], [], [], [], []
era5_i,  era5_ii,   era5_iii,   era5_iv,   era5_v   = [], [], [], [], []
regcm_i, regcm_ii,  regcm_iii,  regcm_iv,  regcm_v  = [], [], [], [], []

for c_i in count_i:
	inmet_i.append(inmet_smn[c_i])
	cru_i.append(cru[c_i])
	cpc_i.append(cpc[c_i])
	era5_i.append(era5[c_i])
	regcm_i.append(regcm5[c_i])

for c_ii in count_ii:
	inmet_ii.append(inmet_smn[c_ii])
	cru_ii.append(cru[c_ii])
	cpc_ii.append(cpc[c_ii])
	era5_ii.append(era5[c_ii])
	regcm_ii.append(regcm5[c_ii])
	
for c_iii in count_iii:
	inmet_iii.append(inmet_smn[c_iii])
	cru_iii.append(cru[c_iii])
	cpc_iii.append(cpc[c_iii])
	era5_iii.append(era5[c_iii])
	regcm_iii.append(regcm5[c_iii])
	
for c_iv in count_iv:
	inmet_iv.append(inmet_smn[c_iv])
	cru_iv.append(cru[c_iv])
	cpc_iv.append(cpc[c_iv])
	era5_iv.append(era5[c_iv])
	regcm_iv.append(regcm5[c_iv])
	
for c_v in count_v:
	inmet_v.append(inmet_smn[c_v])
	cru_v.append(cru[c_v])
	cpc_v.append(cpc[c_v])
	era5_v.append(era5[c_v])
	regcm_v.append(regcm5[c_v])
	
inmet_c_i = np.nanmean(inmet_i, axis=0)
cru_c_i   = np.nanmean(cru_i, axis=0)
cpc_c_i   = np.nanmean(cpc_i, axis=0)
era5_c_i  = np.nanmean(era5_i, axis=0)
regcm_c_i = np.nanmean(regcm_i, axis=0)

inmet_c_ii = np.nanmean(inmet_ii, axis=0)
cru_c_ii   = np.nanmean(cru_ii, axis=0)
cpc_c_ii   = np.nanmean(cpc_ii, axis=0)
era5_c_ii  = np.nanmean(era5_ii, axis=0)
regcm_c_ii = np.nanmean(regcm_ii, axis=0)

inmet_c_iii = np.nanmean(inmet_iii, axis=0)
cru_c_iii   = np.nanmean(cru_iii, axis=0)
cpc_c_iii   = np.nanmean(cpc_iii, axis=0)
era5_c_iii  = np.nanmean(era5_iii, axis=0)
regcm_c_iii = np.nanmean(regcm_iii, axis=0)

inmet_c_iv = np.nanmean(inmet_iv, axis=0)
cru_c_iv   = np.nanmean(cru_iv, axis=0)
cpc_c_iv   = np.nanmean(cpc_iv, axis=0)
era5_c_iv  = np.nanmean(era5_iv, axis=0)
regcm_c_iv = np.nanmean(regcm_iv, axis=0)

inmet_c_v = np.nanmean(inmet_v, axis=0)
cru_c_v   = np.nanmean(cru_v, axis=0)
cpc_c_v   = np.nanmean(cpc_v, axis=0)
era5_c_v  = np.nanmean(era5_v, axis=0)
regcm_c_v = np.nanmean(regcm_v, axis=0)

inmet_c_ii_iii_iv = np.nanmean([inmet_c_ii, inmet_c_iii, inmet_c_iv], axis=0)
cru_c_ii_iii_iv   = np.nanmean([cru_c_ii,   cru_c_iii,   cru_c_iv], axis=0)
cpc_c_ii_iii_iv   = np.nanmean([cpc_c_ii,   cpc_c_iii,   cpc_c_iv], axis=0)
gpcp_c_ii_iii_iv  = np.nanmean([gpcp_c_ii,  gpcp_c_iii,  gpcp_c_iv], axis=0)
era5_c_ii_iii_iv  = np.nanmean([era5_c_ii,  era5_c_iii,  era5_c_iv], axis=0)
regcm_c_ii_iii_iv = np.nanmean([regcm_c_ii, regcm_c_iii, regcm_c_iv], axis=0)

# Plot figure
fig = plt.figure()
time = np.arange(0.5, 12 + 0.5)
font_size = 8
# Plot figure
fig = plt.figure(figsize=(6, 9))
time = np.arange(0.5, 12 + 0.5)
font_size = 8

ax = fig.add_subplot(3, 1, 1)
plt.plot(time, inmet_c_i, linewidth=1., color='blue',   markersize=2, markerfacecolor='white', marker='^', label='INMET+SMN')
plt.plot(time, cru_c_i,   linewidth=1., color='red',  markersize=2, markerfacecolor='white', marker='+', label='CRU')
plt.plot(time, cpc_c_i,   linewidth=1., color='magenta', markersize=2, markerfacecolor='white', marker='*', label='CPC')
plt.plot(time, era5_c_i,  linewidth=1., color='green', markersize=2, markerfacecolor='white', marker='s', label='ERA5')
plt.plot(time, regcm_c_i, linewidth=1., color='black',   markersize=2, markerfacecolor='white', marker='o', label='RegCM5')
plt.title('(a) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(0, 10)
plt.yticks(np.arange(0, 11, 1), fontsize=font_size)
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'), fontsize=font_size)
plt.grid(linestyle='--')
plt.legend(loc=1, ncol=3, fontsize=font_size, shadow=True)

ax = fig.add_subplot(3, 1, 2)
plt.plot(time, inmet_c_ii_iii_iv, linewidth=1., color='blue',   markersize=2, markerfacecolor='white', marker='^', label='INMET+SMN')
plt.plot(time, cru_c_ii_iii_iv,   linewidth=1., color='red',  markersize=2, markerfacecolor='white', marker='+', label='CRU')
plt.plot(time, cpc_c_ii_iii_iv,   linewidth=1., color='magenta', markersize=2, markerfacecolor='white', marker='*', label='CPC')
plt.plot(time, era5_c_ii_iii_iv,  linewidth=1., color='green', markersize=2, markerfacecolor='white', marker='s', label='ERA5')
plt.plot(time, regcm_c_ii_iii_iv, linewidth=1., color='black',  markersize=2, markerfacecolor='white', marker='o', label='RegCM5')
plt.title('(b) Cluster II-III-IV', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.ylim(0, 12)
plt.yticks(np.arange(0, 13, 1), fontsize=font_size)
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'), fontsize=font_size)
plt.grid(linestyle='--')

ax = fig.add_subplot(3, 1, 3)
plt.plot(time, inmet_c_v, linewidth=1., color='blue',   markersize=2, markerfacecolor='white', marker='^', label='INMET+SMN')
plt.plot(time, cru_c_v,   linewidth=1., color='red',  markersize=2, markerfacecolor='white', marker='+', label='CRU')
plt.plot(time, cpc_c_v,   linewidth=1., color='magenta', markersize=2, markerfacecolor='white', marker='*', label='CPC')
plt.plot(time, era5_c_v,  linewidth=1., color='green', markersize=2, markerfacecolor='white', marker='s', label='ERA5')
plt.plot(time, regcm_c_v, linewidth=1., color='black',   markersize=2, markerfacecolor='white', marker='o', label='RegCM5')
plt.title('(c) Cluster V', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('Months', fontsize=font_size, fontweight='bold')
plt.ylim(0, 8)
plt.yticks(np.arange(0, 9, 1), fontsize=font_size)
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'), fontsize=font_size)
plt.grid(linestyle='--')

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_annual_cycle_{0}_{1}_RegCM5_2018-2021.png'.format(var, domain)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
