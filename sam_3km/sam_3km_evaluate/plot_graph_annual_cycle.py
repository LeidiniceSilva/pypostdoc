# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot annual cycle"

import os
import netCDF4
import numpy as np
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

		yy=inmet[i][2]
		xx=inmet[i][3]

		arq_i  = xr.open_dataset('{0}/OBS/WS-SA/INMET/nc/hourly/{1}/'.format(path, param_i) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[station][0]))
		data_i = arq_i['pre']
		time_i = data_i.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_i  = time_i.groupby('time.month').mean('time')
		mean_i.append(var_i.values*24)

		arq_ii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'pre_{0}_CRU_mon_{1}_lonlat.nc'.format(domain, dt))
		data_ii = arq_ii['pre']
		data_ii = data_ii.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).inmet(('lat','lon'))
		time_ii = data_ii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_ii  = time_ii.groupby('time.month').mean('time')
		mean_ii.append(var_ii.values)

		arq_iii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'precip_{0}_CPC_mon_{1}_lonlat.nc'.format(domain, dt))
		data_iii = arq_iii['precip']
		data_iii = data_iii.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).inmet(('lat','lon'))
		time_iii = data_iii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_iii  = time_iii.groupby('time.month').mean('time')
		mean_iii.append(var_iii.values)
		
		arq_iv  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'tp_{0}_ERA5_mon_{1}_lonlat.nc'.format(domain, dt))
		data_iv = arq_iv['tp']
		data_iv = data_iv.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).inmet(('lat','lon'))
		time_iv = data_iv.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_iv  = time_iv.groupby('time.month').mean('time')
		mean_iv.append(var_iv.values)
		
		arq_v  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/'.format(path) + 'pr_{0}_RegCM5_mon_{1}_lonlat.nc'.format(domain, dt))
		data_v = arq_v[param_ii]
		data_v = data_v.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).inmet(('lat','lon'))
		time_v = data_v.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_v  = time_v.groupby('time.month').mean('time')
		mean_v.append(var_v.values)
			
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v
	

def import_situ_ii(param_i, param_ii, domain, dataset):
	
	mean_i, mean_ii, mean_iii, mean_iv, mean_v = [], [], [], [], []
	
	for station in range(1, 73):
		print(station, smn_i[station][0])

		yy=smn_i[i][1]
		xx=smn_i[i][2]

		arq_i  = xr.open_dataset('{0}/OBS/WS-SA/SMN/hourly/nc/'.format(path) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(param_i, smn_i[station][0]))
		data_i = arq_i[param_i]
		time_i = data_i.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_i  = time_i.groupby('time.month').mean('time')
		mean_i.append(var_i.values*24)

		arq_ii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'pre_{0}_CRU_mon_{1}_lonlat.nc'.format(domain, dt))
		data_ii = arq_ii['pre']
		data_ii = data_ii.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).inmet(('lat','lon'))
		time_ii = data_ii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_ii  = time_ii.groupby('time.month').mean('time')
		mean_ii.append(var_ii.values)

		arq_iii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'precip_{0}_CPC_mon_{1}_lonlat.nc'.format(domain, dt))
		data_iii = arq_iii['precip']
		data_iii = data_iii.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).inmet(('lat','lon'))
		time_iii = data_iii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_iii  = time_iii.groupby('time.month').mean('time')
		mean_iii.append(var_iii.values)
		
		arq_iv  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'tp_{0}_ERA5_mon_{1}_lonlat.nc'.format(domain, dt))
		data_iv = arq_iv['tp']
		data_iv = data_iv.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).inmet(('lat','lon'))
		time_iv = data_iv.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_iv  = time_iv.groupby('time.month').mean('time')
		mean_iv.append(var_iv.values)
		
		arq_v  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/'.format(path) + 'pr_{0}_RegCM5_mon_{1}_lonlat.nc'.format(domain, dt))
		data_v = arq_v[param_ii]
		data_v = data_v.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).inmet(('lat','lon'))
		time_v = data_v.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_v  = time_v.groupby('time.month').mean('time')
		mean_v.append(var_v.values)
			
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v
	

def import_situ_iii(param_i, param_ii, domain, dataset):
	
	mean_i, mean_ii, mean_iii, mean_iv, mean_v = [], [], [], [], []
	
	skip_list = [86,87,88,89,90,91,95,96,97,98,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,
	117,118,119,120,121,122,123,124,125,126,127,128]	
		
	for station in range(1, 129):
		print(station, smn_ii[station][0])

		if station in skip_list:
			continue

		yy=smn_ii[i][1]
		xx=smn_ii[i][2]
		
		arq_i  = xr.open_dataset('{0}/OBS/WS-SA/SMN/daily/nc/{1}/'.format(path, param_i) + '{0}_{1}_D_1979-01-01_2021-12-31.nc'.format(param_i, smn_ii[station][0]))
		data_i = arq_i[param_i]
		time_i = data_i.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_i  = time_i.groupby('time.month').mean('time')
		mean_i.append(var_i.values)

		arq_ii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'pre_{0}_CRU_mon_{1}_lonlat.nc'.format(domain, dt))
		data_ii = arq_ii['pre']
		data_ii = data_ii.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).inmet(('lat','lon'))
		time_ii = data_ii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_ii  = time_ii.groupby('time.month').mean('time')
		mean_ii.append(var_ii.values)

		arq_iii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'precip_{0}_CPC_mon_{1}_lonlat.nc'.format(domain, dt))
		data_iii = arq_iii['precip']
		data_iii = data_iii.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).inmet(('lat','lon'))
		time_iii = data_iii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_iii  = time_iii.groupby('time.month').mean('time')
		mean_iii.append(var_iii.values)
		
		arq_iv  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'tp_{0}_ERA5_mon_{1}_lonlat.nc'.format(domain, dt))
		data_iv = arq_iv['tp']
		data_iv = data_iv.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).inmet(('lat','lon'))
		time_iv = data_iv.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_iv  = time_iv.groupby('time.month').mean('time')
		mean_iv.append(var_iv.values)
		
		arq_v  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/'.format(path) + 'pr_{0}_RegCM5_mon_{1}_lonlat.nc'.format(domain, dt))
		data_v = arq_v[param_ii]
		data_v = data_v.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).inmet(('lat','lon'))
		time_v = data_v.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_v  = time_v.groupby('time.month').mean('time')
		mean_v.append(var_v.values)
			
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v



# Import model and obs dataset
# Import dataset
clim_i_x, clim_ii_x, clim_iii_x, clim_iv_x, clim_v_x = import_inmet()			
clim_i_y, clim_ii_y, clim_iii_y, clim_iv_y, clim_v_y = import_smn_i()			
clim_i_z, clim_ii_z, clim_iii_z, clim_iv_z, clim_v_z = import_smn_ii()			

inmet_smn = clim_i_x + clim_i_y + clim_i_z
cru = clim_ii_x + clim_ii_y + clim_ii_z
cpc = clim_iii_x + clim_iii_y + clim_iii_z
era5 = clim_iv_x + clim_iv_y + clim_iv_z
regcm5 = clim_v_x + clim_v_y + clim_v_z





# Plot figure
fig = plt.figure()
time = np.arange(0.5, 12 + 0.5)
font_size = 8

if var == 'pr':
	plt1 = plt.plot(time, cru, time, cpc, time, gpcp, time, era5, time, regcm)
	plt.title(u'a) SESA', loc='left', fontweight='bold', fontsize=8)
	l1, l2, l3, l4, l5 = plt1
	plt.setp(l1, linewidth=1.5, linestyle='-', color='blue')
	plt.setp(l2, linewidth=1.5, linestyle='-', color='green')
	plt.setp(l3, linewidth=1.5, linestyle='-', color='orange',)
	plt.setp(l4, linewidth=1.5, linestyle='-', color='black')
	plt.setp(l5, linewidth=1.5, linestyle='-', color='red')   
	plt.ylim(0, 10)
	plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=font_size)
	plt.yticks(np.arange(0, 11, 1), fontsize=font_size)
	plt.xlabel('Months', fontsize=font_size, fontweight='bold')
	plt.ylabel('Precipitation (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
	plt.grid(linestyle='--')
	plt.axvline(4.5, linewidth=1., linestyle='-', color='black')
	plt.axvline(8.5, linewidth=1., linestyle='-', color='black')
	plt.legend(plt1, ['CRU', 'CPC', 'GPCP', 'ERA5', 'RegCM5'], fontsize=font_size, ncol=1, loc=1, shadow=True)
else:
	plt1 = plt.plot(time, cru, time, era5, time, regcm)
	plt.title(u'a) SESA', loc='left', fontweight='bold', fontsize=8)
	l1, l2, l3 = plt1
	plt.setp(l1, linewidth=1.5, linestyle='-', color='blue')
	plt.setp(l2, linewidth=1.5, linestyle='-', color='black')
	plt.setp(l3, linewidth=1.5, linestyle='-', color='red')   
	plt.ylim(8, 30)
	plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=font_size)
	plt.yticks(np.arange(8, 32, 2), fontsize=font_size)
	plt.xlabel('Months', fontsize=font_size, fontweight='bold')
	plt.ylabel('Temperature (°C)', fontsize=font_size, fontweight='bold')
	plt.grid(linestyle='--')
	plt.axvline(4.5, linewidth=1., linestyle='-', color='black')
	plt.axvline(8.5, linewidth=1., linestyle='-', color='black')
	plt.legend(plt1, ['CRU', 'ERA5', 'RegCM5'], fontsize=font_size, ncol=1, loc=1, shadow=True)

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_annual_cycle_{0}_{1}_RegCM5_2018-2021.png'.format(var, domain)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
