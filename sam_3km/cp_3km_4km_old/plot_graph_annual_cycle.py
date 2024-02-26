# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot bias maps"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

var = 'pr'
path = '/marconi/home/userexternal/mdasilva'


def import_obs(param, domain, dataset, period):

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/{1}_{2}_{3}_{4}_2018-2021_lonlat.nc'.format(path, param, domain, dataset, period)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	
	value  = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	clim = []		
	for mon in range(0, 11 + 1):
		clim.append(np.nanmean(value[mon::12], axis=0))
					
	return clim


def import_cp_3km(param, domain, dataset, period):

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/{1}_{2}_{3}_{4}_2018-2021_lonlat.nc'.format(path, param, domain, dataset, period)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	
	if param == 'pr':
		value  = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
		clim = []		
		for mon in range(0, 11 + 1):
			clim.append(np.nanmean(value[mon::12], axis=0))
	else:
		value  = np.nanmean(np.nanmean(var[:][:,0,:,:], axis=1), axis=1)
		clim = []		
		for mon in range(0, 11 + 1):
			clim.append(np.nanmean(value[mon::12], axis=0))
			
	return clim


def import_cp_4km(param, domain, dataset, period):

	arq   = '{0}/user/mdasilva/CSAM-4/post_evaluate/{1}_{2}_{3}_{4}_2018-2021_lonlat.nc'.format(path, param, domain, dataset, period)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value  = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	value  = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	clim = []		
	for mon in range(0, 11 + 1):
		clim.append(np.nanmean(value[mon::12], axis=0))
						
	return clim
	

# Import model and obs dataset
if var == 'pr':
	dict_var = {'pr': ['pre', 'precip', 'sat_gauge_precip', 'tp']}
	cru = import_obs(dict_var[var][0], 'SESA-3km', 'CRU', 'mon')
	cpc = import_obs(dict_var[var][1], 'SESA-3km', 'CPC', 'mon')
	gpcp = import_obs(dict_var[var][2], 'SESA-3km', 'GPCP', 'mon')
	era5 = import_obs(dict_var[var][3], 'SESA-3km', 'ERA5', 'mon')
	cp_3km = import_cp_3km(var, 'SESA-3km', 'RegCM5', 'mon')
	cp_4km = import_cp_4km(var, 'SESA-4', 'RegCM5', 'mon')

else:
	dict_var = {'tas': ['tmp', 't2m']}
	cru = import_obs(dict_var[var][0], 'SESA-3km', 'CRU', 'mon')
	era5 = import_obs(dict_var[var][1], 'SESA-3km', 'ERA5', 'mon')
	cp_3km = import_cp_3km(var, 'SESA-3km', 'RegCM5', 'mon')
	cp_4km = import_cp_4km(var, 'SESA-4', 'RegCM5', 'mon')
	
# Plot figure
fig = plt.figure()
time = np.arange(0.5, 12 + 0.5)
font_size = 8

if var == 'pr':
	plt1 = plt.plot(time, cru, time, cpc, time, gpcp, time, era5, time, cp_3km, time, cp_4km)
	plt.title(u'a) SESA', loc='left', fontweight='bold', fontsize=8)
	l1, l2, l3, l4, l5, l6 = plt1
	plt.setp(l1, linewidth=1., linestyle='-', markersize=3, marker='o', markerfacecolor='white', color='gray')
	plt.setp(l2, linewidth=1., linestyle='-', markersize=3, marker='o', markerfacecolor='white', color='green')
	plt.setp(l3, linewidth=1., linestyle='-', markersize=3, marker='o', markerfacecolor='white', color='orange')
	plt.setp(l4, linewidth=1., linestyle='-', markersize=3, marker='o', markerfacecolor='white',  color='black')
	plt.setp(l5, linewidth=1., linestyle='-', markersize=3, marker='o', markerfacecolor='white', color='red')  
	plt.setp(l6, linewidth=1., linestyle='-', markersize=3, marker='o', markerfacecolor='white', color='blue')  
	plt.ylim(0, 10)
	plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=font_size)
	plt.yticks(np.arange(0, 11, 1), fontsize=font_size)
	plt.xlabel('Months', fontsize=font_size, fontweight='bold')
	plt.ylabel('Precipitation (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
	plt.grid(linestyle='--')
	plt.axvline(4.5, linewidth=1., linestyle='-', color='black')
	plt.axvline(8.5, linewidth=1., linestyle='-', color='black')
	plt.legend(plt1, ['CRU', 'CPC', 'GPCP', 'ERA5', 'RegCM5-3km', 'RegCM5-4km'], fontsize=font_size, ncol=1, loc=1, shadow=True)
else:
	plt1 = plt.plot(time, cru, time, era5, time, cp_3km, time, cp_4km,)
	plt.title(u'a) SESA', loc='left', fontweight='bold', fontsize=8)
	l1, l2, l3, l4 = plt1
	plt.setp(l1, linewidth=1., linestyle='-', markersize=3, marker='o', markerfacecolor='white', color='gray')
	plt.setp(l2, linewidth=1., linestyle='-', markersize=3, marker='o', markerfacecolor='white',  color='black')
	plt.setp(l3, linewidth=1., linestyle='-', markersize=3, marker='o', markerfacecolor='white', color='red')  
	plt.setp(l4, linewidth=1., linestyle='-', markersize=3, marker='o', markerfacecolor='white', color='blue')  
	plt.ylim(8, 30)
	plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=font_size)
	plt.yticks(np.arange(8, 32, 2), fontsize=font_size)
	plt.xlabel('Months', fontsize=font_size, fontweight='bold')
	plt.ylabel('Temperature (°C)', fontsize=font_size, fontweight='bold')
	plt.grid(linestyle='--')
	plt.axvline(4.5, linewidth=1., linestyle='-', color='black')
	plt.axvline(8.5, linewidth=1., linestyle='-', color='black')
	plt.legend(plt1, ['CRU', 'ERA5', 'RegCM5-3km', 'RegCM5-4km'], fontsize=font_size, ncol=1, loc=1, shadow=True)

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cp_3km-4km'.format(path)
name_out = 'pyplt_graph_annual_cycle_{0}_CP-RegCM5_SAM-3km_2018-2021.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
