# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot annual cycle"

import os
import netCDF4
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.pyplot as plt

var = 'pr'
domain = 'SESA-3km'
path = '/marconi/home/userexternal/mdasilva'


def import_obs(param, domain, dataset, period):

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/obs/{1}_{2}_{3}_{4}_2018-2021_lonlat.nc'.format(path, param, domain, dataset, period)
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

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/{1}_{2}_{3}_{4}_2018-2021_lonlat.nc'.format(path, param, domain, dataset, period)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]

	if param == 'pr':
		value  = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	else:
		value  = np.nanmean(np.nanmean(var[:][:,0,:,:], axis=1), axis=1)

	clim = []		
	for mon in range(0, 11 + 1):
		clim.append(np.nanmean(value[mon::12], axis=0))
		
	return clim
	
	
# Import model and obs dataset
if var == 'pr':
	dict_var = {'pr': ['pre', 'precip', 'sat_gauge_precip', 'tp']}
	cru = import_obs(dict_var[var][0], domain, 'CRU', 'mon')
	cpc = import_obs(dict_var[var][1], domain, 'CPC', 'mon')
	gpcp = import_obs(dict_var[var][2], domain, 'GPCP', 'mon')
	era5 = import_obs(dict_var[var][3], domain, 'ERA5', 'mon')
	regcm = import_cp_3km(var, domain, 'RegCM5', 'mon')
else:
	dict_var = {'tas': ['tmp', 't2m']}
	cru = import_obs(dict_var[var][0], domain, 'CRU', 'mon')
	era5 = import_obs(dict_var[var][1], domain, 'ERA5', 'mon')
	regcm = import_cp_3km(var, domain, 'RegCM5', 'mon')

datasets = {'CRU': cru,
'CPC': cpc,
'GPCP': gpcp,
'ERA5': era5,
'RegCM5': regcm}

datasets_df = pd.DataFrame(datasets,columns=['CRU', 'CPC', 'GPCP', 'ERA5', 'RegCM5'])
datasets_corr = datasets_df.corr()

# Plot figure
fig = plt.figure()
font_size = 10

ax = fig.add_subplot(1, 1, 1)
mask = np.triu(np.ones_like(datasets_corr, dtype=bool))
heatmap = sns.heatmap(datasets_corr, cmap='BrBG', vmin=0.9, vmax=1, annot=True, annot_kws={"size":font_size}, linewidths=.8, ax=ax)
heatmap.set_title('a) Correlation matrix of precipitation', fontdict={'fontsize':font_size}, loc='left', fontweight='bold')
heatmap.set_xlabel('Dataset', fontdict={'fontsize':font_size}, fontweight='bold')
heatmap.set_ylabel('Dataset', fontdict={'fontsize':font_size}, fontweight='bold')
	
# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_annual_cycle_corr_{0}_{1}_RegCM5_2018-2021.png'.format(var, domain)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
