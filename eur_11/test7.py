# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot clim maps"

import os
import netCDF4
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap

var = 'pr'
dt = '20000101'
domain = 'EUR-11'
path = '/marconi/home/userexternal/mdasilva'


def import_obs(param, dataset):

	arq   = '{0}/user/mdasilva/EUR-11/post_evaluate/obs/{1}_{2}_FPS_{3}_{4}_lonlat.nc'.format(path, param, domain, dataset, dt) 
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]

	return mean


def import_rcm(exp, param, dataset):

	arq   = '{0}/user/mdasilva/EUR-11/post_evaluate/rcm/{1}/{2}_{3}_FPS_{4}_{5}_lonlat.nc'.format(path, exp, param, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]

	return mean
	
	
# Import model and obs dataset
dict_var = {'pr': ['rr', 'precipitation', 'precip', ]}

eobs_jan = import_obs(dict_var[var][0], 'EOBS')
cpc_jan = import_obs(dict_var[var][2], 'CPC')
wdm7_jan_v1 = import_rcm('wdm7-Europe_v1', var, 'RegCM5')
wdm7_jan_v2 = import_rcm('wdm7-Europe_v2', var, 'RegCM5')
wdm7_jan_v3 = import_rcm('wdm7-Europe_v3', var, 'RegCM5')
wdm7_jan_v4 = import_rcm('wdm7-Europe_v4', var, 'RegCM5')

# Convert array in list
eobs_jan_list = eobs_jan.flatten()
cpc_jan_list = cpc_jan.flatten()
wdm7_jan_v1_list = wdm7_jan_v1.flatten()
wdm7_jan_v2_list = wdm7_jan_v2.flatten()
wdm7_jan_v3_list = wdm7_jan_v3.flatten()
wdm7_jan_v4_list = wdm7_jan_v4.flatten()

# Round values
eobs_jan_round = np.round(eobs_jan_list,0)
cpc_jan_round = np.round(cpc_jan_list,0)
wdm7_jan_v1_round = np.round(wdm7_jan_v1_list,0)
wdm7_jan_v2_round = np.round(wdm7_jan_v2_list,0)
wdm7_jan_v3_round = np.round(wdm7_jan_v3_list,0)
wdm7_jan_v4_round = np.round(wdm7_jan_v4_list,0)

# Filter 0 mm/day
filter_eobs = eobs_jan_round[eobs_jan_round > 0.]
filter_cpc = cpc_jan_round[cpc_jan_round > 0.]
filter_wdm7_v1 = wdm7_jan_v1_round[wdm7_jan_v1_round > 0.]
filter_wdm7_v2 = wdm7_jan_v2_round[wdm7_jan_v2_round > 0.]
filter_wdm7_v3 = wdm7_jan_v3_round[wdm7_jan_v3_round > 0.]
filter_wdm7_v4 = wdm7_jan_v4_round[wdm7_jan_v4_round > 0.]

# Compute frequency
x_pdf_eobs, pdf_eobs = np.unique(filter_eobs, return_counts=True) 
x_pdf_cpc, pdf_cpc = np.unique(filter_cpc, return_counts=True) 
x_pdf_wdm7_v1, pdf_wdm7_v1 = np.unique(filter_wdm7_v1, return_counts=True) 
x_pdf_wdm7_v2, pdf_wdm7_v2 = np.unique(filter_wdm7_v2, return_counts=True) 
x_pdf_wdm7_v3, pdf_wdm7_v3 = np.unique(filter_wdm7_v3, return_counts=True) 
x_pdf_wdm7_v4, pdf_wdm7_v4 = np.unique(filter_wdm7_v4, return_counts=True) 


# Plot figure
fig = plt.figure()
font_size = 8

ax = fig.add_subplot(1, 1, 1)  
plt.plot(x_pdf_eobs, pdf_eobs, marker='o', markersize=4, mfc='red', mec='red', linestyle='None', label='EOB')
plt.plot(x_pdf_cpc, pdf_cpc, marker='o', markersize=4, mfc='black', mec='black', linestyle='None', label='CPC')
plt.plot(x_pdf_wdm7_v1, pdf_wdm7_v1, marker='o', markersize=4, mfc='blue', mec='blue', alpha=0.65, linestyle='None', label='WDM7_v1(ctrl)')
plt.plot(x_pdf_wdm7_v2, pdf_wdm7_v2, marker='o', markersize=4, mfc='orange', mec='orange', alpha=0.65, linestyle='None', label='WDM7_v2(more ccn)')
plt.plot(x_pdf_wdm7_v3, pdf_wdm7_v3, marker='o', markersize=4, mfc='magenta', mec='magenta', alpha=0.65, linestyle='None', label='WDM7_v3(less ccn)')
plt.plot(x_pdf_wdm7_v4, pdf_wdm7_v4, marker='o', markersize=4, mfc='green', mec='green', alpha=0.65, linestyle='None', label='WDM7_v4(less ccn2)')

plt.title('(a)', loc='left', fontsize=font_size, fontweight='bold') 
plt.yscale('log')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.xlabel('Precipitation (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.legend(loc=1, ncol=2, fontsize=font_size, shadow=True)

# Path out to save figure
path_out = '{0}/user/mdasilva/EUR-11/figs'.format(path)
name_out = 'pyplt_pdf_{0}_{1}_RegCM5_WDM7_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

