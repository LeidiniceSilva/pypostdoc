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

	arq   = '{0}/user/mdasilva/EUR-11/postproc/obs/{1}_{2}_FPS_{3}_{4}_lonlat.nc'.format(path, param, domain, dataset, dt) 
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]

	return mean


def import_rcm(param, dataset):

	arq   = '{0}/user/mdasilva/EUR-11/postproc/rcm/{1}_{2}_FPS_{3}_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]

	return mean
	
	
# Import model and obs dataset
dict_var = {'pr': ['precip', 'rr', 'pr']}

cpc_jan = import_obs(dict_var[var][0], 'CPC')
noto_jan = import_rcm(var, 'NoTo-Europe')
wsm5_jan = import_rcm(var, 'WSM5-Europe')
wsm7_jan = import_rcm(var, 'WSM7-Europe')
wdm7_jan = import_rcm(var, 'WDM7-Europe')

# Convert array in list
cpc_jan_list = cpc_jan.flatten()
noto_jan_list = noto_jan.flatten()
wsm5_jan_list = wsm5_jan.flatten()
wsm7_jan_list = wsm7_jan.flatten()
wdm7_jan_list = wdm7_jan.flatten()

# Round values
cpc_jan_round = np.round(cpc_jan_list,0)
noto_jan_round = np.round(noto_jan_list,0)
wsm5_jan_round = np.round(wsm5_jan_list,0)
wsm7_jan_round = np.round(wsm7_jan_list,0)
wdm7_jan_round = np.round(wdm7_jan_list,0)

# Filter 0 mm/day
cpc_jan_filter = cpc_jan_round[cpc_jan_round > 0.]
noto_jan_filter = noto_jan_round[noto_jan_round > 0.]
wsm5_jan_filter = wsm5_jan_round[wsm5_jan_round > 0.]
wsm7_jan_filter = wsm7_jan_round[wsm7_jan_round > 0.]
wdm7_jan_filter = wdm7_jan_round[wdm7_jan_round > 0.]

# Compute frequency
x_pdf_cpc, pdf_cpc = np.unique(cpc_jan_filter, return_counts=True) 
x_pdf_noto, pdf_noto = np.unique(noto_jan_filter, return_counts=True) 
x_pdf_wsm5, pdf_wsm5 = np.unique(wsm5_jan_filter, return_counts=True) 
x_pdf_wsm7, pdf_wsm7 = np.unique(wsm7_jan_filter, return_counts=True) 
x_pdf_wdm7, pdf_wdm7 = np.unique(wdm7_jan_filter, return_counts=True) 

# Plot figure
fig = plt.figure()
font_size = 10

ax = fig.add_subplot(1, 1, 1)  
plt.plot(x_pdf_cpc, pdf_cpc, marker='o', markersize=3, mfc='black', mec='black', alpha=0.70, linestyle='None', label='CPC')
plt.plot(x_pdf_noto, pdf_noto, marker='o', markersize=3, mfc='blue', mec='blue', alpha=0.70, linestyle='None', label='NoTo')
plt.plot(x_pdf_wsm5, pdf_wsm5, marker='o', markersize=3, mfc='red', mec='red', alpha=0.70, linestyle='None', label='WSM5')
plt.plot(x_pdf_wsm7, pdf_wsm7, marker='o', markersize=3, mfc='magenta', mec='magenta', alpha=0.70, linestyle='None', label='WSM7')
plt.plot(x_pdf_wdm7, pdf_wdm7, marker='o', markersize=3, mfc='green', mec='green', alpha=0.70, linestyle='None', label='WDM7')

plt.title('(a)', loc='left', fontsize=font_size, fontweight='bold') 
plt.yscale('log')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.xlabel('Precipitation (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.legend(loc=1, ncol=2, fontsize=font_size, shadow=True)

# Path out to save figure
path_out = '{0}/user/mdasilva/EUR-11/figs'.format(path)
name_out = 'pyplt_pdf_daily_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

