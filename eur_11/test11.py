# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot pdf"

import os
import netCDF4
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt

var = 'pr'
freq = 'day'
domain = 'EUR-11'
dt = '1970-1970'
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11'


def import_obs(param, dataset):

	arq   = '{0}/postproc/obs/{1}_{2}_FPS_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, freq, dt) 
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	return value


def import_rcm(param, dataset):

	arq   = '{0}/postproc/rcm/{1}_{2}_FPS_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, freq, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	return value


def compute_pdf(data):

	data_list = data.flatten()
	data_round = np.round(data_list,0)
	data_filter = data_round[data_round > 0.]
	x_pdf, pdf = np.unique(data_filter, return_counts=True) 

	return x_pdf, pdf


# Import model and obs dataset
dict_var = {'pr': ['precip', 'rr', 'tp']}

eobs = import_obs(dict_var[var][1], 'EOBS')
noto_v1 = import_rcm(var, 'NoTo-Europe_RegCM5')
noto_v2 = import_rcm(var, 'NoTo-Europe_cordex5_RegCM5')
wdm7 = import_rcm(var, 'WDM7-Europe_RegCM5')
wsm7 = import_rcm(var, 'WSM7-Europe_RegCM5')
wsm5 = import_rcm(var, 'WSM5-Europe_RegCM5')

# Import pdf
x_pdf_eobs, pdf_eobs = compute_pdf(eobs) 
x_pdf_noto_v1, pdf_noto_v1 = compute_pdf(noto_v1) 
x_pdf_noto_v2, pdf_noto_v2 = compute_pdf(noto_v2)
x_pdf_wdm7, pdf_wdm7 = compute_pdf(wdm7) 
x_pdf_wsm7, pdf_wsm7 = compute_pdf(wsm7) 
x_pdf_wsm5, pdf_wsm5 = compute_pdf(wsm5) 

# Plot figure
fig = plt.figure()
font_size = 10

ax = fig.add_subplot(1, 1, 1)  
plt.plot(x_pdf_eobs, pdf_eobs, marker='o', markersize=3, mfc='black', mec='black', alpha=0.70, linestyle='None', label='EOBS')
plt.plot(x_pdf_noto_v1, pdf_noto_v1, marker='o', markersize=3, mfc='blue', mec='blue', alpha=0.70, linestyle='None', label='NoTo (CORE)')
plt.plot(x_pdf_noto_v2, pdf_noto_v2, marker='o', markersize=3, mfc='gray', mec='gray', alpha=0.70, linestyle='None', label='NoTo (CORDEX5)')
plt.plot(x_pdf_wdm7, pdf_wdm7, marker='o', markersize=3, mfc='green', mec='green', alpha=0.70, linestyle='None', label='WDM7')
plt.plot(x_pdf_wsm7, pdf_wsm7, marker='o', markersize=3, mfc='magenta', mec='magenta', alpha=0.70, linestyle='None', label='WSM7')
plt.plot(x_pdf_wsm5, pdf_wsm5, marker='o', markersize=3, mfc='red', mec='red', alpha=0.70, linestyle='None', label='WSM5')

plt.title('(a)', loc='left', fontsize=font_size, fontweight='bold') 
plt.yscale('log')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.xlabel('Precipitation (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.legend(loc=1, ncol=2, fontsize=font_size, shadow=True)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_graph_pdf_{0}_{1}_RegCM5_{2}_{3}_v1-v2.png'.format(var, domain, freq, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()



