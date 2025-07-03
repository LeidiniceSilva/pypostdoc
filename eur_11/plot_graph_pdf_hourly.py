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
freq = '1hr'
domain = 'EUR-11'
dt = '2000-2009'
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11'


def import_obs(param, dataset):

	arq   = '{0}/postproc/obs/{1}_{2}_{3}-FPS_{4}_{5}_lonlat.nc'.format(path, param, dataset, domain, freq, dt) 
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	return value


def import_rcm(param, dataset):

	arq   = '{0}/postproc/rcm/{1}_{2}-FPS_{3}_{4}_lonlat.nc'.format(path, param, dataset, freq, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	return value


def compute_pdf(data):

	step=1
	rain_min=0.1
	rain_max=500.0
	bins = np.arange(rain_min, rain_max + step, step)
	rain_hist = np.zeros(len(bins) - 1)
	for t in range(data.shape[0]):
		values = data[t, :, :]
		wet = values[values >= rain_min]
		hist, _ = np.histogram(wet, bins=bins)
		rain_hist += hist

	total = np.sum(rain_hist)
	rain_hist[rain_hist < 1] = np.nan  
	pdf = rain_hist / (total * step)

	bin_centers = (bins[:-1] + bins[1:]) / 2

	return pdf


# Import model and obs dataset
dict_var = {'pr': ['precip', 'rr', 'tp']}

era5 = import_obs(dict_var[var][2], 'ERA5')
noto = import_rcm(var, 'NoTo-Europe_RegCM5')
wdm7 = import_rcm(var, 'WDM7-Europe_RegCM5')
wsm7 = import_rcm(var, 'WSM7-Europe_RegCM5')
wsm5 = import_rcm(var, 'WSM5-Europe_RegCM5')

era5_pdf = compute_pdf(era5)
noto_pdf = compute_pdf(noto)
wdm7_pdf = compute_pdf(wdm7)
wsm7_pdf = compute_pdf(wsm7)
wsm5_pdf = compute_pdf(wsm5)

# Plot figure
fig = plt.figure()
font_size = 10

ax = fig.add_subplot(1, 1, 1)  
plt.plot(era5_pdf, marker='o', markersize=3, mfc='black', mec='black', alpha=0.70, linestyle='None', label='ERA5')
plt.plot(noto_pdf, pdf_noto, marker='o', markersize=3, mfc='blue', mec='blue', alpha=0.70, linestyle='None', label='NoTo')
plt.plot(wdm7_pdf, pdf_wdm7, marker='o', markersize=3, mfc='green', mec='green', alpha=0.70, linestyle='None', label='WDM7')
plt.plot(wsm7_pdf, pdf_wsm7, marker='o', markersize=3, mfc='magenta', mec='magenta', alpha=0.70, linestyle='None', label='WSM7')
plt.plot(wsm5_pdf, pdf_wsm5, marker='o', markersize=3, mfc='red', mec='red', alpha=0.70, linestyle='None', label='WSM5')

plt.title('(b)', loc='left', fontsize=font_size, fontweight='bold') 
plt.yscale('log')
plt.ylabel('Frequency (Log)', fontsize=font_size, fontweight='bold')
plt.xlabel('Precipitation (mm h$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.legend(loc=1, ncol=2, fontsize=font_size, shadow=True)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_graph_pdf_{0}_{1}_RegCM5_{2}_{3}.png'.format(var, domain, freq, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()



