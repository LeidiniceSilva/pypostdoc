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

var = 'pr'
freq= 'daily'
domain = 'CSAM-3'
idt, fdt = '2000', '2009'

if freq == 'hourly':
	dt = '1hr_{0}-{1}'.format('2000', '2005')
	legend = 'Precipitation (mm h$^-$$^1$)'

else:
	dt = 'day_{0}-{1}'.format('2000', '2005')
	legend = 'Precipitation (mm d$^-$$^1$)'
				
dict_ws = {
	   0: [-68.1193, -16.4897, 'La Paz'],
	   1: [-47.9000, -16.0000, 'Brazilia'],
	   2: [-70.6693, -33.4489, 'Santiago'],
	   3: [-57.5759, -25.2637, 'Asuncion'],
	   4: [-58.4004, -34.6051, 'Buenos Aires'],
	   5: [-56.0000, -34.0000, 'Montevideo']
	   }

path = '/marconi/home/userexternal/mdasilva'		

		
def import_obs(param, dataset):

	if freq == 'hourly':
		dt = '1hr_{0}-{1}'.format(idt, fdt)
	else:
		dt = 'day_{0}-{1}'.format(idt, fdt)
		
	ts = []
	for i in range(0, 6):
		yy=dict_ws[i][1]
		xx=dict_ws[i][0]

		arq  = xr.open_dataset('{0}/user/mdasilva/CORDEX/post_evaluate/obs/{1}_{2}_{3}_lonlat.nc'.format(path, param, dataset, dt))
		data = arq[param]
		var  = data.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		ts.append(var.values)

	return ts


def import_rcm(param, dataset):

	if freq == 'hourly':
		dt = '1hr_{0}-{1}'.format('2000', '2005')
	else:
		dt = 'day_{0}-{1}'.format('2000', '2005')
		
	ts = []
	for i in range(0, 6):
		yy=dict_ws[i][1]
		xx=dict_ws[i][0]
		
		arq  = xr.open_dataset('{0}/user/mdasilva/CORDEX/post_evaluate/rcm/{1}_{2}_{3}_lonlat.nc'.format(path, param, dataset, dt))
		data = arq[param]
		var  = data.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		ts.append(var.values)
		
	return ts


def import_pdf(database):

	x_pdf_, pdf_ = [], []
	for i in range(0, 6):

		round_values = np.round(database[i],0)
		filter_values = round_values[round_values > 0.]
		x_pdf, pdf = np.unique(filter_values, return_counts=True)
		x_pdf_.append(x_pdf)
		pdf_.append(pdf)
	
	return x_pdf_, pdf_
			

# Import model and obs dataset
dict_var = {'pr': ['pre', 'precip', 'cmorph', 'precipitation', 'pr']}

# Import model and obs dataset
if freq == 'hourly':
	ts_cmorph = import_obs(dict_var[var][2], 'CSAM-3_CMORPH')
	ts_era5   = import_obs(dict_var[var][4], 'CSAM-3_ERA5')
	ts_rcm3   = import_rcm(dict_var[var][4], 'CSAM-3_RegCM5')
else:
	ts_cpc    = import_obs(dict_var[var][1], 'CSAM-3_CPC')
	ts_cmorph = import_obs(dict_var[var][2], 'CSAM-3_CMORPH')
	ts_mswep  = import_obs(dict_var[var][3], 'CSAM-3_MSWEP')
	ts_era5   = import_obs(dict_var[var][4], 'CSAM-3_ERA5')
	ts_rcm3   = import_rcm(dict_var[var][4], 'CSAM-3_RegCM5')


# Compute frequency
x_pdf_cpc, pdf_cpc = import_pdf(ts_cpc)
x_pdf_cmorph, pdf_cmorph = import_pdf(ts_cmorph)
x_pdf_mswep, pdf_mswep = import_pdf(ts_mswep)
x_pdf_era5, pdf_era5 = import_pdf(ts_era5)
x_pdf_rcm3, pdf_rcm3 = import_pdf(ts_rcm3)


# Plot figure
fig = plt.figure(figsize=(14, 10))
time = np.arange(0.5, 12 + 0.5)
font_size = 8

ax = fig.add_subplot(3, 3, 1)
plt.plot(x_pdf_rcm3[0], pdf_rcm3[0],     marker='o', markersize=3, mfc='black',  mec='black',  alpha=0.75, linestyle='None', label='CPM3')
plt.plot(x_pdf_cpc[0], pdf_cpc[0],       marker='o', markersize=3, mfc='green',  mec='green',  alpha=0.75, linestyle='None', label='CPC')
plt.plot(x_pdf_cmorph[0], pdf_cmorph[0], marker='o', markersize=3, mfc='red',    mec='red',    alpha=0.75, linestyle='None', label='CMORPH')
plt.plot(x_pdf_mswep[0], pdf_mswep[0],   marker='o', markersize=3, mfc='violet', mec='violet', alpha=0.75, linestyle='None', label='MSWEP')
plt.plot(x_pdf_era5[0], pdf_era5[0],     marker='o', markersize=3, mfc='blue',   mec='blue',   alpha=0.75, linestyle='None', label='ERA5')
plt.title('(a) {0}'.format(dict_ws[0][2]), loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.legend(loc=1, ncol=2, fontsize=font_size, shadow=True)

ax = fig.add_subplot(3, 3, 3)
plt.plot(x_pdf_rcm3[1], pdf_rcm3[1],     marker='o', markersize=3, mfc='black',  mec='black',  alpha=0.75, linestyle='None', label='CPM3')
plt.plot(x_pdf_cpc[1], pdf_cpc[1],       marker='o', markersize=3, mfc='green',  mec='green',  alpha=0.75, linestyle='None', label='CPC')
plt.plot(x_pdf_cmorph[1], pdf_cmorph[1], marker='o', markersize=3, mfc='red',    mec='red',    alpha=0.75, linestyle='None', label='CMORPH')
plt.plot(x_pdf_mswep[1], pdf_mswep[1],   marker='o', markersize=3, mfc='violet', mec='violet', alpha=0.75, linestyle='None', label='MSWEP')
plt.plot(x_pdf_era5[1], pdf_era5[1],     marker='o', markersize=3, mfc='blue',   mec='blue',   alpha=0.75, linestyle='None', label='ERA5')
plt.title('(b) {0}'.format(dict_ws[1][2]), loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.ylabel('{0}', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 4)
plt.plot(x_pdf_rcm3[2], pdf_rcm3[2],     marker='o', markersize=3, mfc='black',  mec='black',  alpha=0.75, linestyle='None', label='CPM3')
plt.plot(x_pdf_cpc[2], pdf_cpc[2],       marker='o', markersize=3, mfc='green',  mec='green',  alpha=0.75, linestyle='None', label='CPC')
plt.plot(x_pdf_cmorph[2], pdf_cmorph[2], marker='o', markersize=3, mfc='red',    mec='red',    alpha=0.75, linestyle='None', label='CMORPH')
plt.plot(x_pdf_mswep[2], pdf_mswep[2],   marker='o', markersize=3, mfc='violet', mec='violet', alpha=0.75, linestyle='None', label='MSWEP')
plt.plot(x_pdf_era5[2], pdf_era5[2],     marker='o', markersize=3, mfc='blue',   mec='blue',   alpha=0.75, linestyle='None', label='ERA5')
plt.title('(c) {0}'.format(dict_ws[2][2]), loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.yscale('log')

ax = fig.add_subplot(3, 3, 6)
plt.plot(x_pdf_rcm3[3], pdf_rcm3[3],     marker='o', markersize=3, mfc='black',  mec='black',  alpha=0.75, linestyle='None', label='CPM3')
plt.plot(x_pdf_cpc[3], pdf_cpc[3],       marker='o', markersize=3, mfc='green',  mec='green',  alpha=0.75, linestyle='None', label='CPC')
plt.plot(x_pdf_cmorph[3], pdf_cmorph[3], marker='o', markersize=3, mfc='red',    mec='red',    alpha=0.75, linestyle='None', label='CMORPH')
plt.plot(x_pdf_mswep[3], pdf_mswep[3],   marker='o', markersize=3, mfc='violet', mec='violet', alpha=0.75, linestyle='None', label='MSWEP')
plt.plot(x_pdf_era5[3], pdf_era5[3],     marker='o', markersize=3, mfc='blue',   mec='blue',   alpha=0.75, linestyle='None', label='ERA5')
plt.title('(d) {0}'.format(dict_ws[3][2]), loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 7)
plt.plot(x_pdf_rcm3[4], pdf_rcm3[4],     marker='o', markersize=3, mfc='black',  mec='black',  alpha=0.75, linestyle='None', label='CPM3')
plt.plot(x_pdf_cpc[4], pdf_cpc[4],       marker='o', markersize=3, mfc='green',  mec='green',  alpha=0.75, linestyle='None', label='CPC')
plt.plot(x_pdf_cmorph[4], pdf_cmorph[4], marker='o', markersize=3, mfc='red',    mec='red',    alpha=0.75, linestyle='None', label='CMORPH')
plt.plot(x_pdf_mswep[4], pdf_mswep[4],   marker='o', markersize=3, mfc='violet', mec='violet', alpha=0.75, linestyle='None', label='MSWEP')
plt.plot(x_pdf_era5[4], pdf_era5[4],     marker='o', markersize=3, mfc='blue',   mec='blue',   alpha=0.75, linestyle='None', label='ERA5')
plt.title('(e) {0}'.format(dict_ws[4][2]), loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('{0}'.format(legend), fontsize=font_size, fontweight='bold')
plt.yscale('log')

ax = fig.add_subplot(3, 3, 9)
plt.plot(x_pdf_rcm3[5], pdf_rcm3[5],     marker='o', markersize=3, mfc='black',  mec='black',  alpha=0.75, linestyle='None', label='CPM3')
plt.plot(x_pdf_cpc[5], pdf_cpc[5],       marker='o', markersize=3, mfc='green',  mec='green',  alpha=0.75, linestyle='None', label='CPC')
plt.plot(x_pdf_cmorph[5], pdf_cmorph[5], marker='o', markersize=3, mfc='red',    mec='red',    alpha=0.75, linestyle='None', label='CMORPH')
plt.plot(x_pdf_mswep[5], pdf_mswep[5],   marker='o', markersize=3, mfc='violet', mec='violet', alpha=0.75, linestyle='None', label='MSWEP')
plt.plot(x_pdf_era5[5], pdf_era5[5],     marker='o', markersize=3, mfc='blue',   mec='blue',   alpha=0.75, linestyle='None', label='ERA5')
plt.title('(f) {0}'.format(dict_ws[5][2]), loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('{0}'.format(legend), fontsize=font_size, fontweight='bold')
plt.yscale('log')

# Path out to save figure
path_out = '{0}/user/mdasilva/CORDEX/figs'.format(path)
name_out = 'pyplt_graph_pdf_function_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
