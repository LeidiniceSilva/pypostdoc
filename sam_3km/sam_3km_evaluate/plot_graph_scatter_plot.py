# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot annual cycle"

import os
import netCDF4
import scipy.stats
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet

var = 'pr'
domain = 'SAM-3km'
idt, fdt = '2018', '2021'
dt = '{0}-{1}'.format(idt, fdt)
path = '/marconi/home/userexternal/mdasilva'

skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,7793,100,105,106,107,112,117,124,135,137,139,
149,152,155,158,168,174,177,183,186,199,204,210,212,224,226,239,240,248,249,253,254,276,277,280,293,298,
303,305,306,308,319,334,335,341,343,359,362,364,384,393,396,398,399,400,402,413,416,417,422,423,426,427,
443,444,446,451,453,457,458,467,474,479,483,488,489,490,495,505,509,513,514,516,529,534,544,559,566]
	
	
def import_situ(param_i, param_ii, domain):
	
	zz, mean_i, mean_ii = [], [], []
	for station in range(1, 567):
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue
		zz.append(inmet[station][4])

		arq_i  = xr.open_dataset('{0}/OBS/WS-SA/INMET/nc/hourly/{1}/'.format(path, param_i) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(param_i, inmet[station][0]))
		data_i = arq_i[param_i]
		time_i = data_i.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_i  = time_i.resample(time='1D').sum()
		var_i  = np.nanmean(var_i.values)
		mean_i.append(var_i)

		arq_ii    = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/'.format(path) + '{0}_{1}_RegCM5_day_{2}_lonlat.nc'.format(param_ii, domain, dt))
		data_ii   = arq_ii[param_ii]
		lonlat_ii = data_ii.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		time_ii   = lonlat_ii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_ii    = np.nanmean(time_ii.values)
		mean_ii.append(var_ii)
				
	return zz, mean_i, mean_ii
	
	
	
# Import model and obs dataset
dict_var = {'pr': ['pre', 'precip', 'sat_gauge_precip', 'tp']}
alt_i, inmet_i, regcm_i = import_situ(dict_var[var][0], var, domain)

#r, p = scipy.stats.pearsonr(inmet_i, regcm_i)

# Plot figure
fig = plt.figure()
font_size = 8

# Create scatter plot with colorbar based on the third variable
sc = plt.scatter(inmet_i, regcm_i, c=alt_i, cmap='terrain', marker='o', edgecolors='k')

#plt.text(2, 46, 'Corr = {0}'.format(r), fontsize=8)
#plt.text(2, 43, 'p_value = {0}'.format(p), fontsize=8)

# Add colorbar
cbar = plt.colorbar(sc)
cbar.set_label('Elevation (meters)', fontweight='bold', fontsize=8)

# Add diagonal line
plt.plot([0, 22], [0, 22], color='red', linestyle='--')

plt.title('a) SAM-3km', loc='left', fontweight='bold', fontsize=8)
plt.xlabel('INMET (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.ylabel('RegCM5 (mm d$^-$$^1$)', fontsize=font_size, fontweight='bold')
plt.xlim(0, 20)
plt.ylim(0, 20)
plt.xticks(np.arange(0, 22, 2), fontsize=font_size)
plt.yticks(np.arange(0, 22, 2), fontsize=font_size)
plt.grid(linestyle='--')

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_scatter_plot_{0}_{1}_RegCM5_2018-2021.png'.format(var, domain)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
