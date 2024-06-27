# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot pdf"

import os
import netCDF4
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i

var = 'pr'
freq= 'hourly'
domain = 'SAM-3km'
idt, fdt = '2018', '2021'
dt = '1hr_{0}-{1}'.format(idt, fdt)
legend = '(mm h$^-$$^1$)'

path = '/marconi/home/userexternal/mdasilva'
	
	
def import_situ_i():

	mean_i, mean_ii, mean_iii = [], [], []

	skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,7793,100,105,106,107,112,117,124,135,137,139,
	149,152,155,158,168,174,177,183,186,199,204,210,212,224,226,239,240,248,249,253,254,276,277,280,293,298,
	303,305,306,308,319,334,335,341,343,359,362,364,384,393,396,398,399,400,402,413,416,417,422,423,426,427,
	443,444,446,451,453,457,458,467,474,479,483,488,489,490,495,505,509,513,514,516,529,534,544,559,566]
	
	for station in range(1, 567): # 567
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
		var_i  = time_i.values
		mean_i.append(var_i)

		arq_ii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'tp_{0}_ERA5_{1}_lonlat.nc'.format(domain, dt))
		data_ii = arq_ii['tp']
		data_ii = data_ii.sel(lat=slice(yy-0.25,yy+0.25),lon=slice(xx-0.25,xx+0.25)).mean(('lat','lon'))
		time_ii = data_ii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_ii  = time_ii.values
		mean_ii.append(var_ii)

		arq_iii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/'.format(path) + 'pr_{0}_RegCM5_{1}_lonlat.nc'.format(domain, dt))
		data_iii = arq_iii['pr']
		data_iii = data_iii.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		time_iii = data_iii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_iii  = time_iii.values
		mean_iii.append(var_iii)
			
	return mean_i, mean_ii, mean_iii
	

def import_situ_ii():
	
	mean_i, mean_ii, mean_iii = [], [], []
	
	for station in range(1, 73): # 73
		print(station, smn_i[station][0])

		yy=smn_i[station][1]
		xx=smn_i[station][2]

		arq_i  = xr.open_dataset('{0}/OBS/WS-SA/SMN/hourly/nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[station][0]))
		data_i = arq_i['pre']
		time_i = data_i.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_i  = time_i.values
		mean_i.append(var_i)
		
		arq_ii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/obs/'.format(path) + 'tp_{0}_ERA5_{1}_lonlat.nc'.format(domain, dt))
		data_ii = arq_ii['tp']
		data_ii = data_ii.sel(lat=slice(yy-0.25,yy+0.25),lon=slice(xx-0.25,xx+0.25)).mean(('lat','lon'))
		time_ii = data_ii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_ii  = time_ii.values
		mean_ii.append(var_ii)

		arq_iii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/'.format(path) + 'pr_{0}_RegCM5_{1}_lonlat.nc'.format(domain, dt))
		data_iii = arq_iii['pr']
		data_iii = data_iii.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		time_iii = data_iii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_iii  = time_iii.values
		mean_iii.append(var_iii)
			
	return mean_i, mean_ii, mean_iii
	
	
# Import model and obs dataset
clim_i_x, clim_ii_x, clim_iii_x = import_situ_i()			
clim_i_y, clim_ii_y, clim_iii_y = import_situ_ii()			

inmet_smn = clim_i_x 
era5      = clim_ii_x 
regcm5    = clim_iii_x 

inmet_smn_ = clim_i_y 
era5_      = clim_ii_y
regcm5_    = clim_iii_y

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
		
inmet_smn_i, inmet_smn_ii, inmet_smn_iii, inmet_smn_iv, inmet_smn_v  = [], [], [], [], []
era5_i,      era5_ii,      era5_iii,      era5_iv,      era5_v       = [], [], [], [], []
regcm5_i,    regcm5_ii,    regcm5_iii,    regcm5_iv,    regcm5_v     = [], [], [], [], []

for c_i in count_i:
	inmet_smn_i.append(inmet_smn[c_i])
	era5_i.append(era5[c_i])
	regcm5_i.append(regcm5[c_i])

for c_ii in count_ii:
	inmet_smn_ii.append(inmet_smn[c_ii])
	era5_ii.append(era5[c_ii])
	regcm5_ii.append(regcm5[c_ii])
	
for c_iii in count_iii:
	inmet_smn_iii.append(inmet_smn[c_iii])
	era5_iii.append(era5[c_iii])
	regcm5_iii.append(regcm5[c_iii])
	
for c_iv in count_iv:
	inmet_smn_iv.append(inmet_smn[c_iv])
	era5_iv.append(era5[c_iv])
	regcm5_iv.append(regcm5[c_iv])
	
for c_v in count_v:
	inmet_smn_v.append(inmet_smn[c_v])
	era5_v.append(era5[c_v])
	regcm5_v.append(regcm5[c_v])

inmet_smn_array_i = np.array(inmet_smn_i)
inmet_smn_c_i     = inmet_smn_array_i.flatten()
era5_array_i      = np.array(era5_i)
era5_c_i          = era5_array_i.flatten()
regcm5_array_i    = np.array(regcm5_i)
regcm5_c_i        = regcm5_array_i.flatten()

inmet_smn_array_ii = np.array(inmet_smn_ii)
inmet_smn_c_ii     = inmet_smn_array_ii.flatten()
era5_array_ii      = np.array(era5_ii)
era5_c_ii          = era5_array_ii.flatten()
regcm5_array_ii    = np.array(regcm5_ii)
regcm5_c_ii        = regcm5_array_ii.flatten()

inmet_smn_array_iii = np.array(inmet_smn_iii)
inmet_smn_c_iii     = inmet_smn_array_iii.flatten()
era5_array_iii      = np.array(era5_iii)
era5_c_iii          = era5_array_iii.flatten()
regcm5_array_iii    = np.array(regcm5_iii)
regcm5_c_iii        = regcm5_array_iii.flatten()

inmet_smn_array_iv = np.array(inmet_smn_iv)
inmet_smn_c_iv     = inmet_smn_array_iv.flatten()
era5_array_iv      = np.array(era5_iv)
era5_c_iv          = era5_array_iv.flatten()
regcm5_array_iv    = np.array(regcm5_iv)
regcm5_c_iv        = regcm5_array_iv.flatten()

inmet_smn_c_ii_iii_iv = np.concatenate((inmet_smn_c_ii, inmet_smn_c_iii, inmet_smn_c_iv))
era5_c_ii_iii_iv      = np.concatenate((era5_c_ii, era5_c_iii, era5_c_iv))
regcm5_c_ii_iii_iv    = np.concatenate((regcm5_c_ii, regcm5_c_iii, regcm5_c_iv))

inmet_smn_array_v = np.array(inmet_smn_v)
inmet_smn_c_v     = inmet_smn_array_v.flatten()
era5_array_v      = np.array(era5_v)
era5_c_v          = era5_array_v.flatten()
regcm5_array_v    = np.array(regcm5_v)
regcm5_c_v        = regcm5_array_v.flatten()

inmet_smn_array_vi = np.array(inmet_smn_)
inmet_smn_c_vi     = inmet_smn_array_vi.flatten()
era5_array_vi      = np.array(era5_)
era5_c_vi          = era5_array_vi.flatten()
regcm5_array_vi    = np.array(regcm5_)
regcm5_c_vi        = regcm5_array_vi.flatten()

# Round values to each cluster
round_inmet_smn_c_i = np.round(inmet_smn_c_i,0)
round_era5_c_i      = np.round(era5_c_i,0)
round_regcm5_c_i    = np.round(regcm5_c_i,0)

round_inmet_smn_c_ii_iii_iv = np.round(inmet_smn_c_ii_iii_iv,0)
round_era5_c_ii_iii_iv      = np.round(era5_c_ii_iii_iv,0)
round_regcm5_c_ii_iii_iv    = np.round(regcm5_c_ii_iii_iv,0)

round_inmet_smn_c_v = np.round(inmet_smn_c_v,0)
round_era5_c_v      = np.round(era5_c_v,0)
round_regcm5_c_v    = np.round(regcm5_c_v,0)

round_inmet_smn_c_vi = np.round(inmet_smn_c_vi,0)
round_era5_c_vi      = np.round(era5_c_vi,0)
round_regcm5_c_vi    = np.round(regcm5_c_vi,0)

# Filter 0 mm/day
filter_inmet_smn_c_i = round_inmet_smn_c_i[round_inmet_smn_c_i > 0.]
filter_era5_c_i      = round_era5_c_i[round_era5_c_i > 0.]
filter_regcm5_c_i    = round_regcm5_c_i[round_regcm5_c_i > 0.]

filter_inmet_smn_c_ii_iii_iv = round_inmet_smn_c_ii_iii_iv[round_inmet_smn_c_ii_iii_iv > 0.]
filter_era5_c_ii_iii_iv      = round_era5_c_ii_iii_iv[round_era5_c_ii_iii_iv > 0.]
filter_regcm5_c_ii_iii_iv    = round_regcm5_c_ii_iii_iv[round_regcm5_c_ii_iii_iv > 0.]

filter_inmet_smn_c_v = round_inmet_smn_c_v[round_inmet_smn_c_v > 0.]
filter_era5_c_v      = round_era5_c_v[round_era5_c_v > 0.]
filter_regcm5_c_v    = round_regcm5_c_v[round_regcm5_c_v > 0.]

filter_inmet_smn_c_vi = round_inmet_smn_c_vi[round_inmet_smn_c_vi > 0.]
filter_era5_c_vi      = round_era5_c_vi[round_era5_c_vi > 0.]
filter_regcm5_c_vi    = round_regcm5_c_vi[round_regcm5_c_vi > 0.]

# Compute frequency
x_pdf_inmet_c_i, pdf_inmet_c_i = np.unique(filter_inmet_smn_c_i, return_counts=True) 
x_pdf_era5_c_i,  pdf_era5_c_i  = np.unique(filter_era5_c_i,      return_counts=True) 
x_pdf_regcm_c_i, pdf_regcm_c_i = np.unique(filter_regcm5_c_i,    return_counts=True) 

x_pdf_inmet_c_ii_iii_iv, pdf_inmet_c_ii_iii_iv = np.unique(filter_inmet_smn_c_ii_iii_iv, return_counts=True) 
x_pdf_era5_c_ii_iii_iv,  pdf_era5_c_ii_iii_iv  = np.unique(filter_era5_c_ii_iii_iv,      return_counts=True) 
x_pdf_regcm_c_ii_iii_iv, pdf_regcm_c_ii_iii_iv = np.unique(filter_regcm5_c_ii_iii_iv,    return_counts=True) 

x_pdf_inmet_c_v, pdf_inmet_c_v = np.unique(filter_inmet_smn_c_v, return_counts=True) 
x_pdf_era5_c_v,  pdf_era5_c_v  = np.unique(filter_era5_c_v,      return_counts=True) 
x_pdf_regcm_c_v, pdf_regcm_c_v = np.unique(filter_regcm5_c_v,    return_counts=True) 

x_pdf_inmet_c_vi, pdf_inmet_c_vi = np.unique(filter_inmet_smn_c_vi, return_counts=True) 
x_pdf_era5_c_vi,  pdf_era5_c_vi  = np.unique(filter_era5_c_vi,      return_counts=True) 
x_pdf_regcm_c_vi, pdf_regcm_c_vi = np.unique(filter_regcm5_c_vi,    return_counts=True)

# Plot figure
fig = plt.figure(figsize=(8, 6))
time = np.arange(0.5, 12 + 0.5)
font_size = 8

ax = fig.add_subplot(2, 2, 1)  
plt.plot(x_pdf_inmet_c_i, pdf_inmet_c_i, marker='.', markersize=4, mfc='blue',    mec='blue',    alpha=0.75, linestyle='None', label='INMET+SMN')
plt.plot(x_pdf_era5_c_i,  pdf_era5_c_i,  marker='.', markersize=4, mfc='green',   mec='green',   alpha=0.75, linestyle='None', label='ERA5')
plt.plot(x_pdf_regcm_c_i, pdf_regcm_c_i, marker='.', markersize=4, mfc='black',   mec='black',   alpha=0.75, linestyle='None', label='RegCM5')
plt.title('(a) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.legend(loc=1, ncol=2, fontsize=font_size, shadow=True)

ax = fig.add_subplot(2, 2, 2)  
plt.plot(x_pdf_inmet_c_ii_iii_iv, pdf_inmet_c_ii_iii_iv, marker='.', markersize=4, mfc='blue',    alpha=0.75, mec='blue',    linestyle='None', label='INMET+SMN')
plt.plot(x_pdf_era5_c_ii_iii_iv,  pdf_era5_c_ii_iii_iv,  marker='.', markersize=4, mfc='green',   alpha=0.75, mec='green',   linestyle='None', label='ERA5')
plt.plot(x_pdf_regcm_c_ii_iii_iv, pdf_regcm_c_ii_iii_iv, marker='.', markersize=4, mfc='black',   alpha=0.75, mec='black',   linestyle='None', label='RegCM5')
plt.title('(b) Cluster II-III-IV', loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')

ax = fig.add_subplot(2, 2, 3)  
plt.plot(x_pdf_inmet_c_v, pdf_inmet_c_v, marker='.', markersize=4, mfc='blue',    mec='blue',    alpha=0.75, linestyle='None', label='INMET+SMN')
plt.plot(x_pdf_era5_c_v,  pdf_era5_c_v,  marker='.', markersize=4, mfc='green',   mec='green',   alpha=0.75, linestyle='None', label='ERA5')
plt.plot(x_pdf_regcm_c_v, pdf_regcm_c_v, marker='.', markersize=4, mfc='black',   mec='black',   alpha=0.75, linestyle='None', label='RegCM5')
plt.title('(c) Cluster V', loc='left', fontsize=font_size, fontweight='bold') 
plt.yscale('log')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.xlabel('{0}'.format(legend), fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 2, 4)  
plt.plot(x_pdf_inmet_c_vi, pdf_inmet_c_vi, marker='.', markersize=4, mfc='blue',    mec='blue',    alpha=0.75, linestyle='None', label='INMET+SMN')
plt.plot(x_pdf_era5_c_vi,  pdf_era5_c_vi,  marker='.', markersize=4, mfc='green',   mec='green',   alpha=0.75, linestyle='None', label='ERA5')
plt.plot(x_pdf_regcm_c_vi, pdf_regcm_c_vi, marker='.', markersize=4, mfc='black',   mec='black',   alpha=0.75, linestyle='None', label='RegCM5')
plt.title('(d) Cluster VI', loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.xlabel('{0}'.format(legend), fontsize=font_size, fontweight='bold')

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_pdf_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
