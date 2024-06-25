# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 15, 2024"
__description__ = "This script plot pdf function"

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from dict_inmet_stations_dc import inmet
from dict_smn_stations_dc import smn

path = '/marconi/home/userexternal/mdasilva'

var = 'pr'
freq = '1hr'

if freq == '1hr':
	dt = '1hr_20180101-20211231'
	legend = '(mm h$^-$$^1$)'
else:
	dt = 'day_20180101-20211231'
	legend = '(mm d$^-$$^1$)'


def import_inmet():

	mean_     = []
	mean_i    = []
	mean_ii   = []
	mean_iii  = []
	mean_iv   = []
	mean_v    = []
	mean_vi   = []
	mean_vii  = []

	# Select lat and lon 
	for i in range(1, 100):
		yy=inmet[i][2]
		xx=inmet[i][3]
				
		print('Reading weather station:', i, inmet[i][0])	
		# reading regcm ictp		
		d_ = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/'.format(path) + 'pr_SESA_RegCM5_1hr_2018-2021_lonlat.nc')
		d_ = d_.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ = d_.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		d_ = d_.values
		mean_.append(d_)
		
		# reading regcm usp 
		d_i = xr.open_dataset('{0}/user/mdasilva/CSAM-4i/RegCM4/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_USP-RegCM471_v2_{0}.nc'.format(dt))
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.values
		mean_i.append(d_i)
			
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/user/mdasilva/CSAM-4i/RegCM5/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_{0}.nc'.format(dt))
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.values
		mean_ii.append(d_ii)
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/user/mdasilva/CSAM-4i/RegCM5/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_{0}.nc'.format(dt))
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.values
		mean_iii.append(d_iii)
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/user/mdasilva/CSAM-4i/WRF-NCAR/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_NCAR-WRF415_v1_{0}.nc'.format(dt))
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.values
		mean_iv.append(d_iv)
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/user/mdasilva/CSAM-4i/WRF-UCAN/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_{0}.nc'.format(dt))
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.values
		mean_v.append(d_v)
		
		# Reading inmet 
		d_vi = xr.open_dataset('{0}/OBS/BDMET/database/nc/hourly/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.values
		mean_vi.append(d_vi)

		# Reading gpm 
		d_vii = xr.open_dataset('{0}/OBS/GPM/'.format(path) + 'precipitation_SESA_1hr_2018-2021.nc'.format(inmet[i][0]))
		d_vii = d_vii.precipitation.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=yy, lon=xx, method='nearest')
		d_vii = d_vii.values
		mean_vii.append(d_vii)
				
	return mean_, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii


def import_smn():
	
	mean_     = []
	mean_i    = []
	mean_ii   = []
	mean_iii  = []
	mean_iv   = []
	mean_v    = []
	mean_vi   = []
	mean_vii  = []
	
	# Select lat and lon 
	for i in range(1, 73):
		yy=smn[i][1]
		xx=smn[i][2]
		
		print('Reading weather station:', i, smn[i][0])	
		# reading regcm ictp		
		d_ = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/'.format(path) + 'pr_SESA_RegCM5_1hr_2018-2021_lonlat.nc')
		d_ = d_.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ = d_.sel(lat=slice(yy-0.03,yy+0.03),lon=slice(xx-0.03,xx+0.03)).mean(('lat','lon'))
		d_ = d_.values
		mean_.append(d_)
		
		# reading regcm usp 
		d_i = xr.open_dataset('{0}/user/mdasilva/CSAM-4i/RegCM4/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_USP-RegCM471_v2_{0}.nc'.format(dt))
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.values
		mean_i.append(d_i)
			
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/user/mdasilva/CSAM-4i/RegCM5/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_{0}.nc'.format(dt))
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.values
		mean_ii.append(d_ii)
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/user/mdasilva/CSAM-4i/RegCM5/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_{0}.nc'.format(dt))
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.values
		mean_iii.append(d_iii)
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/user/mdasilva/CSAM-4i/WRF-NCAR/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_NCAR-WRF415_v1_{0}.nc'.format(dt))
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.values
		mean_iv.append(d_iv)
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/user/mdasilva/CSAM-4i/WRF-UCAN/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_{0}.nc'.format(dt))
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.values
		mean_v.append(d_v)
						
		# Reading smn 
		d_vi = xr.open_dataset('{0}/OBS/SMN/smn_nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.values
		mean_vi.append(d_vi)
		
		# Reading gpm 
		d_vii = xr.open_dataset('{0}/OBS/GPM/'.format(path) + 'precipitation_SESA_1hr_2018-2021.nc'.format(inmet[i][0]))
		d_vii = d_vii.precipitation.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=yy, lon=xx, method='nearest')
		d_vii = d_vii.values
		mean_vii.append(d_vii)
				
	return mean_, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii
	
	
# Import dataset
clim_x, clim_i_x, clim_ii_x, clim_iii_x, clim_iv_x, clim_v_x, clim_vi_x, clim_vii_x = import_inmet()			
clim_y, clim_i_y, clim_ii_y, clim_iii_y, clim_iv_y, clim_v_y, clim_vi_y, clim_vii_y = import_smn()
			
reg_usp      = clim_i_x   + clim_i_y 
reg_ictp_i   = clim_ii_x  + clim_ii_y 
reg_ictp_ii  = clim_iii_x + clim_iii_y 
reg_ictp_iii = clim_x     + clim_y 
wrf_ncar     = clim_iv_x  + clim_iv_y  
wrf_ucan     = clim_v_x   + clim_v_y  
inmet_smn    = clim_vi_x  + clim_vi_y  
gpm          = clim_vii_x + clim_vii_y  

list_hc = [3, 3, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 3, 3, 3, 3, 4, 3, 3, 3, 1, 1, 3, 1, 4, 3,
3, 3, 3, 1, 4, 3, 1, 3, 3, 4, 3, 3, 3, 1, 4, 3, 1, 1, 3, 1, 4, 3, 4, 1, 3, 3, 3, 4, 1, 2, 4, 1, 1, 1, 1, 1, 
3, 2, 1, 1, 2, 3, 1, 3, 2, 1, 2, 1, 2, 1, 2, 2, 2, 2, 2, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 2, 0, 0, 
0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0]

print(len(reg_usp))
print(len(reg_ictp_i))
print(len(reg_ictp_ii))
print(len(reg_ictp_iii))
print(len(wrf_ncar))
print(len(wrf_ucan))
print(len(inmet_smn))
print(len(gpm))
print(len(list_hc))

count_i = []
count_ii = []
count_iii = []
count_iv = []
count_v = []

# Select cluster
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

reg_usp_i,      reg_usp_ii,      reg_usp_iii,      reg_usp_iv,      reg_usp_v      = [], [], [], [], []
reg_ictp_i_i,   reg_ictp_i_ii,   reg_ictp_i_iii,   reg_ictp_i_iv,   reg_ictp_i_v   = [], [], [], [], []
reg_ictp_ii_i,  reg_ictp_ii_ii , reg_ictp_ii_iii,  reg_ictp_ii_iv,  reg_ictp_ii_v  = [], [], [], [], []
reg_ictp_iii_i, reg_ictp_iii_ii, reg_ictp_iii_iii, reg_ictp_iii_iv, reg_ictp_iii_v = [], [], [], [], []
wrf_ncar_i,     wrf_ncar_ii,     wrf_ncar_iii,     wrf_ncar_iv,     wrf_ncar_v     = [], [], [], [], []
wrf_ucan_i,     wrf_ucan_ii,     wrf_ucan_iii,     wrf_ucan_iv,     wrf_ucan_v     = [], [], [], [], []
inmet_smn_i,    inmet_smn_ii,    inmet_smn_iii,    inmet_smn_iv,    inmet_smn_v    = [], [], [], [], []
gpm_i,          gpm_ii,          gpm_iii,          gpm_iv,          gpm_v          = [], [], [], [], []

for c_i in count_i:
	reg_usp_i.append(reg_usp[c_i])
	reg_ictp_i_i.append(reg_ictp_i[c_i])
	reg_ictp_ii_i.append(reg_ictp_ii[c_i])
	reg_ictp_iii_i.append(reg_ictp_iii[c_i])
	wrf_ncar_i.append(wrf_ncar[c_i])
	wrf_ucan_i.append(wrf_ucan[c_i])
	inmet_smn_i.append(inmet_smn[c_i])
	gpm_i.append(gpm[c_i])

for c_ii in count_ii:
	reg_usp_ii.append(reg_usp[c_ii])
	reg_ictp_i_ii.append(reg_ictp_i[c_ii])
	reg_ictp_ii_ii.append(reg_ictp_ii[c_ii])
	reg_ictp_iii_ii.append(reg_ictp_iii[c_ii])
	wrf_ncar_ii.append(wrf_ncar[c_ii])
	wrf_ucan_ii.append(wrf_ucan[c_ii])
	inmet_smn_ii.append(inmet_smn[c_ii])
	gpm_ii.append(gpm[c_ii])
	
for c_iii in count_iii:
	reg_usp_iii.append(reg_usp[c_iii])
	reg_ictp_i_iii.append(reg_ictp_i[c_iii])
	reg_ictp_ii_iii.append(reg_ictp_ii[c_iii])
	reg_ictp_iii_iii.append(reg_ictp_iii[c_iii])
	wrf_ncar_iii.append(wrf_ncar[c_iii])
	wrf_ucan_iii.append(wrf_ucan[c_iii])
	inmet_smn_iii.append(inmet_smn[c_iii])
	gpm_iii.append(gpm[c_iii])
	
for c_iv in count_iv:
	reg_usp_iv.append(reg_usp[c_iv])
	reg_ictp_i_iv.append(reg_ictp_i[c_iv])
	reg_ictp_ii_iv.append(reg_ictp_ii[c_iv])
	reg_ictp_iii_iv.append(reg_ictp_iii[c_iv])
	wrf_ncar_iv.append(wrf_ncar[c_iv])
	wrf_ucan_iv.append(wrf_ucan[c_iv])
	inmet_smn_iv.append(inmet_smn[c_iv])
	gpm_iv.append(gpm[c_iv])
	
for c_v in count_v:
	reg_usp_v.append(reg_usp[c_v])
	reg_ictp_i_v.append(reg_ictp_i[c_v])
	reg_ictp_ii_v.append(reg_ictp_ii[c_v])
	reg_ictp_iii_v.append(reg_ictp_iii[c_v])
	wrf_ncar_v.append(wrf_ncar[c_v])
	wrf_ucan_v.append(wrf_ucan[c_v])
	inmet_smn_v.append(inmet_smn[c_v])
	gpm_v.append(gpm[c_v])

reg_usp_i   = np.array(reg_usp_i)
reg_usp_c_i = reg_usp_i.flatten()

reg_ictp_i_i   = np.array(reg_ictp_i_i)
reg_ictp_i_c_i = reg_ictp_i_i.flatten()

reg_ictp_ii_i   = np.array(reg_ictp_ii_i)
reg_ictp_ii_c_i = reg_ictp_ii_i.flatten()

reg_ictp_iii_i   = np.array(reg_ictp_iii_i)
reg_ictp_iii_c_i = reg_ictp_iii_i.flatten()

wrf_ncar_i   = np.array(wrf_ncar_i)
wrf_ncar_c_i = wrf_ncar_i.flatten()

wrf_ucan_i   = np.array(wrf_ucan_i)
wrf_ucan_c_i = wrf_ucan_i.flatten()

inmet_smn_i   = np.array(inmet_smn_i)
inmet_smn_c_i = inmet_smn_i.flatten()

gpm_i   = np.array(gpm_i)
gpm_c_i = gpm_i.flatten()



reg_usp_ii   = np.array(reg_usp_ii)
reg_usp_c_ii = reg_usp_ii.flatten()

reg_ictp_i_ii   = np.array(reg_ictp_i_ii)
reg_ictp_i_c_ii = reg_ictp_i_ii.flatten()

reg_ictp_ii_ii   = np.array(reg_ictp_ii_ii)
reg_ictp_ii_c_ii = reg_ictp_ii_ii.flatten()

reg_ictp_iii_ii   = np.array(reg_ictp_iii_ii)
reg_ictp_iii_c_ii = reg_ictp_iii_ii.flatten()

wrf_ncar_ii   = np.array(wrf_ncar_ii)
wrf_ncar_c_ii = wrf_ncar_ii.flatten()

wrf_ucan_ii   = np.array(wrf_ucan_ii)
wrf_ucan_c_ii = wrf_ucan_ii.flatten()

inmet_smn_ii   = np.array(inmet_smn_ii)
inmet_smn_c_ii = inmet_smn_ii.flatten()

gpm_ii   = np.array(gpm_ii)
gpm_c_ii = gpm_ii.flatten()



reg_usp_iii   = np.array(reg_usp_iii)
reg_usp_c_iii = reg_usp_iii.flatten()

reg_ictp_i_iii   = np.array(reg_ictp_i_iii)
reg_ictp_i_c_iii = reg_ictp_i_iii.flatten()

reg_ictp_ii_iii   = np.array(reg_ictp_ii_iii)
reg_ictp_ii_c_iii = reg_ictp_ii_iii.flatten()

reg_ictp_iii_iii   = np.array(reg_ictp_iii_iii)
reg_ictp_iii_c_iii = reg_ictp_iii_iii.flatten()

wrf_ncar_iii   = np.array(wrf_ncar_iii)
wrf_ncar_c_iii = wrf_ncar_iii.flatten()

wrf_ucan_iii   = np.array(wrf_ucan_iii)
wrf_ucan_c_iii = wrf_ucan_iii.flatten()

inmet_smn_iii   = np.array(inmet_smn_iii)
inmet_smn_c_iii = inmet_smn_iii.flatten()

gpm_iii   = np.array(gpm_iii)
gpm_c_iii = gpm_iii.flatten()



reg_usp_iv   = np.array(reg_usp_iv)
reg_usp_c_iv = reg_usp_iv.flatten()

reg_ictp_i_iv   = np.array(reg_ictp_i_iv)
reg_ictp_i_c_iv = reg_ictp_i_iv.flatten()

reg_ictp_ii_iv   = np.array(reg_ictp_ii_iv)
reg_ictp_ii_c_iv = reg_ictp_ii_iv.flatten()

reg_ictp_iii_iv   = np.array(reg_ictp_iii_iv)
reg_ictp_iii_c_iv = reg_ictp_iii_iv.flatten()

wrf_ncar_iv   = np.array(wrf_ncar_iv)
wrf_ncar_c_iv = wrf_ncar_iv.flatten()

wrf_ucan_iv   = np.array(wrf_ucan_iv)
wrf_ucan_c_iv = wrf_ucan_iv.flatten()

inmet_smn_iv   = np.array(inmet_smn_iv)
inmet_smn_c_iv = inmet_smn_iv.flatten()

gpm_iv   = np.array(gpm_iv)
gpm_c_iv = gpm_iv.flatten()



reg_usp_v   = np.array(reg_usp_v)
reg_usp_c_v = reg_usp_v.flatten()

reg_ictp_i_v   = np.array(reg_ictp_i_v)
reg_ictp_i_c_v = reg_ictp_i_v.flatten()

reg_ictp_ii_v   = np.array(reg_ictp_ii_v)
reg_ictp_ii_c_v = reg_ictp_ii_v.flatten()

reg_ictp_iii_v   = np.array(reg_ictp_iii_v)
reg_ictp_iii_c_v = reg_ictp_iii_v.flatten()

wrf_ncar_v   = np.array(wrf_ncar_v)
wrf_ncar_c_v = wrf_ncar_v.flatten()

wrf_ucan_v   = np.array(wrf_ucan_v)
wrf_ucan_c_v = wrf_ucan_v.flatten()

inmet_smn_v   = np.array(inmet_smn_v)
inmet_smn_c_v = inmet_smn_v.flatten()

gpm_v   = np.array(gpm_v)
gpm_c_v = gpm_v.flatten()

# Round values to each cluster
round_reg_usp_c_i      = np.round(reg_usp_c_i,0)
round_reg_ictp_i_c_i   = np.round(reg_ictp_i_c_i,0)
round_reg_ictp_ii_c_i  = np.round(reg_ictp_ii_c_i,0)
round_reg_ictp_iii_c_i = np.round(reg_ictp_iii_c_i,0)
round_wrf_ncar_c_i     = np.round(wrf_ncar_c_i,0)
round_wrf_ucan_c_i     = np.round(wrf_ucan_c_i,0)
round_inmet_smn_c_i    = np.round(inmet_smn_c_i,0)
round_gpm_c_i          = np.round(gpm_c_i,0)

round_reg_usp_c_ii      = np.round(reg_usp_c_ii,0)
round_reg_ictp_i_c_ii   = np.round(reg_ictp_i_c_ii,0)
round_reg_ictp_ii_c_ii  = np.round(reg_ictp_ii_c_ii,0)
round_reg_ictp_iii_c_ii = np.round(reg_ictp_iii_c_ii,0)
round_wrf_ncar_c_ii     = np.round(wrf_ncar_c_ii,0)
round_wrf_ucan_c_ii     = np.round(wrf_ucan_c_ii,0)
round_inmet_smn_c_ii    = np.round(inmet_smn_c_ii,0)
round_gpm_c_ii          = np.round(gpm_c_ii,0)

round_reg_usp_c_iii      = np.round(reg_usp_c_iii,0)
round_reg_ictp_i_c_iii   = np.round(reg_ictp_i_c_iii,0)
round_reg_ictp_ii_c_iii  = np.round(reg_ictp_ii_c_iii,0)
round_reg_ictp_iii_c_iii = np.round(reg_ictp_iii_c_iii,0)
round_wrf_ncar_c_iii     = np.round(wrf_ncar_c_iii,0)
round_wrf_ucan_c_iii     = np.round(wrf_ucan_c_iii,0)
round_inmet_smn_c_iii    = np.round(inmet_smn_c_iii,0)
round_gpm_c_iii          = np.round(gpm_c_iii,0)

round_reg_usp_c_iv      = np.round(reg_usp_c_iv,0)
round_reg_ictp_i_c_iv   = np.round(reg_ictp_i_c_iv,0)
round_reg_ictp_ii_c_iv  = np.round(reg_ictp_ii_c_iv,0)
round_reg_ictp_iii_c_iv = np.round(reg_ictp_iii_c_iv,0)
round_wrf_ncar_c_iv     = np.round(wrf_ncar_c_iv,0)
round_wrf_ucan_c_iv     = np.round(wrf_ucan_c_iv,0)
round_inmet_smn_c_iv    = np.round(inmet_smn_c_iv,0)
round_gpm_c_iv          = np.round(gpm_c_iv,0)

round_reg_usp_c_v      = np.round(reg_usp_c_v,0)
round_reg_ictp_i_c_v   = np.round(reg_ictp_i_c_v,0)
round_reg_ictp_ii_c_v  = np.round(reg_ictp_ii_c_v,0)
round_reg_ictp_iii_c_v = np.round(reg_ictp_iii_c_v,0)
round_wrf_ncar_c_v     = np.round(wrf_ncar_c_v,0)
round_wrf_ucan_c_v     = np.round(wrf_ucan_c_v,0)
round_inmet_smn_c_v    = np.round(inmet_smn_c_v,0)
round_gpm_c_v          = np.round(gpm_c_v,0)

# Filter 0 mm/day
filter_reg_usp_c_i      = round_reg_usp_c_i[round_reg_usp_c_i > 0.]
filter_reg_ictp_i_c_i   = round_reg_ictp_i_c_i[round_reg_ictp_i_c_i > 0.]
filter_reg_ictp_ii_c_i  = round_reg_ictp_ii_c_i[round_reg_ictp_ii_c_i > 0.]
filter_reg_ictp_iii_c_i = round_reg_ictp_iii_c_i[round_reg_ictp_iii_c_i > 0.]
filter_wrf_ncar_c_i     = round_wrf_ncar_c_i[round_wrf_ncar_c_i > 0.]
filter_wrf_ucan_c_i     = round_wrf_ucan_c_i[round_wrf_ucan_c_i > 0.]
filter_inmet_smn_c_i    = round_inmet_smn_c_i[round_inmet_smn_c_i > 0.]
filter_gpm_c_i          = round_gpm_c_i[round_gpm_c_i > 0.]

filter_reg_usp_c_ii      = round_reg_usp_c_ii[round_reg_usp_c_ii > 0.]
filter_reg_ictp_i_c_ii   = round_reg_ictp_i_c_ii[round_reg_ictp_i_c_ii > 0.]
filter_reg_ictp_ii_c_ii  = round_reg_ictp_ii_c_ii[round_reg_ictp_ii_c_ii > 0.]
filter_reg_ictp_iii_c_ii = round_reg_ictp_iii_c_ii[round_reg_ictp_iii_c_ii > 0.]
filter_wrf_ncar_c_ii     = round_wrf_ncar_c_ii[round_wrf_ncar_c_ii > 0.]
filter_wrf_ucan_c_ii     = round_wrf_ucan_c_ii[round_wrf_ucan_c_ii > 0.]
filter_inmet_smn_c_ii    = round_inmet_smn_c_ii[round_inmet_smn_c_ii > 0.]
filter_gpm_c_ii          = round_gpm_c_ii[round_gpm_c_ii > 0.]

filter_reg_usp_c_iii      = round_reg_usp_c_iii[round_reg_usp_c_iii > 0.]
filter_reg_ictp_i_c_iii   = round_reg_ictp_i_c_iii[round_reg_ictp_i_c_iii > 0.]
filter_reg_ictp_ii_c_iii  = round_reg_ictp_ii_c_iii[round_reg_ictp_ii_c_iii > 0.]
filter_reg_ictp_iii_c_iii = round_reg_ictp_iii_c_iii[round_reg_ictp_iii_c_iii > 0.]
filter_wrf_ncar_c_iii     = round_wrf_ncar_c_iii[round_wrf_ncar_c_iii > 0.]
filter_wrf_ucan_c_iii     = round_wrf_ucan_c_iii[round_wrf_ucan_c_iii > 0.]
filter_inmet_smn_c_iii    = round_inmet_smn_c_iii[round_inmet_smn_c_iii > 0.]
filter_gpm_c_iii          = round_gpm_c_iii[round_gpm_c_iii > 0.]
 
filter_reg_usp_c_iv      = round_reg_usp_c_iv[round_reg_usp_c_iv > 0.]
filter_reg_ictp_i_c_iv   = round_reg_ictp_i_c_iv[round_reg_ictp_i_c_iv > 0.]
filter_reg_ictp_ii_c_iv  = round_reg_ictp_ii_c_iv[round_reg_ictp_ii_c_iv > 0.]
filter_reg_ictp_iii_c_iv = round_reg_ictp_iii_c_iv[round_reg_ictp_iii_c_iv > 0.]
filter_wrf_ncar_c_iv     = round_wrf_ncar_c_iv[round_wrf_ncar_c_iv > 0.]
filter_wrf_ucan_c_iv     = round_wrf_ucan_c_iv[round_wrf_ucan_c_iv > 0.]
filter_inmet_smn_c_iv    = round_inmet_smn_c_iv[round_inmet_smn_c_iv > 0.]
filter_gpm_c_iv          = round_gpm_c_iv[round_gpm_c_iv > 0.]

filter_reg_usp_c_v      = round_reg_usp_c_v[round_reg_usp_c_v > 0.]
filter_reg_ictp_i_c_v   = round_reg_ictp_i_c_v[round_reg_ictp_i_c_v > 0.]
filter_reg_ictp_ii_c_v  = round_reg_ictp_ii_c_v[round_reg_ictp_ii_c_v > 0.]
filter_reg_ictp_iii_c_v = round_reg_ictp_iii_c_v[round_reg_ictp_iii_c_v > 0.]
filter_wrf_ncar_c_v     = round_wrf_ncar_c_v[round_wrf_ncar_c_v > 0.]
filter_wrf_ucan_c_v     = round_wrf_ucan_c_v[round_wrf_ucan_c_v > 0.]
filter_inmet_smn_c_v    = round_inmet_smn_c_v[round_inmet_smn_c_v > 0.]
filter_gpm_c_v          = round_gpm_c_v[round_gpm_c_v > 0.]

# Compute frequency
x_freq_reg_usp_c_i,      freq_reg_usp_c_i      = np.unique(filter_reg_usp_c_i,      return_counts=True) 
x_freq_reg_ictp_i_c_i,   freq_reg_ictp_i_c_i   = np.unique(filter_reg_ictp_i_c_i,   return_counts=True) 
x_freq_reg_ictp_ii_c_i,  freq_reg_ictp_ii_c_i  = np.unique(filter_reg_ictp_ii_c_i,  return_counts=True) 
x_freq_reg_ictp_iii_c_i, freq_reg_ictp_iii_c_i = np.unique(filter_reg_ictp_iii_c_i, return_counts=True) 
x_freq_wrf_ncar_c_i,    freq_wrf_ncar_c_i      = np.unique(filter_wrf_ncar_c_i,     return_counts=True) 
x_freq_wrf_ucan_c_i,    freq_wrf_ucan_c_i      = np.unique(filter_wrf_ucan_c_i,     return_counts=True) 
x_freq_inmet_smn_c_i,   freq_inmet_smn_c_i     = np.unique(filter_inmet_smn_c_i,    return_counts=True) 
x_freq_gpm_c_i,         freq_gpm_c_i           = np.unique(filter_gpm_c_i,          return_counts=True) 

x_freq_reg_usp_c_ii,      freq_reg_usp_c_ii      = np.unique(filter_reg_usp_c_ii,      return_counts=True) 
x_freq_reg_ictp_i_c_ii,   freq_reg_ictp_i_c_ii   = np.unique(filter_reg_ictp_i_c_ii,   return_counts=True) 
x_freq_reg_ictp_ii_c_ii,  freq_reg_ictp_ii_c_ii  = np.unique(filter_reg_ictp_ii_c_ii,  return_counts=True) 
x_freq_reg_ictp_iii_c_ii, freq_reg_ictp_iii_c_ii = np.unique(filter_reg_ictp_iii_c_ii, return_counts=True) 
x_freq_wrf_ncar_c_ii,     freq_wrf_ncar_c_ii     = np.unique(filter_wrf_ncar_c_ii,     return_counts=True) 
x_freq_wrf_ucan_c_ii,     freq_wrf_ucan_c_ii     = np.unique(filter_wrf_ucan_c_ii,     return_counts=True) 
x_freq_inmet_smn_c_ii,    freq_inmet_smn_c_ii    = np.unique(filter_inmet_smn_c_ii,    return_counts=True) 
x_freq_gpm_c_ii,          freq_gpm_c_ii          = np.unique(filter_gpm_c_ii,          return_counts=True) 

x_freq_reg_usp_c_iii,      freq_reg_usp_c_iii      = np.unique(filter_reg_usp_c_iii,      return_counts=True) 
x_freq_reg_ictp_i_c_iii,   freq_reg_ictp_i_c_iii   = np.unique(filter_reg_ictp_i_c_iii,   return_counts=True) 
x_freq_reg_ictp_ii_c_iii,  freq_reg_ictp_ii_c_iii  = np.unique(filter_reg_ictp_ii_c_iii,  return_counts=True) 
x_freq_reg_ictp_iii_c_iii, freq_reg_ictp_iii_c_iii = np.unique(filter_reg_ictp_iii_c_iii, return_counts=True) 
x_freq_wrf_ncar_c_iii,     freq_wrf_ncar_c_iii     = np.unique(filter_wrf_ncar_c_iii,     return_counts=True) 
x_freq_wrf_ucan_c_iii,     freq_wrf_ucan_c_iii     = np.unique(filter_wrf_ucan_c_iii,     return_counts=True) 
x_freq_inmet_smn_c_iii,    freq_inmet_smn_c_iii    = np.unique(filter_inmet_smn_c_iii,    return_counts=True) 
x_freq_gpm_c_iii,          freq_gpm_c_iii          = np.unique(filter_gpm_c_iii,          return_counts=True) 

x_freq_reg_usp_c_iv,      freq_reg_usp_c_iv      = np.unique(filter_reg_usp_c_iv,      return_counts=True) 
x_freq_reg_ictp_i_c_iv,   freq_reg_ictp_i_c_iv   = np.unique(filter_reg_ictp_i_c_iv,   return_counts=True) 
x_freq_reg_ictp_ii_c_iv,  freq_reg_ictp_ii_c_iv  = np.unique(filter_reg_ictp_ii_c_iv,  return_counts=True) 
x_freq_reg_ictp_iii_c_iv, freq_reg_ictp_iii_c_iv = np.unique(filter_reg_ictp_iii_c_iv, return_counts=True) 
x_freq_wrf_ncar_c_iv,     freq_wrf_ncar_c_iv     = np.unique(filter_wrf_ncar_c_iv,     return_counts=True) 
x_freq_wrf_ucan_c_iv,     freq_wrf_ucan_c_iv     = np.unique(filter_wrf_ucan_c_iv,     return_counts=True) 
x_freq_inmet_smn_c_iv,    freq_inmet_smn_c_iv    = np.unique(filter_inmet_smn_c_iv,    return_counts=True) 
x_freq_gpm_c_iv,          freq_gpm_c_iv          = np.unique(filter_gpm_c_iv,          return_counts=True) 

x_freq_reg_usp_c_v,      freq_reg_usp_c_v      = np.unique(filter_reg_usp_c_v,      return_counts=True) 
x_freq_reg_ictp_i_c_v,   freq_reg_ictp_i_c_v   = np.unique(filter_reg_ictp_i_c_v,   return_counts=True) 
x_freq_reg_ictp_ii_c_v,  freq_reg_ictp_ii_c_v  = np.unique(filter_reg_ictp_ii_c_v,  return_counts=True) 
x_freq_reg_ictp_iii_c_v, freq_reg_ictp_iii_c_v = np.unique(filter_reg_ictp_iii_c_v, return_counts=True) 
x_freq_wrf_ncar_c_v,     freq_wrf_ncar_c_v     = np.unique(filter_wrf_ncar_c_v,     return_counts=True) 
x_freq_wrf_ucan_c_v,     freq_wrf_ucan_c_v     = np.unique(filter_wrf_ucan_c_v,     return_counts=True) 
x_freq_inmet_smn_c_v,    freq_inmet_smn_c_v    = np.unique(filter_inmet_smn_c_v,    return_counts=True) 
x_freq_gpm_c_v,          freq_gpm_c_v          = np.unique(filter_gpm_c_v,          return_counts=True) 

# Plot figure
fig = plt.figure(figsize=(10, 10))
font_size = 10

ax = fig.add_subplot(3, 2, 1)  
plt.plot(x_freq_gpm_c_i,          freq_gpm_c_i,          marker='.', markersize=4, mfc='gray',   mec='gray',   linestyle='None', label='GPMsat')
plt.plot(x_freq_inmet_smn_c_i,    freq_inmet_smn_c_i,    marker='.', markersize=4, mfc='black',  mec='black',  linestyle='None', label='station')
plt.plot(x_freq_reg_ictp_i_c_i,   freq_reg_ictp_i_c_i,   marker='.', markersize=4, mfc='blue',   mec='blue',   linestyle='None', label='RegCM5-pbl1')
plt.plot(x_freq_reg_ictp_ii_c_i,  freq_reg_ictp_ii_c_i,  marker='.', markersize=4, mfc='cyan',   mec='cyan',   linestyle='None', label='RegCM5-pbl2')
plt.plot(x_freq_reg_ictp_iii_c_i, freq_reg_ictp_iii_c_i, marker='.', markersize=4, mfc='violet', mec='violet', linestyle='None', label='RegCM5-3km')
plt.plot(x_freq_reg_usp_c_i,      freq_reg_usp_c_i,      marker='.', markersize=4, mfc='green',  mec='green',  linestyle='None', label='USP-RegCM471')
plt.plot(x_freq_wrf_ncar_c_i,     freq_wrf_ncar_c_i,     marker='.', markersize=4, mfc='red',    mec='red',    linestyle='None', label='NCAR-WRF413')
plt.plot(x_freq_wrf_ucan_c_i,     freq_wrf_ucan_c_i,     marker='.', markersize=4, mfc='orange', mec='orange', linestyle='None', label='UCAN-WRF433')
plt.title('(a) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.legend(loc=1, ncol=1, frameon=False)

ax = fig.add_subplot(3, 2, 2)
plt.plot(x_freq_gpm_c_ii,          freq_gpm_c_ii,          marker='.', markersize=4, mfc='gray',   mec='gray',   linestyle='None', label='GPMsat')
plt.plot(x_freq_inmet_smn_c_ii,    freq_inmet_smn_c_ii,    marker='.', markersize=4, mfc='black',  mec='black',  linestyle='None', label='station')
plt.plot(x_freq_reg_ictp_i_c_ii,   freq_reg_ictp_i_c_ii,   marker='.', markersize=4, mfc='blue',   mec='blue',   linestyle='None', label='RegCM5-pbl1')
plt.plot(x_freq_reg_ictp_ii_c_ii,  freq_reg_ictp_ii_c_ii,  marker='.', markersize=4, mfc='cyan',   mec='cyan',   linestyle='None', label='RegCM5-pbl2')
plt.plot(x_freq_reg_ictp_iii_c_ii, freq_reg_ictp_iii_c_ii, marker='.', markersize=4, mfc='violet', mec='violet', linestyle='None', label='RegCM5-3km')
plt.plot(x_freq_reg_usp_c_ii,      freq_reg_usp_c_ii,      marker='.', markersize=4, mfc='green',  mec='green',  linestyle='None', label='USP-RegCM471')
plt.plot(x_freq_wrf_ncar_c_ii,     freq_wrf_ncar_c_ii,     marker='.', markersize=4, mfc='red',    mec='red',    linestyle='None', label='NCAR-WRF413')
plt.plot(x_freq_wrf_ucan_c_ii,     freq_wrf_ucan_c_ii,     marker='.', markersize=4, mfc='orange', mec='orange', linestyle='None', label='UCAN-WRF433')
plt.title('(b) Cluster II', loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')

ax = fig.add_subplot(3, 2, 3)
plt.plot(x_freq_gpm_c_iii,          freq_gpm_c_iii,          marker='.', markersize=4, mfc='gray',   mec='gray',   linestyle='None', label='GPMsat')
plt.plot(x_freq_inmet_smn_c_iii,    freq_inmet_smn_c_iii,    marker='.', markersize=4, mfc='black',  mec='black',  linestyle='None', label='station')
plt.plot(x_freq_reg_ictp_i_c_iii,   freq_reg_ictp_i_c_iii,   marker='.', markersize=4, mfc='blue',   mec='blue',   linestyle='None', label='RegCM5-pbl1')
plt.plot(x_freq_reg_ictp_ii_c_iii,  freq_reg_ictp_ii_c_iii,  marker='.', markersize=4, mfc='cyan',   mec='cyan',   linestyle='None', label='RegCM5-pbl2')
plt.plot(x_freq_reg_ictp_iii_c_iii, freq_reg_ictp_iii_c_iii, marker='.', markersize=4, mfc='violet', mec='violet', linestyle='None', label='RegCM5-3km')
plt.plot(x_freq_reg_usp_c_iii,      freq_reg_usp_c_iii,      marker='.', markersize=4, mfc='green',  mec='green',  linestyle='None', label='USP-RegCM471')
plt.plot(x_freq_wrf_ncar_c_iii,     freq_wrf_ncar_c_iii,     marker='.', markersize=4, mfc='red',    mec='red',    linestyle='None', label='NCAR-WRF413')
plt.plot(x_freq_wrf_ucan_c_iii,     freq_wrf_ucan_c_iii,     marker='.', markersize=4, mfc='orange', mec='orange', linestyle='None', label='UCAN-WRF433')
plt.yscale('log')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 2, 4)
plt.plot(x_freq_gpm_c_iv,          freq_gpm_c_iv,          marker='.', markersize=4, mfc='gray',   mec='gray',   linestyle='None', label='GPMsat')
plt.plot(x_freq_inmet_smn_c_iv,    freq_inmet_smn_c_iv,    marker='.', markersize=4, mfc='black',  mec='black',  linestyle='None', label='station')
plt.plot(x_freq_reg_ictp_i_c_iv,   freq_reg_ictp_i_c_iv,   marker='.', markersize=4, mfc='blue',   mec='blue',   linestyle='None', label='RegCM5-pbl1')
plt.plot(x_freq_reg_ictp_ii_c_iv,  freq_reg_ictp_ii_c_iv,  marker='.', markersize=4, mfc='cyan',   mec='cyan',   linestyle='None', label='RegCM5-pbl2')
plt.plot(x_freq_reg_ictp_iii_c_iv, freq_reg_ictp_iii_c_iv, marker='.', markersize=4, mfc='violet', mec='violet', linestyle='None', label='RegCM5-3km')
plt.plot(x_freq_reg_usp_c_iv,      freq_reg_usp_c_iv,      marker='.', markersize=4, mfc='green',  mec='green',  linestyle='None', label='USP-RegCM471')
plt.plot(x_freq_wrf_ncar_c_iv,     freq_wrf_ncar_c_iv,     marker='.', markersize=4, mfc='red',    mec='red',    linestyle='None', label='NCAR-WRF413')
plt.plot(x_freq_wrf_ucan_c_iv,     freq_wrf_ucan_c_iv,     marker='.', markersize=4, mfc='orange', mec='orange', linestyle='None', label='UCAN-WRF433')
plt.title('(d) Cluster IV', loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.xlabel('Intensity of precipitation {0}'.format(legend), fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 2, 5)
plt.plot(x_freq_gpm_c_v,          freq_gpm_c_v,          marker='.', markersize=4, mfc='gray',   mec='gray',   linestyle='None', label='GPMsat')
plt.plot(x_freq_inmet_smn_c_v,    freq_inmet_smn_c_v,    marker='.', markersize=4, mfc='black',  mec='black',  linestyle='None', label='station')
plt.plot(x_freq_reg_ictp_i_c_v,   freq_reg_ictp_i_c_v,   marker='.', markersize=4, mfc='blue',   mec='blue',   linestyle='None', label='RegCM5-pbl1')
plt.plot(x_freq_reg_ictp_ii_c_v,  freq_reg_ictp_ii_c_v,  marker='.', markersize=4, mfc='cyan',   mec='cyan',   linestyle='None', label='RegCM5-pbl2')
plt.plot(x_freq_reg_ictp_iii_c_v, freq_reg_ictp_iii_c_v, marker='.', markersize=4, mfc='violet', mec='violet', linestyle='None', label='RegCM5-3km')
plt.plot(x_freq_reg_usp_c_v,      freq_reg_usp_c_v,      marker='.', markersize=4, mfc='green',  mec='green',  linestyle='None', label='USP-RegCM471')
plt.plot(x_freq_wrf_ncar_c_v,     freq_wrf_ncar_c_v,     marker='.', markersize=4, mfc='red',    mec='red',    linestyle='None', label='NCAR-WRF413')
plt.plot(x_freq_wrf_ucan_c_v,     freq_wrf_ucan_c_v,     marker='.', markersize=4, mfc='orange', mec='orange', linestyle='None', label='UCAN-WRF433')
plt.title('(e) Cluster V', loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.xlabel('Intensity of precipitation {0}'.format(legend), fontsize=font_size, fontweight='bold')

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cp_3km-4km'.format(path)
name_out = 'pyplt_graph_pdf_pr_CSAM-i_CP-RegCM5_1hr_2018-2021.png'.format(var, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
exit()
