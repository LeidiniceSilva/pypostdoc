# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot bias maps"

import os
import netCDF4
import numpy as np
import xarray as xr
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat

from netCDF4 import Dataset
from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii
from import_climate_tools import compute_mbe
from cartopy import config
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

var = 'sfcWindmax'
domain = 'SAM-3km'
idt, fdt = '2018', '2021'
dt = '{0}-{1}'.format(idt, fdt)

font_size = 8
path = '/leonardo/home/userexternal/mdasilva/leonardo_work'
	

def import_situ_i(param_i, param_ii, domain, dataset):
	
	yy, xx = [], []
	mean_i, mean_ii = [], []
	
	skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,7793,100,105,106,107,112,117,124,135,137,139,
	149,152,155,158,168,174,177,183,186,199,204,210,212,224,226,239,240,248,249,253,254,276,277,280,293,298,
	303,305,306,308,319,334,335,341,343,359,362,364,384,393,396,398,399,400,402,413,416,417,422,423,426,427,
	443,444,446,451,453,457,458,467,474,479,483,488,489,490,495,505,509,513,514,516,529,534,544,559,566]
	
	for station in range(1, 567):
		print(station, inmet[station][1])
		
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue
			
		yy.append(inmet[station][2])
		xx.append(inmet[station][3])

		arq_i  = xr.open_dataset('{0}/FPS_SESA/database/obs/inmet/inmet_nc/hourly/{1}/'.format(path, param_i) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(param_i, inmet[station][0]))
		data_i = arq_i[param_i]
		time_i = data_i.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_i  = time_i.groupby('time.season').mean(dim='time')
		
		if param_i == 'pre':
			mean_i.append(var_i.values*24)
		else:
			mean_i.append(var_i.values)

		arq_ii  = xr.open_dataset('{0}/SAM-3km/postproc/evaluate/rcm/'.format(path) + '{0}_{1}_{2}_mon_{3}_lonlat.nc'.format(param_ii, domain, dataset, dt))
		data_ii = arq_ii[param_ii]
		data_ii = data_ii.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		time_ii = data_ii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_ii  = time_ii.groupby('time.season').mean(dim='time')
		mean_ii.append(var_ii.values)
		
	return yy, xx, mean_i, mean_ii
	

def import_situ_ii(param_i, param_ii, domain, dataset):
	
	yy, xx = [], []
	mean_i, mean_ii = [], []
	
	for station in range(1, 73):
		print(station, smn_i[station][0])

		yy.append(smn_i[station][1])
		xx.append(smn_i[station][2])

		arq_i  = xr.open_dataset('{0}/FPS_SESA/database/obs/smn_i/smn_nc/'.format(path) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(param_i, smn_i[station][0]))
		data_i = arq_i[param_i]
		time_i = data_i.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_i  = time_i.groupby('time.season').mean(dim='time')
		mean_i.append(var_i.values*24)

		arq_ii  = xr.open_dataset('{0}/SAM-3km/postproc/evaluate/rcm/'.format(path) + '{0}_{1}_{2}_mon_{3}_lonlat.nc'.format(param_ii, domain, dataset, dt))
		data_ii = arq_ii[param_ii]
		data_ii = data_ii.sel(lat=slice(smn_i[station][1]-0.03,smn_i[station][1]+0.03),lon=slice(smn_i[station][2]-0.03,smn_i[station][2]+0.03)).mean(('lat','lon'))
		time_ii = data_ii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_ii  = time_ii.groupby('time.season').mean(dim='time')
		mean_ii.append(var_ii.values)
		
	return yy, xx, mean_i, mean_ii
	

def import_situ_iii(param_i, param_ii, domain, dataset):
	
	yy, xx = [], []
	mean_i, mean_ii = [], []

	skip_list = (86,87,88,89,90,91,95,96,97,98,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128)	
		
	for station in range(1, 129):
		print(station, smn_ii[station][0])

		if station in skip_list:
			continue
													
		yy.append(smn_ii[station][1])
		xx.append(smn_ii[station][2])

		arq_i  = xr.open_dataset('{0}/FPS_SESA/database/obs/smn_ii/smn_nc/{1}/'.format(path, param_i) + '{0}_{1}_D_1979-01-01_2021-12-31.nc'.format(param_i, smn_ii[station][0]))
		data_i = arq_i[param_i]
		time_i = data_i.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_i  = time_i.groupby('time.season').mean(dim='time')
		mean_i.append(var_i.values)

		arq_ii  = xr.open_dataset('{0}/SAM-3km/postproc/evaluate/rcm/'.format(path) + '{0}_{1}_{2}_mon_{3}_lonlat.nc'.format(param_ii, domain, dataset, dt))
		data_ii = arq_ii[param_ii]
		data_ii = data_ii.sel(lat=slice(smn_ii[station][1]-0.03,smn_ii[station][1]+0.03),lon=slice(smn_ii[station][2]-0.03,smn_ii[station][2]+0.03)).mean(('lat','lon'))
		time_ii = data_ii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_ii  = time_ii.groupby('time.season').mean(dim='time')
		mean_ii.append(var_ii.values)
		
	return yy, xx, mean_i, mean_ii


def import_obs(param, domain, dataset, season):

	arq   = '{0}/SAM-3km/postproc/evaluate/obs/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]

	if param == 'msnlwrf':
		mean = var[:][0,:,:]*(-1)
	elif param == 'pev':
		mean_i = var[:][0,:,:]*(-1000)
		mean_i = np.abs(mean_i)
		mask_nc = Dataset('{0}/SAM-3km/postproc/evaluate/rcm/sea_land_mask_lonlat.nc'.format(path))
		lsm = mask_nc.variables['lsm'][:]
		ocean_mask = (lsm == 0)
		mean_ = np.where(ocean_mask[None, :, :] == 1, np.nan, mean_i)
		mean = mean_[0,0,:,:]	
	else:
		mean = var[:][0,:,:]
	
	return lat, lon, mean


def import_rcm(param, domain, dataset, season):

	arq   = '{0}/SAM-3km/postproc/evaluate/rcm/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]

	if param == 'tas':
		mean = var[:][0,0,:,:]
	elif param == 'evspsblpot':
		mask_nc = Dataset('{0}/SAM-3km/postproc/evaluate/rcm/sea_land_mask_lonlat.nc'.format(path))
		lsm = mask_nc.variables['lsm'][:]
		ocean_mask = (lsm == 0)
		mean_ = np.where(ocean_mask[None, :, :] == 1, np.nan, var[:][0,:,:])
		mean = mean_[0,0,:,:]	
	else:
		mean = var[:][0,:,:]

	return lat, lon, mean


def configure_subplot(ax):

	ax.set_extent([-80, -34, -38, -8], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(-80,-34,12), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(-38,-8,6), crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.tick_params(labelsize=font_size)
	ax.add_feature(cfeat.BORDERS)
	ax.coastlines()	
	ax.grid(c='k', ls='--', alpha=0.4)


# Import model and obs dataset
dict_var = {'pr': ['pre', 'pre', 'precip', 'cmorph', 'tp'],
'tas': ['tmp', 'tmp', 't2m'],
'clt': ['tcc'],
'cll': ['lcc'],
'clm': ['mcc'],
'clh': ['hcc'],
'sfcWindmax': ['u10max', 'v10max'],
'evspsblpot': ['pev'],
'rsnl': ['msnlwrf']}

if var == 'pr':
	lat_i, lon_i, inmet_v1, regcm_v1 = import_situ_i(dict_var[var][0], var, domain, 'RegCM5')
	lat_ii, lon_ii, smn_v2, regcm_v2 = import_situ_ii(dict_var[var][0], var, domain, 'RegCM5')
	lat_iii, lon_iii, smn_v3, regcm_v3 = import_situ_iii(dict_var[var][0], var, domain, 'RegCM5')

	lat_yy = lat_i + lat_ii + lat_iii
	lon_xx = lon_i + lon_ii + lon_iii
	inmet_smn = inmet_v1 + smn_v2 + smn_v3
	regcm_latlon = regcm_v1 + regcm_v2 + regcm_v3

	mbe_djf_regcm_inmet_smn, mbe_mam_regcm_inmet_smn, mbe_jja_regcm_inmet_smn, mbe_son_regcm_inmet_smn = [], [], [], []
	for i in range(0, len(regcm_latlon)):
		mbe_djf_regcm_inmet_smn.append(compute_mbe(regcm_latlon[i][0], inmet_smn[i][0]))
		mbe_mam_regcm_inmet_smn.append(compute_mbe(regcm_latlon[i][2], inmet_smn[i][2]))
		mbe_jja_regcm_inmet_smn.append(compute_mbe(regcm_latlon[i][1], inmet_smn[i][1]))
		mbe_son_regcm_inmet_smn.append(compute_mbe(regcm_latlon[i][3], inmet_smn[i][3]))

	lat, lon, cru_djf = import_obs(dict_var[var][1], domain, 'CRU', 'DJF')
	lat, lon, cru_mam = import_obs(dict_var[var][1], domain, 'CRU', 'MAM')
	lat, lon, cru_jja = import_obs(dict_var[var][1], domain, 'CRU', 'JJA')
	lat, lon, cru_son = import_obs(dict_var[var][1], domain, 'CRU', 'SON')

	lat, lon, cpc_djf = import_obs(dict_var[var][2], domain, 'CPC', 'DJF')
	lat, lon, cpc_mam = import_obs(dict_var[var][2], domain, 'CPC', 'MAM')
	lat, lon, cpc_jja = import_obs(dict_var[var][2], domain, 'CPC', 'JJA')
	lat, lon, cpc_son = import_obs(dict_var[var][2], domain, 'CPC', 'SON')

	lat, lon, cmorph_djf = import_obs(dict_var[var][3], domain, 'CMORPH', 'DJF')
	lat, lon, cmorph_mam = import_obs(dict_var[var][3], domain, 'CMORPH', 'MAM')
	lat, lon, cmorph_jja = import_obs(dict_var[var][3], domain, 'CMORPH', 'JJA')
	lat, lon, cmorph_son = import_obs(dict_var[var][3], domain, 'CMORPH', 'SON')

	lat, lon, era5_djf = import_obs(dict_var[var][4], domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_obs(dict_var[var][4], domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_obs(dict_var[var][4], domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_obs(dict_var[var][4], domain, 'ERA5', 'SON')

	lat, lon, regcm_djf = import_rcm(var, domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_rcm(var, domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_rcm(var, domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_rcm(var, domain, 'RegCM5', 'SON')
	
	mbe_djf_regcm_cru = compute_mbe(regcm_djf, cru_djf)
	mbe_mam_regcm_cru = compute_mbe(regcm_mam, cru_mam)
	mbe_jja_regcm_cru = compute_mbe(regcm_jja, cru_jja)
	mbe_son_regcm_cru = compute_mbe(regcm_son, cru_son)	

	mbe_djf_regcm_cpc = compute_mbe(regcm_djf, cpc_djf)
	mbe_mam_regcm_cpc = compute_mbe(regcm_mam, cpc_mam)
	mbe_jja_regcm_cpc = compute_mbe(regcm_jja, cpc_jja)
	mbe_son_regcm_cpc = compute_mbe(regcm_son, cpc_son)	
	
	mbe_djf_regcm_cmorph = compute_mbe(regcm_djf, cmorph_djf)
	mbe_mam_regcm_cmorph = compute_mbe(regcm_mam, cmorph_mam)
	mbe_jja_regcm_cmorph = compute_mbe(regcm_jja, cmorph_jja)
	mbe_son_regcm_cmorph = compute_mbe(regcm_son, cmorph_son)
		
	mbe_djf_regcm_era5 = compute_mbe(regcm_djf, era5_djf)
	mbe_mam_regcm_era5 = compute_mbe(regcm_mam, era5_mam)
	mbe_jja_regcm_era5 = compute_mbe(regcm_jja, era5_jja)
	mbe_son_regcm_era5 = compute_mbe(regcm_son, era5_son)		
elif var == 'tas':
	lat_i, lon_i, inmet_i, regcm_i = import_situ_i(dict_var[var][0], var, domain, 'RegCM5')
	mbe_djf_regcm_inmet, mbe_mam_regcm_inmet, mbe_jja_regcm_inmet, mbe_son_regcm_inmet = [], [], [], []

	for i in range(0, 298):
		mbe_djf_regcm_inmet.append(compute_mbe(regcm_i[i][0], inmet_i[i][0]))
		mbe_mam_regcm_inmet.append(compute_mbe(regcm_i[i][2], inmet_i[i][2]))
		mbe_jja_regcm_inmet.append(compute_mbe(regcm_i[i][1], inmet_i[i][1]))
		mbe_son_regcm_inmet.append(compute_mbe(regcm_i[i][3], inmet_i[i][3]))

	lat, lon, cru_djf = import_obs(dict_var[var][1], domain, 'CRU', 'DJF')
	lat, lon, cru_mam = import_obs(dict_var[var][1], domain, 'CRU', 'MAM')
	lat, lon, cru_jja = import_obs(dict_var[var][1], domain, 'CRU', 'JJA')
	lat, lon, cru_son = import_obs(dict_var[var][1], domain, 'CRU', 'SON')

	lat, lon, era5_djf = import_obs(dict_var[var][2], domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_obs(dict_var[var][2], domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_obs(dict_var[var][2], domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_obs(dict_var[var][2], domain, 'ERA5', 'SON')

	lat, lon, regcm_djf = import_rcm(var, domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_rcm(var, domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_rcm(var, domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_rcm(var, domain, 'RegCM5', 'SON')
	
	mbe_djf_regcm_cru = compute_mbe(regcm_djf, cru_djf)
	mbe_mam_regcm_cru = compute_mbe(regcm_mam, cru_mam)
	mbe_jja_regcm_cru = compute_mbe(regcm_jja, cru_jja)
	mbe_son_regcm_cru = compute_mbe(regcm_son, cru_son)	
		
	mbe_djf_regcm_era5 = compute_mbe(regcm_djf, era5_djf)
	mbe_mam_regcm_era5 = compute_mbe(regcm_mam, era5_mam)
	mbe_jja_regcm_era5 = compute_mbe(regcm_jja, era5_jja)
	mbe_son_regcm_era5 = compute_mbe(regcm_son, era5_son)
elif var == 'sfcWindmax':
	lat, lon, u_era5_djf = import_obs(dict_var[var][0], domain, 'ERA5', 'DJF')
	lat, lon, u_era5_mam = import_obs(dict_var[var][0], domain, 'ERA5', 'MAM')
	lat, lon, u_era5_jja = import_obs(dict_var[var][0], domain, 'ERA5', 'JJA')
	lat, lon, u_era5_son = import_obs(dict_var[var][0], domain, 'ERA5', 'SON')

	lat, lon, v_era5_djf = import_obs(dict_var[var][1], domain, 'ERA5', 'DJF')
	lat, lon, v_era5_mam = import_obs(dict_var[var][1], domain, 'ERA5', 'MAM')
	lat, lon, v_era5_jja = import_obs(dict_var[var][1], domain, 'ERA5', 'JJA')
	lat, lon, v_era5_son = import_obs(dict_var[var][1], domain, 'ERA5', 'SON')

	uv_era5_djf = np.sqrt(u_era5_djf**2 + v_era5_djf**2)
	uv_era5_mam = np.sqrt(u_era5_mam**2 + v_era5_mam**2)
	uv_era5_jja = np.sqrt(u_era5_jja**2 + v_era5_jja**2)
	uv_era5_son = np.sqrt(u_era5_son**2 + v_era5_son**2)

	lat, lon, uv_regcm_djf = import_rcm(var, domain, 'RegCM5', 'DJF')
	lat, lon, uv_regcm_mam = import_rcm(var, domain, 'RegCM5', 'MAM')
	lat, lon, uv_regcm_jja = import_rcm(var, domain, 'RegCM5', 'JJA')
	lat, lon, uv_regcm_son = import_rcm(var, domain, 'RegCM5', 'SON')

	mbe_djf_regcm_era5 = compute_mbe(uv_regcm_djf, uv_era5_djf)
	mbe_mam_regcm_era5 = compute_mbe(uv_regcm_mam, uv_era5_mam)
	mbe_jja_regcm_era5 = compute_mbe(uv_regcm_jja, uv_era5_jja)
	mbe_son_regcm_era5 = compute_mbe(uv_regcm_son, uv_era5_son)
else:
	lat, lon, era5_djf = import_obs(dict_var[var][0], domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_obs(dict_var[var][0], domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_obs(dict_var[var][0], domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_obs(dict_var[var][0], domain, 'ERA5', 'SON')

	lat, lon, regcm_djf = import_rcm(var, domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_rcm(var, domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_rcm(var, domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_rcm(var, domain, 'RegCM5', 'SON')

	mbe_djf_regcm_era5 = compute_mbe(regcm_djf, era5_djf)
	mbe_mam_regcm_era5 = compute_mbe(regcm_mam, era5_mam)
	mbe_jja_regcm_era5 = compute_mbe(regcm_jja, era5_jja)
	mbe_son_regcm_era5 = compute_mbe(regcm_son, era5_son)

# Plot figure   
dict_plot = {'pr': ['Bias of  precipitation (mm d$^-$$^1$)', np.arange(-10, 11, 1), cm.BrBG],
'tas': ['Bias of air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'clt': ['Bias of total cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'cll': ['Bias of low cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'clm': ['Bias of medium cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'clh': ['Bias of high cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'sfcWindmax': ['Bias of maximum wind speed at 10 meters (m s$^-$$^1$)', np.arange(-5, 5.5, 0.5), cm.bwr],
'evspsblpot': ['Bias of potential evapotranspiration (mm d$^-$$^1$)', np.arange(-5, 5.5, 0.5), cm.bwr],
'rsnl': ['Bias of surface net upward longwave flux (W mm$^-$$^2$)', np.arange(-80, 90, 10), cm.RdBu_r]}

if var == 'pr':
	fig, axes = plt.subplots(4, 5, figsize=(14, 8), subplot_kw={'projection': ccrs.PlateCarree()})

	ax1 = axes[0, 0]
	plt_map = ax1.scatter(lon_xx, lat_yy, 4, mbe_djf_regcm_inmet_smn, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=10) 
	ax1.set_title(u'(a) RegCM5-INMET DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax1)

	ax2 = axes[0, 1]
	plt_map = ax2.contourf(lon, lat, mbe_djf_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax2.set_title(u'(b) RegCM5-CRU DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax2)

	ax3 = axes[0, 2]
	plt_map = ax3.contourf(lon, lat, mbe_djf_regcm_cpc, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax3.set_title(u'(c) RegCM5-CPC DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax3)

	ax4 = axes[0, 3] 
	plt_map = ax4.contourf(lon, lat, mbe_djf_regcm_cmorph, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax4.set_title(u'(d) RegCM5-CMORPH DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax4)

	ax5 = axes[0, 4]
	plt_map = ax5.contourf(lon, lat, mbe_djf_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax5.set_title(u'(e) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax5)

	ax6 = axes[1, 0]
	plt_map = ax6.scatter(lon_xx, lat_yy, 4, mbe_mam_regcm_inmet_smn, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=10) 
	ax6.set_title(u'(f) RegCM5-INMET MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax6)

	ax7 = axes[1, 1]
	plt_map = ax7.contourf(lon, lat, mbe_mam_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax7.set_title(u'(g) RegCM5-CRU MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax7)

	ax8 = axes[1, 2]
	plt_map = ax8.contourf(lon, lat, mbe_mam_regcm_cpc, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax8.set_title(u'(h) RegCM5-CPC MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax8)

	ax9 = axes[1, 3] 
	plt_map = ax9.contourf(lon, lat, mbe_mam_regcm_cmorph, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax9.set_title(u'(i) RegCM5-CMORPH MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax9)

	ax10 = axes[1, 4]
	plt_map = ax10.contourf(lon, lat, mbe_mam_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax10.set_title(u'(j) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax10)

	ax11 = axes[2, 0]
	plt_map = ax11.scatter(lon_xx, lat_yy, 4, mbe_jja_regcm_inmet_smn, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=10) 
	ax11.set_title(u'(k) RegCM5-INMET JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax11)

	ax12 = axes[2, 1]
	plt_map = ax12.contourf(lon, lat, mbe_jja_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax12.set_title(u'(l) RegCM5-CRU JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax12)

	ax13 = axes[2, 2]
	plt_map = ax13.contourf(lon, lat, mbe_jja_regcm_cpc, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax13.set_title(u'(m) RegCM5-CPC JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax13)

	ax14 = axes[2, 3] 
	plt_map = ax14.contourf(lon, lat, mbe_jja_regcm_cmorph, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax14.set_title(u'(n) RegCM5-CMORPH JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax14)

	ax15 = axes[2, 4]
	plt_map = ax15.contourf(lon, lat, mbe_jja_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax15.set_title(u'(o) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax15)

	ax16 = axes[3, 0]
	plt_map = ax16.scatter(lon_xx, lat_yy, 4, mbe_son_regcm_inmet_smn, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=10) 
	ax16.set_title(u'(p) RegCM5-INMET SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax16)

	ax17 = axes[3, 1]
	plt_map = ax17.contourf(lon, lat, mbe_son_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax17.set_title(u'(q) RegCM5-CRU SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax17)

	ax18 = axes[3, 2]
	plt_map = ax18.contourf(lon, lat, mbe_son_regcm_cpc, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax18.set_title(u'(r) RegCM5-CPC SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax18)

	ax19 = axes[3, 3] 
	plt_map = ax19.contourf(lon, lat, mbe_son_regcm_cmorph, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax19.set_title(u'(s) RegCM5-CMORPH SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax19)

	ax20 = axes[3, 4]
	plt_map = ax20.contourf(lon, lat, mbe_son_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax20.set_title(u'(t) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax20)
	
elif var == 'tas':
	fig, axes = plt.subplots(4, 3, figsize=(8, 8), subplot_kw={'projection': ccrs.PlateCarree()})

	ax1 = axes[0, 0]
	plt_map = ax1.scatter(lon_i, lat_i, 4, mbe_djf_regcm_inmet, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=10) 
	ax1.set_title(u'(a) RegCM5-INMET DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax1)

	ax2 = axes[0, 1]
	plt_map = ax2.contourf(lon, lat, mbe_djf_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax2.set_title(u'(b) RegCM5-CRU DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax2)

	ax3 = axes[0, 2]
	plt_map = ax3.contourf(lon, lat, mbe_djf_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax3.set_title(u'(c) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax3)

	ax4 = axes[1, 0]
	plt_map = ax4.scatter(lon_i, lat_i, 4, mbe_mam_regcm_inmet, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=10) 
	ax4.set_title(u'(d) RegCM5-INMET MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax4)

	ax5 = axes[1, 1]
	plt_map = ax5.contourf(lon, lat, mbe_mam_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax5.set_title(u'(e) RegCM5-CRU MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax5)

	ax6 = axes[1, 2]
	plt_map = ax6.contourf(lon, lat, mbe_mam_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax6.set_title(u'(f) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax6)

	ax7 = axes[2, 0]
	plt_map = ax7.scatter(lon_i, lat_i, 4, mbe_jja_regcm_inmet, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=10) 
	ax7.set_title(u'(g) RegCM5-INMET JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax7)

	ax8 = axes[2, 1]
	plt_map = ax8.contourf(lon, lat, mbe_jja_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax8.set_title(u'(h) RegCM5-CRU JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax8)

	ax9 = axes[2, 2]
	plt_map = ax9.contourf(lon, lat, mbe_jja_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax9.set_title(u'(i) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax9)

	ax10 = axes[3, 0]
	plt_map = ax10.scatter(lon_i, lat_i, 4, mbe_son_regcm_inmet, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=10) 
	ax10.set_title(u'(j) RegCM5-INMET SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax10)

	ax11 = axes[3, 1]
	plt_map = ax11.contourf(lon, lat, mbe_son_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax11.set_title(u'(k) RegCM5-CRU SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax11)

	ax12 = axes[3, 2]
	plt_map = ax12.contourf(lon, lat, mbe_son_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax12.set_title(u'(l) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax12)

else:
	fig, axes = plt.subplots(4, 1, figsize=(5, 12), subplot_kw={'projection': ccrs.PlateCarree()})
	axes = axes.flatten()

	ax1 = axes[0]
	plt_map = ax1.contourf(lon, lat, mbe_djf_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax1.set_title(u'(a) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax1)

	ax2 = axes[1]
	plt_map = ax2.contourf(lon, lat, mbe_mam_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax2.set_title(u'(b) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax2)

	ax3 = axes[2] 
	plt_map = ax3.contourf(lon, lat, mbe_jja_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax3.set_title(u'(c) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax3)

	ax4 = axes[3]
	plt_map = ax4.contourf(lon, lat, mbe_son_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax4.set_title(u'(d) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax4)

# Set colobar
if var == 'clt' or var == 'cll' or var == 'clm' or var == 'clh' or var == 'evspsblpot' or var == 'rsnl':
	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.84, 0.3, 0.03, 0.4]))
	cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)
else:
	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.3, 0.015, 0.4]))
	cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
