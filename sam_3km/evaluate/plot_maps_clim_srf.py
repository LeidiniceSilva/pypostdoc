# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot clim maps"

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
from cartopy import config
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

var = 'pr'
domain = 'SAM-3km'
idt, fdt = '2018', '2021'
dt = '{0}-{1}'.format(idt, fdt)

font_size = 6
path = '/leonardo/home/userexternal/mdasilva/leonardo_work'

		
def import_situ_i(param_i):
	
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
		var_i  = time_i.groupby('time.season').mean(dim='time', skipna=True)
		
		if param_i == 'pre':
			mean_i.append(var_i.values*24)
		else:
			mean_i.append(var_i.values)

	return yy, xx, mean_i
	

def import_situ_ii(param_i):
	
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
		
	return yy, xx, mean_i
	

def import_situ_iii(param_i):
	
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
		
	return yy, xx, mean_i


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

	
# Import model and obs dataset
dict_var = {'pr': ['pre', 'pre', 'precip', 'cmorph', 'tp'],
'tas': ['tmp', 'tmp', 't2m'],
'clt': ['tcc'],
'cll': ['lcc'],
'clm': ['mcc'],
'clh': ['hcc'],
'evspsblpot': ['pev'],
'rsnl': ['msnlwrf'],
'cape': ['cape'],
'cin': ['cin']}

lat, lon, regcm_djf = import_rcm(var, domain, 'RegCM5', 'DJF')
lat, lon, regcm_mam = import_rcm(var, domain, 'RegCM5', 'MAM')
lat, lon, regcm_jja = import_rcm(var, domain, 'RegCM5', 'JJA')
lat, lon, regcm_son = import_rcm(var, domain, 'RegCM5', 'SON')

if var == 'pr':
	lat_i, lon_i, inmet_i = import_situ_i(dict_var[var][0])
	lat_ii, lon_ii, smn_ii = import_situ_ii(dict_var[var][0])
	lat_iii, lon_iii, smn_iii = import_situ_iii(dict_var[var][0])

	lat_yy = lat_i + lat_ii + lat_iii
	lon_xx = lon_i + lon_ii + lon_iii
	weather_stations = inmet_i + smn_ii + smn_iii

	ws_djf, ws_mam, ws_jja, ws_son = [], [], [], []
	for i in range(0, len(lat_yy)):
		ws_djf.append(weather_stations[i][0])
		ws_mam.append(weather_stations[i][2])
		ws_jja.append(weather_stations[i][1])
		ws_son.append(weather_stations[i][3])

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

elif var == 'tas':
	lat_i, lon_i, inmet_i = import_situ_i(dict_var[var][0])

	lat_yy = lat_i
	lon_xx = lon_i
	weather_stations = inmet_i

	ws_djf, ws_mam, ws_jja, ws_son = [], [], [], []
	for i in range(0, len(lat_yy)):
		ws_djf.append(weather_stations[i][0])
		ws_mam.append(weather_stations[i][2])
		ws_jja.append(weather_stations[i][1])
		ws_son.append(weather_stations[i][3])

	lat, lon, cru_djf = import_obs(dict_var[var][1], domain, 'CRU', 'DJF')
	lat, lon, cru_mam = import_obs(dict_var[var][1], domain, 'CRU', 'MAM')
	lat, lon, cru_jja = import_obs(dict_var[var][1], domain, 'CRU', 'JJA')
	lat, lon, cru_son = import_obs(dict_var[var][1], domain, 'CRU', 'SON')

	lat, lon, era5_djf = import_obs(dict_var[var][2], domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_obs(dict_var[var][2], domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_obs(dict_var[var][2], domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_obs(dict_var[var][2], domain, 'ERA5', 'SON')
else:
	lat, lon, era5_djf = import_obs(dict_var[var][0], domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_obs(dict_var[var][0], domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_obs(dict_var[var][0], domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_obs(dict_var[var][0], domain, 'ERA5', 'SON')

	lat, lon, regcm_djf = import_rcm(var, domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_rcm(var, domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_rcm(var, domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_rcm(var, domain, 'RegCM5', 'SON')

# Plot figure 
def configure_subplot(ax):

	lon_min = np.round(np.min(lon), 1)
	lon_max = np.round(np.max(lon), 1)
	lat_min = np.round(np.min(lat), 1)
	lat_max = np.round(np.max(lat), 1)

	ax.set_extent([np.min(lon), np.max(lon), np.min(lat), np.max(lat)], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(lon_min,lon_max,10), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(lat_min,lat_max,5), crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.tick_params(labelsize=font_size)
	ax.add_feature(cfeat.BORDERS, linewidth=0.5)
	ax.coastlines(linewidth=0.5)	
	ax.grid(c='k', ls='--', alpha=0.4)

color=['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']

dict_plot = {'pr': ['Precipitation (mm d$^-$$^1$)', np.arange(0, 18, 1), matplotlib.colors.ListedColormap(color)],
'tas': ['Air temperature (°C)', np.arange(-6, 34, 2), cm.jet],
'clt': ['Total cloud cover (0-1)', np.arange(0, 1.05, 0.05), cm.Greys],
'cll': ['Low cloud cover (0-1)', np.arange(0, 1.05, 0.05), cm.Greys],
'clm': ['Medium cloud cover (0-1)', np.arange(0, 1.05, 0.05), cm.Greys],
'clh': ['High cloud cover (0-1)', np.arange(0, 1.05, 0.05), cm.Greys],
'evspsblpot': ['Potential evapotranspiration (mm d$^-$$^1$)', np.arange(0, 10.5, 0.5), cm.jet],
'rsnl': ['Surface net upward longwave flux (W mm$^-$$^2$)', np.arange(0, 220, 10), cm.rainbow],
'CAPE': ['Convective Available Potential Energy (J kg$^-$$^1$)', np.arange(0, 2550, 50), cm.jet],
'CIN': ['Convective inhibition (J kg$^-$$^1$)', np.arange(0, 675, 25), cm.jet]}

if var == 'pr':
	fig, axes = plt.subplots(4, 5, figsize=(15, 8), subplot_kw={'projection': ccrs.PlateCarree()})

	ax1 = axes[0, 0]
	plt_map = ax1.contourf(lon, lat, regcm_djf, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt_sct = ax1.scatter(lon_xx, lat_yy, 4, ws_djf, cmap=dict_plot[var][2], edgecolors='black', linewidth=0.5, marker='o', vmin=0, vmax=18) 
	ax1.set_title(u'(a) RegCM5(shaded) INMET(dots) DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax1)

	ax2 = axes[0, 1]
	plt_map = ax2.contourf(lon, lat, cru_djf, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax2.set_title(u'(b) CRU DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax2)

	ax3 = axes[0, 2]
	plt_map = ax3.contourf(lon, lat, cpc_djf, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax3.set_title(u'(c) CPC DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax3)

	ax4 = axes[0, 3] 
	plt_map = ax4.contourf(lon, lat, cmorph_djf, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax4.set_title(u'(d) CMORPH DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax4)

	ax5 = axes[0, 4]
	plt_map = ax5.contourf(lon, lat, era5_djf, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax5.set_title(u'(e) ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax5)

	ax6 = axes[1, 0]
	plt_map = ax6.contourf(lon, lat, regcm_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt_sct = ax6.scatter(lon_xx, lat_yy, 4, ws_mam, cmap=dict_plot[var][2], edgecolors='black', linewidth=0.5, marker='o', vmin=0, vmax=18) 
	ax6.set_title(u'(f) RegCM5(shaded) INMET(dots) MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax6)

	ax7 = axes[1, 1]
	plt_map = ax7.contourf(lon, lat, cru_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax7.set_title(u'(g) CRU MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax7)

	ax8 = axes[1, 2]
	plt_map = ax8.contourf(lon, lat, cpc_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax8.set_title(u'(h) CPC MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax8)

	ax9 = axes[1, 3] 
	plt_map = ax9.contourf(lon, lat, cmorph_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax9.set_title(u'(i) CMORPH MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax9)

	ax10 = axes[1, 4]
	plt_map = ax10.contourf(lon, lat, era5_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax10.set_title(u'(j) ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax10)

	ax11 = axes[2, 0]
	plt_map = ax11.contourf(lon, lat, regcm_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt_sct = ax11.scatter(lon_xx, lat_yy, 4, ws_jja, cmap=dict_plot[var][2], edgecolors='black', linewidth=0.5, marker='o', vmin=0, vmax=18) 
	ax11.set_title(u'(k) RegCM5(shaded) INMET(dots) JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax11)

	ax12 = axes[2, 1]
	plt_map = ax12.contourf(lon, lat, cru_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax12.set_title(u'(l) CRU JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax12)

	ax13 = axes[2, 2]
	plt_map = ax13.contourf(lon, lat, cpc_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax13.set_title(u'(m) CPC JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax13)

	ax14 = axes[2, 3] 
	plt_map = ax14.contourf(lon, lat, cmorph_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax14.set_title(u'(n) CMORPH JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax14)

	ax15 = axes[2, 4]
	plt_map = ax15.contourf(lon, lat, era5_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax15.set_title(u'(o) ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax15)

	ax16 = axes[3, 0]
	plt_map = ax16.contourf(lon, lat, regcm_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt_sct = ax16.scatter(lon_xx, lat_yy, 4, ws_son, cmap=dict_plot[var][2], edgecolors='black', linewidth=0.5, marker='o', vmin=0, vmax=18)
	ax16.set_title(u'(p) RegCM5(shaded) INMET(dots) SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax16)

	ax17 = axes[3, 1]
	plt_map = ax17.contourf(lon, lat, cru_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax17.set_title(u'(q) CRU SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax17)

	ax18 = axes[3, 2]
	plt_map = ax18.contourf(lon, lat, cpc_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax18.set_title(u'(r) CPC SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax18)

	ax19 = axes[3, 3] 
	plt_map = ax19.contourf(lon, lat, cmorph_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax19.set_title(u'(s) CMORPH SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax19)

	ax20 = axes[3, 4]
	plt_map = ax20.contourf(lon, lat, era5_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax20.set_title(u'(t) ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax20)
	
elif var == 'tas':
	fig, axes = plt.subplots(4, 3, figsize=(8, 8), subplot_kw={'projection': ccrs.PlateCarree()})

	ax1 = axes[0, 0]
	plt_map = ax1.contourf(lon, lat, regcm_djf, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt_sct = ax1.scatter(lon_xx, lat_yy, 4, ws_djf, cmap=dict_plot[var][2], edgecolors='black', linewidth=0.5, marker='o', vmin=10, vmax=32)
	ax1.set_title(u'(a) RegCM5(shaded) INMET(dots) DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax1)

	ax2 = axes[0, 1]
	plt_map = ax2.contourf(lon, lat, cru_djf, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax2.set_title(u'(b) CRU DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax2)

	ax3 = axes[0, 2]
	plt_map = ax3.contourf(lon, lat, era5_djf, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax3.set_title(u'(c) ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax3)

	ax4 = axes[1, 0]
	plt_map = ax4.contourf(lon, lat, regcm_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt_sct = ax4.scatter(lon_xx, lat_yy, 4, ws_mam, cmap=dict_plot[var][2], edgecolors='black', linewidth=0.5, marker='o', vmin=10, vmax=32)
	ax4.set_title(u'(d) RegCM5(shaded) INMET(dots) MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax4)

	ax5 = axes[1, 1]
	plt_map = ax5.contourf(lon, lat, cru_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax5.set_title(u'(e) CRU MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax5)

	ax6 = axes[1, 2]
	plt_map = ax6.contourf(lon, lat, era5_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax6.set_title(u'(f) ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax6)

	ax7 = axes[2, 0]
	plt_map = ax7.contourf(lon, lat, regcm_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt_sct = ax7.scatter(lon_xx, lat_yy, 4, ws_jja, cmap=dict_plot[var][2], edgecolors='black', linewidth=0.5, marker='o', vmin=10, vmax=32)
	ax7.set_title(u'(g) RegCM5(shaded) INMET(dots) JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax7)

	ax8 = axes[2, 1]
	plt_map = ax8.contourf(lon, lat, cru_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax8.set_title(u'(h) CRU JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax8)

	ax9 = axes[2, 2]
	plt_map = ax9.contourf(lon, lat, era5_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax9.set_title(u'(i) ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax9)

	ax10 = axes[3, 0]
	plt_map = ax10.contourf(lon, lat, regcm_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt_sct = ax10.scatter(lon_xx, lat_yy, 4, ws_son, cmap=dict_plot[var][2], edgecolors='black', linewidth=0.5, marker='o', vmin=10, vmax=32)
	ax10.set_title(u'(j) RegCM5(shaded) INMET(dots) SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax10)

	ax11 = axes[3, 1]
	plt_map = ax11.contourf(lon, lat, cru_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax11.set_title(u'(k) CRU SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax11)

	ax12 = axes[3, 2]
	plt_map = ax12.contourf(lon, lat, era5_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax12.set_title(u'(l) ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax12)

else:
	fig, axes = plt.subplots(4, 2, figsize=(10, 12), subplot_kw={'projection': ccrs.PlateCarree()})

	ax1 = axes[0, 0]
	plt_map = ax1.contourf(lon, lat, regcm_djf, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax1.set_title(u'(a) RegCM5 DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax1)

	ax2 = axes[0, 1]
	plt_map = ax2.contourf(lon, lat, era5_djf, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax2.set_title(u'(b) ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax2)

	ax3 = axes[1, 0]
	plt_map = ax3.contourf(lon, lat, regcm_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax3.set_title(u'(c) RegCM5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax3)

	ax4 = axes[1, 1]
	plt_map = ax4.contourf(lon, lat, era5_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax4.set_title(u'(d) ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax4)

	ax5 = axes[2, 0]
	plt_map = ax5.contourf(lon, lat, regcm_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax5.set_title(u'(e) RegCM5 DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax5)

	ax6 = axes[2, 1]
	plt_map = ax6.contourf(lon, lat, era5_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax6.set_title(u'(f) ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax6)

	ax7 = axes[3, 0]
	plt_map = ax7.contourf(lon, lat, regcm_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax7.set_title(u'(g) RegCM5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax7)

	ax8 = axes[3, 1]
	plt_map = ax8.contourf(lon, lat, era5_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax8.set_title(u'(h) ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax8)

# Set colobar
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.3, 0.015, 0.4]))
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_maps_clim_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()



