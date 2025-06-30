# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot bias maps"

import os
import netCDF4
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat

from cartopy import config
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from import_climate_tools import compute_mbe

var = 'rlds'
domain = 'CSAM-3'
idt, fdt = '2000', '2009'
dt = '{0}-{1}'.format(idt, fdt)
font_size = 6

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5'

def import_obs(param, domain, dataset, season):

	arq   = '{0}/postproc/evaluate/obs/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]

	if param == 'pev':
		mean = var[:][0,:,:]*(-1000)
	elif param == 'msnlwrf':
		mean = var[:][0,:,:]*(-1)
	else:
		mean = var[:][0,:,:]
	
	return lat, lon, mean


def import_rcm(param, domain, dataset, season):

	arq   = '{0}/postproc/evaluate/rcm/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

	return lat, lon, mean


# Import model and obs dataset
dict_var = {'pr': ['pre', 'precip', 'cmorph', 'precipitation', 'tp'],
'tas': ['tmp', 't2m'],
'tasmax': ['tmx', 'tmax', 'tasmax'],
'tasmin': ['tmn', 'tmin', 'tasmin'],
'clt': ['cld', 'tcc'],
'cll': ['lcc'],
'clm': ['mcc'],
'clh': ['hcc'],
'evspsblpot': ['pev'],
'rlds': ['msdwlwrf']}

if var == 'pr':

	lat, lon, cru_djf = import_obs(dict_var[var][0], domain, 'CRU', 'DJF')
	lat, lon, cru_mam = import_obs(dict_var[var][0], domain, 'CRU', 'MAM')
	lat, lon, cru_jja = import_obs(dict_var[var][0], domain, 'CRU', 'JJA')
	lat, lon, cru_son = import_obs(dict_var[var][0], domain, 'CRU', 'SON')

	lat, lon, cpc_djf = import_obs(dict_var[var][1], domain, 'CPC', 'DJF')
	lat, lon, cpc_mam = import_obs(dict_var[var][1], domain, 'CPC', 'MAM')
	lat, lon, cpc_jja = import_obs(dict_var[var][1], domain, 'CPC', 'JJA')
	lat, lon, cpc_son = import_obs(dict_var[var][1], domain, 'CPC', 'SON')

	lat, lon, cmorph_djf = import_obs(dict_var[var][2], domain, 'CMORPH', 'DJF')
	lat, lon, cmorph_mam = import_obs(dict_var[var][2], domain, 'CMORPH', 'MAM')
	lat, lon, cmorph_jja = import_obs(dict_var[var][2], domain, 'CMORPH', 'JJA')
	lat, lon, cmorph_son = import_obs(dict_var[var][2], domain, 'CMORPH', 'SON')

	lat, lon, mswep_djf = import_obs(dict_var[var][3], domain, 'MSWEP', 'DJF')
	lat, lon, mswep_mam = import_obs(dict_var[var][3], domain, 'MSWEP', 'MAM')
	lat, lon, mswep_jja = import_obs(dict_var[var][3], domain, 'MSWEP', 'JJA')
	lat, lon, mswep_son = import_obs(dict_var[var][3], domain, 'MSWEP', 'SON')

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

	mbe_djf_regcm_mswep = compute_mbe(regcm_djf, mswep_djf)
	mbe_mam_regcm_mswep = compute_mbe(regcm_mam, mswep_mam)
	mbe_jja_regcm_mswep = compute_mbe(regcm_jja, mswep_jja)
	mbe_son_regcm_mswep = compute_mbe(regcm_son, mswep_son)
		
	mbe_djf_regcm_era5 = compute_mbe(regcm_djf, era5_djf)
	mbe_mam_regcm_era5 = compute_mbe(regcm_mam, era5_mam)
	mbe_jja_regcm_era5 = compute_mbe(regcm_jja, era5_jja)
	mbe_son_regcm_era5 = compute_mbe(regcm_son, era5_son)	
	
elif var == 'tas' or var == 'clt':
	lat, lon, cru_djf = import_obs(dict_var[var][0], domain, 'CRU', 'DJF')
	lat, lon, cru_mam = import_obs(dict_var[var][0], domain, 'CRU', 'MAM')
	lat, lon, cru_jja = import_obs(dict_var[var][0], domain, 'CRU', 'JJA')
	lat, lon, cru_son = import_obs(dict_var[var][0], domain, 'CRU', 'SON')

	lat, lon, era5_djf = import_obs(dict_var[var][1], domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_obs(dict_var[var][1], domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_obs(dict_var[var][1], domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_obs(dict_var[var][1], domain, 'ERA5', 'SON')

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

elif var == 'tasmax' or var == 'tasmin':
	lat, lon, cru_djf = import_obs(dict_var[var][0], domain, 'CRU', 'DJF')
	lat, lon, cru_mam = import_obs(dict_var[var][0], domain, 'CRU', 'MAM')
	lat, lon, cru_jja = import_obs(dict_var[var][0], domain, 'CRU', 'JJA')
	lat, lon, cru_son = import_obs(dict_var[var][0], domain, 'CRU', 'SON')

	lat, lon, cpc_djf = import_obs(dict_var[var][1], domain, 'CPC', 'DJF')
	lat, lon, cpc_mam = import_obs(dict_var[var][1], domain, 'CPC', 'MAM')
	lat, lon, cpc_jja = import_obs(dict_var[var][1], domain, 'CPC', 'JJA')
	lat, lon, cpc_son = import_obs(dict_var[var][1], domain, 'CPC', 'SON')

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

	mbe_djf_regcm_cpc = compute_mbe(regcm_djf, cpc_djf)
	mbe_mam_regcm_cpc = compute_mbe(regcm_mam, cpc_mam)
	mbe_jja_regcm_cpc = compute_mbe(regcm_jja, cpc_jja)
	mbe_son_regcm_cpc = compute_mbe(regcm_son, cpc_son)
		
	mbe_djf_regcm_era5 = compute_mbe(regcm_djf, era5_djf)
	mbe_mam_regcm_era5 = compute_mbe(regcm_mam, era5_mam)
	mbe_jja_regcm_era5 = compute_mbe(regcm_jja, era5_jja)
	mbe_son_regcm_era5 = compute_mbe(regcm_son, era5_son)

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
def configure_subplot(ax):

	lon_min = np.round(np.min(lon), 1)
	lon_max = np.round(np.max(lon), 1)
	lat_min = np.round(np.min(lat), 1)
	lat_max = np.round(np.max(lat), 1)
	ax.set_extent([np.min(lon), np.max(lon), np.min(lat), np.max(lat)], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(lon_min,lon_max,10), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(lat_min,lat_max,5), crs=ccrs.PlateCarree())

	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(6)

	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.add_feature(cfeat.BORDERS, linewidth=0.5)
	ax.coastlines(linewidth=0.5)

	if var == 'evspsblpot':
		ax.add_feature(cfeat.OCEAN, facecolor='white', zorder=1) 


dict_plot = {'pr': ['Bias of  precipitation (mm d$^-$$^1$)', np.arange(-10, 11, 1), cm.BrBG],
'tas': ['Bias of air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'tasmax': ['Bias of maximum air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'tasmin': ['Bias of minimum air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'clt': ['Bias of total cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'cll': ['Bias of low cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'clm': ['Bias of medium cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'clh': ['Bias of high cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'evspsblpot': ['Bias of potential evapotranspiration (mm d$^-$$^1$)', np.arange(-5, 5.5, 0.5), cm.bwr],
'rlds': ['Bias of surface downwelling longwave radiation (W mm$^-$$^2$)', np.arange(-60, 55, 5), cm.RdBu_r]}

if var == 'pr':
	fig, axes = plt.subplots(4, 5, figsize=(16, 8), subplot_kw={'projection': ccrs.PlateCarree()})

	ax1 = axes[0, 0]
	plt_map = ax1.contourf(lon, lat, mbe_djf_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax1.set_title(u'(a) RegCM5-CRU DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax1)

	ax2 = axes[0, 1]
	plt_map = ax2.contourf(lon, lat, mbe_djf_regcm_cpc, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax2.set_title(u'(b) RegCM5-CPC DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax2)

	ax3 = axes[0, 2] 
	plt_map = ax3.contourf(lon, lat, mbe_djf_regcm_cmorph, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax3.set_title(u'(c) RegCM5-CMORPH DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax3)

	ax4 = axes[0, 3]
	plt_map = ax4.contourf(lon, lat, mbe_djf_regcm_mswep, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax4.set_title(u'(d) RegCM5-MSWEP DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax4)

	ax5 = axes[0, 4]
	plt_map = ax5.contourf(lon, lat, mbe_djf_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax5.set_title(u'(e) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax5)

	ax6 = axes[1, 0]
	plt_map = ax6.contourf(lon, lat, mbe_mam_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax6.set_title(u'(f) RegCM5-CRU MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax6)

	ax7 = axes[1, 1]
	plt_map = ax7.contourf(lon, lat, mbe_mam_regcm_cpc, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax7.set_title(u'(g) RegCM5-CPC MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax7)

	ax8 = axes[1, 2] 
	plt_map = ax8.contourf(lon, lat, mbe_mam_regcm_cmorph, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax8.set_title(u'(h) RegCM5-CMORPH MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax8)

	ax9 = axes[1, 3]
	plt_map = ax9.contourf(lon, lat, mbe_mam_regcm_mswep, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax9.set_title(u'(i) RegCM5-MSWEP MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax9)

	ax10 = axes[1, 4]
	plt_map = ax10.contourf(lon, lat, mbe_mam_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax10.set_title(u'(j) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax10)

	ax11 = axes[2, 0]
	plt_map = ax11.contourf(lon, lat, mbe_jja_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax11.set_title(u'(k) RegCM5-CRU JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax11)

	ax12 = axes[2, 1]
	plt_map = ax12.contourf(lon, lat, mbe_jja_regcm_cpc, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax12.set_title(u'(l) RegCM5-CPC JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax12)

	ax13 = axes[2, 2] 
	plt_map = ax13.contourf(lon, lat, mbe_jja_regcm_cmorph, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax13.set_title(u'(m) RegCM5-CMORPH JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax13)

	ax14 = axes[2, 3]
	plt_map = ax14.contourf(lon, lat, mbe_jja_regcm_mswep, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax14.set_title(u'(n) RegCM5-MSWEP JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax14)

	ax15 = axes[2, 4]
	plt_map = ax15.contourf(lon, lat, mbe_jja_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax15.set_title(u'(o) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax15)

	ax16 = axes[3, 0]
	plt_map = ax16.contourf(lon, lat, mbe_son_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax16.set_title(u'(p) RegCM5-CRU SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax16)

	ax17 = axes[3, 1]
	plt_map = ax17.contourf(lon, lat, mbe_son_regcm_cpc, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax17.set_title(u'(q) RegCM5-CPC SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax17)

	ax18 = axes[3, 2]
	plt_map = ax18.contourf(lon, lat, mbe_son_regcm_cmorph, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax18.set_title(u'(r) RegCM5-CMORPH SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax18)

	ax19 = axes[3, 3] 
	plt_map = ax19.contourf(lon, lat, mbe_son_regcm_mswep, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax19.set_title(u'(s) RegCM5-MSWEP SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax19)

	ax20 = axes[3, 4]
	plt_map = ax20.contourf(lon, lat, mbe_son_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax20.set_title(u'(t) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax20)

elif var == 'tas' or var == 'clt':
	fig, axes = plt.subplots(4, 2, figsize=(6, 8), subplot_kw={'projection': ccrs.PlateCarree()})

	ax1 = axes[0, 0]
	plt_map = ax1.contourf(lon, lat, mbe_djf_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax1.set_title(u'(a) RegCM5-CRU DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax1)

	ax2 = axes[0, 1]
	plt_map = ax2.contourf(lon, lat, mbe_djf_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax2.set_title(u'(b) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax2)

	ax3 = axes[1, 0]
	plt_map = ax3.contourf(lon, lat, mbe_mam_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax3.set_title(u'(c) RegCM5-CRU MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax3)

	ax4 = axes[1, 1]
	plt_map = ax4.contourf(lon, lat, mbe_mam_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax4.set_title(u'(d) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax4)

	ax5 = axes[2, 0]
	plt_map = ax5.contourf(lon, lat, mbe_jja_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax5.set_title(u'(e) RegCM5-CRU JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax5)

	ax6 = axes[2, 1]
	plt_map = ax6.contourf(lon, lat, mbe_jja_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax6.set_title(u'(f) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax6)

	ax7 = axes[3, 0]
	plt_map = ax7.contourf(lon, lat, mbe_son_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax7.set_title(u'(g) RegCM5-CRU SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax7)

	ax8 = axes[3, 1]
	plt_map = ax8.contourf(lon, lat, mbe_son_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax8.set_title(u'(h) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax8)

elif var == 'tasmax' or var == 'tasmin':
	fig, axes = plt.subplots(4, 3, figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})

	ax1 = axes[0, 0]
	plt_map = ax1.contourf(lon, lat, mbe_djf_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax1.set_title(u'(a) RegCM5-CRU DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax1)

	ax2 = axes[0, 1]
	plt_map = ax2.contourf(lon, lat, mbe_djf_regcm_cpc, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax2.set_title(u'(b) RegCM5-CPC DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax2)

	ax3 = axes[0, 2]
	plt_map = ax3.contourf(lon, lat, mbe_djf_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax3.set_title(u'(c) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax3)

	ax4 = axes[1, 0]
	plt_map = ax4.contourf(lon, lat, mbe_mam_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax4.set_title(u'(d) RegCM5-CRU MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax4)

	ax5 = axes[1, 1]
	plt_map = ax5.contourf(lon, lat, mbe_mam_regcm_cpc, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax5.set_title(u'(e) RegCM5-CPC MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax5)

	ax6 = axes[1, 2]
	plt_map = ax6.contourf(lon, lat, mbe_mam_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax6.set_title(u'(f) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax6)

	ax7 = axes[2, 0]
	plt_map = ax7.contourf(lon, lat, mbe_jja_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax7.set_title(u'(g) RegCM5-CRU JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax7)

	ax8 = axes[2, 1]
	plt_map = ax8.contourf(lon, lat, mbe_jja_regcm_cpc, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax8.set_title(u'(h) RegCM5-CPC JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax8)

	ax9 = axes[2, 2]
	plt_map = ax9.contourf(lon, lat, mbe_jja_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax9.set_title(u'(i) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax9)

	ax10 = axes[3, 0]
	plt_map = ax10.contourf(lon, lat, mbe_son_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax10.set_title(u'(j) RegCM5-CRU SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax10)

	ax11 = axes[3, 1]
	plt_map = ax11.contourf(lon, lat, mbe_son_regcm_cpc, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax11.set_title(u'(k) RegCM5-CPC SON', loc='left', fontsize=font_size, fontweight='bold')
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
cbar = fig.colorbar(plt_map, ax=fig.axes, pad=0.02, aspect=50)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/figs/evaluate'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()



	
	

	
	
