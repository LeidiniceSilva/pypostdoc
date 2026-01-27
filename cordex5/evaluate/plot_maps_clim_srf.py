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
from matplotlib.colors import LinearSegmentedColormap
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

var = 'tas'
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
		mean = np.abs(var[:][0,:,:])
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
'rlds': ['msdwlwrf'],
'cape': ['cape'],
'cin': ['cin']}

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

	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(6)

	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.add_feature(cfeat.BORDERS, linewidth=0.5)
	ax.coastlines(linewidth=0.5)

	if var == 'evspsblpot':
		ax.add_feature(cfeat.OCEAN, facecolor='white', zorder=1) 

pr_colormap=['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']

colors = [
    (1, 1, 1),      # White for the lowest values
    (0.8, 0.9, 1),  # Light blue
    (0.6, 0.7, 1),  # Medium blue
    (0.4, 0.4, 1),  # Dark blue
    (1, 0.8, 0.6),  # Light orange
    (1, 0.5, 0.2),  # Orange
    (1, 0.2, 0),    # Red
    (0.6, 0, 0.6),  # Dark red/purple
    (0.4, 0, 0.4)   # Dark purple
]
cape_colormap = LinearSegmentedColormap.from_list('cape_colormap', colors, N=256)

dict_plot = {'pr': ['Precipitation (mm d$^-$$^1$)', np.arange(0, 18, 1), matplotlib.colors.ListedColormap(pr_colormap)],
'tas': ['Air temperature (°C)', np.arange(-6, 34, 2), cm.jet],
'tasmax': ['Maximum air temperature (°C)', np.arange(-6, 40, 2), cm.jet],
'tasmin': ['Minimum air temperature (°C)', np.arange(-10, 36, 2), cm.jet],
'clt': ['Total cloud cover (0-1)', np.arange(0, 1.05, 0.05), cm.Greys],
'cll': ['Low cloud cover (0-1)', np.arange(0, 1.05, 0.05), cm.Greys],
'clm': ['Medium cloud cover (0-1)', np.arange(0, 1.05, 0.05), cm.Greys],
'clh': ['High cloud cover (0-1)', np.arange(0, 1.05, 0.05), cm.Greys],
'evspsblpot': ['Potential evapotranspiration (mm d$^-$$^1$)', np.arange(0, 10.5, 0.5), cm.jet],
'rlds': ['Surface downwelling longwave radiation (W mm$^-$$^2$)', np.arange(150, 440, 10), cm.rainbow],
'cape': ['Convective available potential energy (J kg$^-$$^1$)', np.arange(0, 1850, 50), cape_colormap],
'cin': ['Convective inhibition (J kg$^-$$^1$)', np.arange(0, 430, 10), cape_colormap]}

if var == 'pr':
	fig, axes = plt.subplots(4, 6, figsize=(20, 8), subplot_kw={'projection': ccrs.PlateCarree()})

	ax1 = axes[0, 0]
	plt_map = ax1.contourf(lon, lat, regcm_djf, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax1.set_title(u'(a) RegCM5 DJF', loc='left', fontsize=font_size, fontweight='bold')
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
	plt_map = ax5.contourf(lon, lat, mswep_djf, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax5.set_title(u'(e) MSWEP DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax5)

	ax6 = axes[0, 5]
	plt_map = ax6.contourf(lon, lat, era5_djf, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax6.set_title(u'(f) ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax6)

	ax7 = axes[1, 0]
	plt_map = ax7.contourf(lon, lat, regcm_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax7.set_title(u'(g) RegCM5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax7)

	ax8 = axes[1, 1]
	plt_map = ax8.contourf(lon, lat, cru_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax8.set_title(u'(h) CRU MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax8)

	ax9 = axes[1, 2] 
	plt_map = ax9.contourf(lon, lat, cpc_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax9.set_title(u'(i) CPC MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax9)

	ax10 = axes[1, 3]
	plt_map = ax10.contourf(lon, lat, cmorph_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax10.set_title(u'(j) CMORPH MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax10)

	ax11 = axes[1, 4]
	plt_map = ax11.contourf(lon, lat, mswep_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax11.set_title(u'(k) MSWEP MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax11)

	ax12 = axes[1, 5]
	plt_map = ax12.contourf(lon, lat, era5_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax12.set_title(u'(l) ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax12)

	ax13 = axes[2, 0]
	plt_map = ax13.contourf(lon, lat, regcm_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax13.set_title(u'(m) RegCM5 JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax13)

	ax14 = axes[2, 1]
	plt_map = ax14.contourf(lon, lat, cru_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax14.set_title(u'(n) CRU JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax14)

	ax15 = axes[2, 2] 
	plt_map = ax15.contourf(lon, lat, cpc_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax15.set_title(u'(o) CPC JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax15)

	ax16 = axes[2, 3]
	plt_map = ax16.contourf(lon, lat, cmorph_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax16.set_title(u'(p) CMORPH JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax16)

	ax17 = axes[2, 4]
	plt_map = ax17.contourf(lon, lat, mswep_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax17.set_title(u'(q) MSWEP JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax17)

	ax18 = axes[2, 5]
	plt_map = ax18.contourf(lon, lat, era5_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax18.set_title(u'(r) ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax18)

	ax19 = axes[3, 0]
	plt_map = ax19.contourf(lon, lat, regcm_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax19.set_title(u'(s) RegCM5 SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax19)

	ax20 = axes[3, 1]
	plt_map = ax20.contourf(lon, lat, cru_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax20.set_title(u'(t) CRU SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax20)

	ax21 = axes[3, 2] 
	plt_map = ax21.contourf(lon, lat, cpc_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax21.set_title(u'(u) CPC SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax21)

	ax22 = axes[3, 3]
	plt_map = ax22.contourf(lon, lat, cmorph_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax22.set_title(u'(v) CMORPH SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax22)

	ax23 = axes[3, 4]
	plt_map = ax23.contourf(lon, lat, mswep_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax23.set_title(u'(w) MSWEP SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax23)

	ax24 = axes[3, 5]
	plt_map = ax24.contourf(lon, lat, era5_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax24.set_title(u'(x) ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax24)
	
elif var == 'tas' or var == 'clt':
	fig, axes = plt.subplots(4, 3, figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})

	ax1 = axes[0, 0]
	plt_map = ax1.contourf(lon, lat, regcm_djf, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax1.set_title(u'(a) RegCM5 DJF', loc='left', fontsize=font_size, fontweight='bold')
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
	ax4.set_title(u'(d) RegCM5 MAM', loc='left', fontsize=font_size, fontweight='bold')
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
	ax7.set_title(u'(g) RegCM5 JJA', loc='left', fontsize=font_size, fontweight='bold')
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
	ax10.set_title(u'(j) RegCM5 SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax10)

	ax11 = axes[3, 1]
	plt_map = ax11.contourf(lon, lat, cru_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax11.set_title(u'(k) CRU SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax11)

	ax12 = axes[3, 2]
	plt_map = ax12.contourf(lon, lat, era5_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax12.set_title(u'(l) ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax12)

elif var == 'tasmax' or var == 'tasmin':
	fig, axes = plt.subplots(4, 4, figsize=(14, 8), subplot_kw={'projection': ccrs.PlateCarree()})

	ax1 = axes[0, 0]
	plt_map = ax1.contourf(lon, lat, regcm_djf, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax1.set_title(u'(a) RegCM5 DJF', loc='left', fontsize=font_size, fontweight='bold')
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
	plt_map = ax4.contourf(lon, lat, era5_djf, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax4.set_title(u'(d) ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax4)

	ax5 = axes[1, 0]
	plt_map = ax5.contourf(lon, lat, regcm_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax5.set_title(u'(e) RegCM5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax5)

	ax6 = axes[1, 1]
	plt_map = ax6.contourf(lon, lat, cru_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax6.set_title(u'(f) CRU MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax6)

	ax7 = axes[1, 2]
	plt_map = ax7.contourf(lon, lat, cpc_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax7.set_title(u'(g) CPC MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax7)

	ax8 = axes[1, 3]
	plt_map = ax8.contourf(lon, lat, era5_mam, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax8.set_title(u'(h) ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax8)

	ax9 = axes[2, 0]
	plt_map = ax9.contourf(lon, lat, regcm_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax9.set_title(u'(i) RegCM5 JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax9)

	ax10 = axes[2, 1]
	plt_map = ax10.contourf(lon, lat, cru_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax10.set_title(u'(j) CRU JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax10)

	ax11 = axes[2, 2]
	plt_map = ax11.contourf(lon, lat, cpc_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax11.set_title(u'(k) CPC JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax11)

	ax12 = axes[2, 3]
	plt_map = ax12.contourf(lon, lat, era5_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax12.set_title(u'(l) ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax12)

	ax13 = axes[3, 0]
	plt_map = ax13.contourf(lon, lat, regcm_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax13.set_title(u'(m) RegCM5 SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax13)

	ax14 = axes[3, 1]
	plt_map = ax14.contourf(lon, lat, cru_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax14.set_title(u'(n) CRU SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax14)

	ax15 = axes[3, 2]
	plt_map = ax15.contourf(lon, lat, cpc_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax15.set_title(u'(o) CPC SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax15)

	ax16 = axes[3, 3]
	plt_map = ax16.contourf(lon, lat, era5_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax16.set_title(u'(p) ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax16)

else:
	fig, axes = plt.subplots(4, 2, figsize=(6, 8), subplot_kw={'projection': ccrs.PlateCarree()})

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
	ax5.set_title(u'(e) RegCM5 JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax5)

	ax6 = axes[2, 1]
	plt_map = ax6.contourf(lon, lat, era5_jja, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax6.set_title(u'(f) ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax6)

	ax7 = axes[3, 0]
	plt_map = ax7.contourf(lon, lat, regcm_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax7.set_title(u'(g) RegCM5 SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax7)

	ax8 = axes[3, 1]
	plt_map = ax8.contourf(lon, lat, era5_son, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax8.set_title(u'(h) ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax8)

# Set colobar
cbar = fig.colorbar(plt_map, ax=fig.axes, pad=0.02, aspect=50)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/figs/evaluate/rcm'.format(path)
name_out = 'pyplt_maps_clim_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


