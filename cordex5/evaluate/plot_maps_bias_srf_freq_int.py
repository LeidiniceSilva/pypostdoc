# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot bias maps"

import os
import netCDF4
import argparse
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat

from cartopy import config
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from import_climate_tools import compute_mbe

parser = argparse.ArgumentParser()
parser.add_argument('--var', required=True, help='Variable name')
parser.add_argument('--stats', required=True, help='Statistic')
parser.add_argument('--freq', required=True, help='Frequency')
parser.add_argument('--domain', required=True, help='Domain name')
parser.add_argument('--idt', required=True, help='Initial year')
parser.add_argument('--fdt', required=True, help='Final year')
args = parser.parse_args()

var = args.var
stats = args.stats
freq = args.freq
domain = args.domain
idt = args.idt
fdt = args.fdt
font_size = 6

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5'


def import_obs(param, domain, dataset, season):

	if freq == 'hourly':
		dt = '1hr_{0}_{1}-{2}_th0.5'.format(season, idt, fdt)
	else:
		dt = '{0}_{1}-{2}'.format(season, idt, fdt)

	arq   = '{0}/postproc/evaluate/obs/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, stats, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,0,:,:]
	
	return lat, lon, mean
	
		
def import_rcm(param, domain, dataset, season):

	if freq == 'hourly':
		dt = '1hr_{0}_{1}-{2}_th0.5'.format(season, idt, fdt)
	else:
		dt = '{0}_{1}-{2}'.format(season, idt, fdt)

	arq   = '{0}/postproc/evaluate/rcm/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, stats, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,0,:,:]
	
	return lat, lon, mean
	

# Import model and obs dataset
if freq == 'hourly':
	lat, lon, cmorph_djf = import_obs('cmorph', domain, 'CMORPH', 'DJF')
	lat, lon, cmorph_mam = import_obs('cmorph', domain, 'CMORPH', 'MAM')
	lat, lon, cmorph_jja = import_obs('cmorph', domain, 'CMORPH', 'JJA')
	lat, lon, cmorph_son = import_obs('cmorph', domain, 'CMORPH', 'SON')

	lat, lon, era5_djf = import_obs('tp', domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_obs('tp', domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_obs('tp', domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_obs('tp', domain, 'ERA5', 'SON')
	
	lat, lon, regcm_djf = import_rcm('pr', domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_rcm('pr', domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_rcm('pr', domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_rcm('pr', domain, 'RegCM5', 'SON')

	mbe_regcm_cmorph_djf = regcm_djf - cmorph_djf
	mbe_regcm_cmorph_mam = regcm_mam - cmorph_mam
	mbe_regcm_cmorph_jja = regcm_jja - cmorph_jja
	mbe_regcm_cmorph_son = regcm_son - cmorph_son
	
	mbe_regcm_era5_djf = regcm_djf - era5_djf
	mbe_regcm_era5_mam = regcm_mam - era5_mam
	mbe_regcm_era5_jja = regcm_jja - era5_jja
	mbe_regcm_era5_son = regcm_son - era5_son
else:
	lat, lon, cpc_djf = import_obs('precip', domain, 'CPC', 'DJF')
	lat, lon, cpc_mam = import_obs('precip', domain, 'CPC', 'MAM')
	lat, lon, cpc_jja = import_obs('precip', domain, 'CPC', 'JJA')
	lat, lon, cpc_son = import_obs('precip', domain, 'CPC', 'SON')

	lat, lon, cmorph_djf = import_obs('cmorph', domain, 'CMORPH', 'DJF')
	lat, lon, cmorph_mam = import_obs('cmorph', domain, 'CMORPH', 'MAM')
	lat, lon, cmorph_jja = import_obs('cmorph', domain, 'CMORPH', 'JJA')
	lat, lon, cmorph_son = import_obs('cmorph', domain, 'CMORPH', 'SON')

	lat, lon, mswep_djf = import_obs('precipitation', domain, 'MSWEP', 'DJF')
	lat, lon, mswep_mam = import_obs('precipitation', domain, 'MSWEP', 'MAM')
	lat, lon, mswep_jja = import_obs('precipitation', domain, 'MSWEP', 'JJA')
	lat, lon, mswep_son = import_obs('precipitation', domain, 'MSWEP', 'SON')
	
	lat, lon, era5_djf = import_obs('tp', domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_obs('tp', domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_obs('tp', domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_obs('tp', domain, 'ERA5', 'SON')
	
	lat, lon, regcm_djf = import_rcm('pr', domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_rcm('pr', domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_rcm('pr', domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_rcm('pr', domain, 'RegCM5', 'SON')

	mbe_regcm_cpc_djf = regcm_djf - cpc_djf
	mbe_regcm_cpc_mam = regcm_mam - cpc_mam
	mbe_regcm_cpc_jja = regcm_jja - cpc_jja
	mbe_regcm_cpc_son = regcm_son - cpc_son

	mbe_regcm_cmorph_djf = regcm_djf - cmorph_djf
	mbe_regcm_cmorph_mam = regcm_mam - cmorph_mam
	mbe_regcm_cmorph_jja = regcm_jja - cmorph_jja
	mbe_regcm_cmorph_son = regcm_son - cmorph_son
	
	mbe_regcm_mswep_djf = regcm_djf - mswep_djf
	mbe_regcm_mswep_mam = regcm_mam - mswep_mam
	mbe_regcm_mswep_jja = regcm_jja - mswep_jja
	mbe_regcm_mswep_son = regcm_son - mswep_son
	
	mbe_regcm_era5_djf = regcm_djf - era5_djf
	mbe_regcm_era5_mam = regcm_mam - era5_mam
	mbe_regcm_era5_jja = regcm_jja - era5_jja
	mbe_regcm_era5_son = regcm_son - era5_son

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

if stats == 'freq':
	if freq == 'hourly':
		levs = np.arange(-22, 23, 1)
		legend = 'Frequency (%)'
		dt = '1hr_{0}-{1}_th0.5'.format(idt, fdt)
	else:
		levs = np.arange(-44, 46, 2)
		legend = 'Frequency (%)'
		dt = '{0}-{1}'.format(idt, fdt)
else:
	if freq == 'hourly':
		levs = np.arange(-6, 6.5, 0.5)
		legend = 'Intensity (mm h$^-$$^1$)'
		dt = '1hr_{0}-{1}_th0.5'.format(idt, fdt)
	else:
		levs = np.arange(-22, 23, 1)
		legend = 'Intensity (mm d$^-$$^1$)'
		dt = '{0}-{1}'.format(idt, fdt)

if freq == 'hourly':
	fig, axes = plt.subplots(4, 2, figsize=(6, 8), subplot_kw={'projection': ccrs.PlateCarree()})

	ax1 = axes[0, 0]
	plt_map = ax1.contourf(lon, lat, mbe_regcm_cmorph_djf, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax1.set_title(u'(a) RegCM5-CMORPH DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax1)

	ax2 = axes[0, 1]
	plt_map = ax2.contourf(lon, lat, mbe_regcm_era5_djf, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax2.set_title(u'(b) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax2)

	ax3 = axes[1, 0]
	plt_map = ax3.contourf(lon, lat, mbe_regcm_cmorph_mam, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax3.set_title(u'(c) RegCM5-CMORPH MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax3)

	ax4 = axes[1, 1]
	plt_map = ax4.contourf(lon, lat, mbe_regcm_era5_mam, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax4.set_title(u'(d) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax4)

	ax5 = axes[2, 0]
	plt_map = ax5.contourf(lon, lat, mbe_regcm_cmorph_jja, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax5.set_title(u'(e) RegCM5-CMORPH JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax5)

	ax6 = axes[2, 1]
	plt_map = ax6.contourf(lon, lat, mbe_regcm_era5_jja, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax6.set_title(u'(f) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax6)

	ax7 = axes[3, 0]
	plt_map = ax7.contourf(lon, lat, mbe_regcm_cmorph_son, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax7.set_title(u'(g) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax7)

	ax8 = axes[3, 1]
	plt_map = ax8.contourf(lon, lat, mbe_regcm_era5_son, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax8.set_title(u'(h) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax8)	
else:
	fig, axes = plt.subplots(4, 4, figsize=(14, 8), subplot_kw={'projection': ccrs.PlateCarree()})

	ax1 = axes[0, 0]
	plt_map = ax1.contourf(lon, lat, mbe_regcm_cpc_djf, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax1.set_title(u'(a) RegCM5-CPC DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax1)

	ax2 = axes[0, 1]
	plt_map = ax2.contourf(lon, lat, mbe_regcm_cmorph_djf, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax2.set_title(u'(b) RegCM5-CMORPH DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax2)

	ax3 = axes[0, 2]
	plt_map = ax3.contourf(lon, lat, mbe_regcm_mswep_djf, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax3.set_title(u'(c) RegCM5-MSWEP DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax3)

	ax4 = axes[0, 3]
	plt_map = ax4.contourf(lon, lat, mbe_regcm_era5_djf, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax4.set_title(u'(d) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax4)

	ax5 = axes[1, 0]
	plt_map = ax5.contourf(lon, lat, mbe_regcm_cpc_mam, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax5.set_title(u'(e) RegCM5-CPC MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax5)

	ax6 = axes[1, 1]
	plt_map = ax6.contourf(lon, lat, mbe_regcm_cmorph_mam, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax6.set_title(u'(f) RegCM5-CMORPH MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax6)

	ax7 = axes[1, 2]
	plt_map = ax7.contourf(lon, lat, mbe_regcm_mswep_mam, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax7.set_title(u'(g) RegCM5-MSWEP MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax7)

	ax8 = axes[1, 3]
	plt_map = ax8.contourf(lon, lat, mbe_regcm_era5_mam, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax8.set_title(u'(h) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax8)	

	ax9 = axes[2, 0]
	plt_map = ax9.contourf(lon, lat, mbe_regcm_cpc_jja, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax9.set_title(u'(i) RegCM5-CPC JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax9)

	ax10 = axes[2, 1]
	plt_map = ax10.contourf(lon, lat, mbe_regcm_cmorph_jja, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax10.set_title(u'(j) RegCM5-CMORPH JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax10)

	ax11 = axes[2, 2]
	plt_map = ax11.contourf(lon, lat, mbe_regcm_mswep_jja, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax11.set_title(u'(g) RegCM5-MSWEP JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax11)

	ax12 = axes[2, 3]
	plt_map = ax12.contourf(lon, lat, mbe_regcm_era5_jja, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax12.set_title(u'(k) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax12)	

	ax13 = axes[3, 0]
	plt_map = ax13.contourf(lon, lat, mbe_regcm_cpc_mam, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax13.set_title(u'(l) RegCM5-CPC MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax13)

	ax14 = axes[3, 1]
	plt_map = ax14.contourf(lon, lat, mbe_regcm_cmorph_mam, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax14.set_title(u'(m) RegCM5-CMORPH MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax14)

	ax15 = axes[3, 2]
	plt_map = ax15.contourf(lon, lat, mbe_regcm_mswep_mam, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax15.set_title(u'(n) RegCM5-MSWEP MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax15)

	ax16 = axes[3, 3]
	plt_map = ax16.contourf(lon, lat, mbe_regcm_era5_mam, transform=ccrs.PlateCarree(), levels=levs, cmap=cm.BrBG, extend='neither') 
	ax16.set_title(u'(o) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax16)	

# Set colobar
cbar = fig.colorbar(plt_map, ax=fig.axes, pad=0.02, aspect=50)
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/figs/evaluate'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_{2}_RegCM5_{3}.png'.format(var, stats, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
