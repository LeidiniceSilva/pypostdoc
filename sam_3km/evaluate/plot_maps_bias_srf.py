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

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii

var = 'rsnl'
domain = 'SAM-3km'
idt, fdt = '2018', '2021'
dt = '{0}-{1}'.format(idt, fdt)

path = '/leonardo/home/userexternal/mdasilva/leonardo_work'
	
	
def import_obs(param, domain, dataset, season):

	arq   = '{0}/SAM-3km/postproc/evaluate/obs/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]

	if param == 'pev' or param == 'msnlwrf':
		mean = var[:][0,:,:]*(-1)
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
	else:
		mean = var[:][0,:,:]

	return lat, lon, mean
	

def configure_subplot(ax):

	ax.set_xticks(np.arange(-78,-32,12), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(-38,-6,6), crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.tick_params(axis='x', labelsize=6, labelcolor='black')
	ax.tick_params(axis='y', labelsize=6, labelcolor='black')
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.coastlines()


# Import model and obs dataset
dict_var = {
'pr': ['pre', 'pre', 'precip', 'sat_gauge_precip', 'tp'],
'tas': ['tmp', 'tmp', 't2m'],
'clt': ['tcc'],
'cll': ['lcc'],
'clm': ['mcc'],
'clh': ['hcc'],
'evspsblpot': ['pev'],
'rsnl': ['msnlwrf'],
'rsns': ['msnswrf']
}

if var == 'pr':
	lat, lon, cru_djf = import_obs(dict_var[var][1], domain, 'CRU', 'DJF')
	lat, lon, cru_mam = import_obs(dict_var[var][1], domain, 'CRU', 'MAM')
	lat, lon, cru_jja = import_obs(dict_var[var][1], domain, 'CRU', 'JJA')
	lat, lon, cru_son = import_obs(dict_var[var][1], domain, 'CRU', 'SON')

	lat, lon, cpc_djf = import_obs(dict_var[var][2], domain, 'CPC', 'DJF')
	lat, lon, cpc_mam = import_obs(dict_var[var][2], domain, 'CPC', 'MAM')
	lat, lon, cpc_jja = import_obs(dict_var[var][2], domain, 'CPC', 'JJA')
	lat, lon, cpc_son = import_obs(dict_var[var][2], domain, 'CPC', 'SON')

	lat, lon, gpcp_djf = import_obs(dict_var[var][3], domain, 'GPCP', 'DJF')
	lat, lon, gpcp_mam = import_obs(dict_var[var][3], domain, 'GPCP', 'MAM')
	lat, lon, gpcp_jja = import_obs(dict_var[var][3], domain, 'GPCP', 'JJA')
	lat, lon, gpcp_son = import_obs(dict_var[var][3], domain, 'GPCP', 'SON')

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
	
	mbe_djf_regcm_gpcp = compute_mbe(regcm_djf, gpcp_djf)
	mbe_mam_regcm_gpcp = compute_mbe(regcm_mam, gpcp_mam)
	mbe_jja_regcm_gpcp = compute_mbe(regcm_jja, gpcp_jja)
	mbe_son_regcm_gpcp = compute_mbe(regcm_son, gpcp_son)
		
	mbe_djf_regcm_era5 = compute_mbe(regcm_djf, era5_djf)
	mbe_mam_regcm_era5 = compute_mbe(regcm_mam, era5_mam)
	mbe_jja_regcm_era5 = compute_mbe(regcm_jja, era5_jja)
	mbe_son_regcm_era5 = compute_mbe(regcm_son, era5_son)		

elif var == 'tas':
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
font_size = 6

dict_plot = {
'pr': ['Bias of  precipitation (mm d$^-$$^1$)', np.arange(-10, 11, 1), cm.BrBG],
'tas': ['Bias of air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'clt': ['Bias of total cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'cll': ['Bias of low cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'clm': ['Bias of medium cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'clh': ['Bias of high cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'evspsblpot': ['Potential evaporation (mm d$^-$$^1$)', np.arange(-10, 11, 1), cm.bwr],
'rsnl': ['Bias of surface net upward longwave flux (W mm$^-$$^2$)', np.arange(-80, 90, 10), cm.RdBu_r],
'rsns': ['Bias of surface net downward shortwave flux (W mm$^-$$^2$)', np.arange(-80, 90, 10), cm.RdBu_r]
}
	
if var == 'pr':
	fig, axes = plt.subplots(4, 4, figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})

	ax1 = axes[0, 0]
	plt_map = ax1.contourf(lon, lat, mbe_djf_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax1.set_title(u'(a) RegCM5-CRU DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax1)

	ax2 = axes[0, 1]
	plt_map = ax2.contourf(lon, lat, mbe_djf_regcm_cpc, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax2.set_title(u'(b) RegCM5-CPC DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax2)

	ax3 = axes[0, 2] 
	plt_map = ax3.contourf(lon, lat, mbe_djf_regcm_gpcp, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax3.set_title(u'(c) RegCM5-GPCP DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax3)

	ax4 = axes[0, 3]
	plt_map = ax4.contourf(lon, lat, mbe_djf_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax4.set_title(u'(d) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax4)

	ax5 = axes[1, 0]
	plt_map = ax5.contourf(lon, lat, mbe_mam_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax5.set_title(u'(e) RegCM5-CRU MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax5)

	ax6 = axes[1, 1]
	plt_map = ax6.contourf(lon, lat, mbe_mam_regcm_cpc, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax6.set_title(u'(f) RegCM5-CPC MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax6)

	ax7 = axes[1, 2] 
	plt_map = ax7.contourf(lon, lat, mbe_mam_regcm_gpcp, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax7.set_title(u'(g) RegCM5-GPCP MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax7)

	ax8 = axes[1, 3]
	plt_map = ax8.contourf(lon, lat, mbe_mam_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax8.set_title(u'(h) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax8)

	ax9 = axes[2, 0]
	plt_map = ax9.contourf(lon, lat, mbe_jja_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax9.set_title(u'(i) RegCM5-CRU JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax9)

	ax10 = axes[2, 1]
	plt_map = ax10.contourf(lon, lat, mbe_jja_regcm_cpc, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax10.set_title(u'(j) RegCM5-CPC JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax10)

	ax11 = axes[2, 2] 
	plt_map = ax11.contourf(lon, lat, mbe_jja_regcm_gpcp, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax11.set_title(u'(k) RegCM5-GPCP JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax11)

	ax12 = axes[2, 3]
	plt_map = ax12.contourf(lon, lat, mbe_jja_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax12.set_title(u'(l) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax12)

	ax13 = axes[3, 0]
	plt_map = ax13.contourf(lon, lat, mbe_son_regcm_cru, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax13.set_title(u'(m) RegCM5-CRU SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax13)

	ax14 = axes[3, 1]
	plt_map = ax14.contourf(lon, lat, mbe_son_regcm_cpc, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax14.set_title(u'(n) RegCM5-CPC SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax14)

	ax15 = axes[3, 2] 
	plt_map = ax15.contourf(lon, lat, mbe_son_regcm_gpcp, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax15.set_title(u'(o) RegCM5-GPCP SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax15)

	ax16 = axes[3, 3]
	plt_map = ax16.contourf(lon, lat, mbe_son_regcm_era5, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	ax16.set_title(u'(p) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax16)
	
elif var == 'tas':
	fig, axes = plt.subplots(4, 2, figsize=(6, 10), subplot_kw={'projection': ccrs.PlateCarree()})

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

else:
	fig, axes = plt.subplots(4, 1, figsize=(5, 10), subplot_kw={'projection': ccrs.PlateCarree()})
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
if var == 'clt' or var == 'cll' or var == 'clm' or var == 'clh' or var == 'evspsblpot' or var == 'rsnl' or var == 'rsns':
	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.84, 0.3, 0.03, 0.4]))
	cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)
else:
	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.94, 0.3, 0.015, 0.4]))
	cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
exit()
