# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Feb 02, 2026"
__description__ = "This script plot maps"

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
dt = '197001-02'
domain = 'SAM-12'
latlon = [-105, -16, -57, 18]

exp_i = 'RegCM5-ERA5_ICTP'
exp_ii = 'RegCM5-ERA5_USP'
exp_iii = 'RegCM5-ERA5_USP-dt60'
exp_iv = 'RegCM5-ERA5_USP-garoa'

dict_var = {'pr': ['tp'],
'tas': ['t2m']}

font_size = 6
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/{0}'.format(domain)


def import_obs(param, dataset):

	arq   = '{0}/postproc/obs/{1}_{2}_{3}_mon_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = np.nanmean(var[:][:,:,:], axis=0)
	print(mean.shape)

	return lat, lon, mean


def import_rcm_i(param, dataset):

	arq   = '{0}/postproc/rcm/{1}_{2}_{3}_mon_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = np.nanmean(var[:][:,:,:], axis=0)
	print(mean.shape)

	return lat, lon, mean


def import_rcm_ii(param, dataset):

	arq   = '{0}/postproc/rcm/{1}_{2}_{3}_mon_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = np.nanmean(var[:][:-1,0,:,:], axis=0)
	print(mean.shape)

	return lat, lon, mean


def configure_subplot(ax):

	ax.set_extent(latlon, crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(latlon[0], latlon[1], 20), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(latlon[2], latlon[3], 20), crs=ccrs.PlateCarree())

	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(6)

	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.add_feature(cfeat.BORDERS, linewidth=0.5)
	ax.coastlines(linewidth=0.5)


# Import model and obs dataset
lat, lon, obs_ = import_obs(dict_var[var][0], 'ERA5')
lat, lon, exp_i_ = import_rcm_i(var, exp_i)
lat, lon, exp_ii_ = import_rcm_ii(var, exp_ii)
lat, lon, exp_iii_ = import_rcm_ii(var, exp_iii)
lat, lon, exp_iv_ = import_rcm_ii(var, exp_iv)

exp_i_obs = exp_i_ - obs_
exp_ii_obs = exp_ii_ - obs_
exp_iii_obs = exp_iii_ - obs_
exp_iv_obs = exp_iv_ - obs_

# Plot figure   
fig, axes = plt.subplots(2, 5, figsize=(14, 5), subplot_kw={'projection': ccrs.PlateCarree()})
fig.delaxes(axes[1, 4])

pr_colormap=['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']

if var == 'pr':
	legend = 'Precipitation (mm d$^-$$^1$)'
	bmin = -10
	bmax = 10
	vmin = 0
	vmax = 18
	int = 1
	colorc = matplotlib.colors.ListedColormap(pr_colormap)
	colorb = cm.BrBG
else:
	legend = 'Air temperature (°C)'
	bmin = -10
	bmax = 12
	vmin = 0
	vmax = 40
	int = 2
	colorc = cm.jet
	colorb = cm.bwr

ax1 = axes[0, 0]
plt_map1 = ax1.contourf(lon, lat, obs_, transform=ccrs.PlateCarree(), levels=np.arange(vmin, vmax, int), cmap=colorc, extend='max') 
ax1.set_title(u'(a) ERA5 {0}'.format(dt), loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

ax2 = axes[0, 1]
plt_map2 = ax2.contourf(lon, lat, exp_i_, transform=ccrs.PlateCarree(), levels=np.arange(vmin, vmax, int), cmap=colorc, extend='max') 
ax2.set_title(u'(b) {0} {1}'.format(exp_i, dt), loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

ax3 = axes[0, 2] 
plt_map3 = ax3.contourf(lon, lat, exp_ii_, transform=ccrs.PlateCarree(), levels=np.arange(vmin, vmax, int), cmap=colorc, extend='max') 
ax3.set_title(u'(c) {0} {1}'.format(exp_ii, dt), loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

ax4 = axes[0, 3]
plt_map4 = ax4.contourf(lon, lat, exp_iii_, transform=ccrs.PlateCarree(), levels=np.arange(vmin, vmax, int), cmap=colorc, extend='max') 
ax4.set_title(u'(d) {0} {1}'.format(exp_iii, dt), loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax4)

ax5 = axes[0, 4] 
plt_map5 = ax5.contourf(lon, lat, exp_iv_, transform=ccrs.PlateCarree(), levels=np.arange(vmin, vmax, int), cmap=colorc, extend='max') 
ax5.set_title(u'(e) {0} {1}'.format(exp_iv, dt), loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax5)

ax6 = axes[1, 0]
plt_map6 = ax6.contourf(lon, lat, exp_i_obs, transform=ccrs.PlateCarree(), levels=np.arange(bmin, bmax, int), cmap=colorb, extend='both') 
ax6.set_title(u'(f) {0} - ERA5 {1}'.format(exp_i, dt), loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax6)

ax7 = axes[1, 1]
plt_map7 = ax7.contourf(lon, lat, exp_ii_obs, transform=ccrs.PlateCarree(), levels=np.arange(bmin, bmax, int), cmap=colorb, extend='both') 
ax7.set_title(u'(g) {0} - ERA5 {1}'.format(exp_ii, dt), loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax7)

ax8 = axes[1, 2]
plt_map8 = ax8.contourf(lon, lat, exp_iii_obs, transform=ccrs.PlateCarree(), levels=np.arange(bmin, bmax, int), cmap=colorb, extend='both') 
ax8.set_title(u'(h) {0} - ERA5 {1}'.format(exp_iii, dt), loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax8)

ax9 = axes[1, 3]
plt_map9 = ax9.contourf(lon, lat, exp_iv_obs, transform=ccrs.PlateCarree(), levels=np.arange(bmin, bmax, int), cmap=colorb, extend='both') 
ax9.set_title(u'(h) {0} - ERA5 {1}'.format(exp_iv, dt), loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax9)

# Set colobar
cbar1_ax = fig.add_axes([0.92, 0.55, 0.02, 0.35]) 
cbar1 = fig.colorbar(plt_map1, cax=cbar1_ax, orientation='vertical')
cbar1.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar1.ax.tick_params(labelsize=font_size)

cbar2_ax = fig.add_axes([0.92, 0.1, 0.02, 0.35])  
cbar2 = fig.colorbar(plt_map6, cax=cbar2_ax, orientation='vertical')
cbar2.set_label('Bias of {0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar2.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_maps_clim_bias_{0}_{1}_RegCM5_{2}_v2.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
