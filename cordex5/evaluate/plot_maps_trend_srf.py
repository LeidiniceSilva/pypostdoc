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
from scipy.stats import linregress

from cartopy import config
from matplotlib.colors import LinearSegmentedColormap
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

var = 'pr'
domain = 'CSAM-3'
idt, fdt = '2000', '2009'
dt = '{0}-{1}'.format(idt, fdt)
font_size = 8

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5'

def import_obs(param, domain, dataset):

	arq   = '{0}/postproc/evaluate/obs/{1}_{2}_{3}.nc'.format(path, param, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean


def import_rcm(param, domain, dataset):

	arq   = '{0}/postproc/evaluate/obs/{1}_{2}_{3}.nc'.format(path, param, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean


def comp_trend_sig(data, min_valid=10, per_decade=False):
	
	nt, ny, nx = data.shape
	t = np.arange(nt)
	trend = np.full((ny, nx), np.nan)
	pval  = np.full((ny, nx), np.nan)

	for i in range(ny):
		for j in range(nx):
			y = data[:, i, j]
			mask = np.isfinite(y)

			if mask.sum() >= min_valid:
				res = linregress(t[mask], y[mask])
				trend[i, j] = res.slope
				pval[i, j]  = res.pvalue

	if per_decade:
		trend = trend * 10.0

	return trend, pval


# Import model and obs dataset
dict_var = {'pr': ['pre', 'tp'],
'tas': ['tmp', 't2m'],
'tasmax': ['tmx', 'tasmax'],
'tasmin': ['tmn', 'tasmin']}

lat, lon, cru_yr = import_obs(dict_var[var][0], domain, 'CRU')
#lat, lon, era5_yr = import_obs(dict_var[var][1], domain, 'ERA5')
#lat, lon, regcm_yr = import_rcm(var, domain, 'RegCM5')

cru_trend, cru_pval = comp_trend_sig(cru_yr, min_valid=10, per_decade=False)
cru_sig = np.ma.masked_where(cru_pval >= 0.05, cru_trend)

# Plot figure   
def configure_subplot(ax):

	lon_min = np.round(np.min(-78.81965), 1)
	lon_max = np.round(np.max(-35.32753), 1)
	lat_min = np.round(np.min(-36.70233), 1)
	lat_max = np.round(np.max(-12.24439), 1)
	ax.set_extent([np.min(-78.81965), np.max(-35.32753), np.min(-36.70233), np.max(-12.24439)], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(lon_min,lon_max,10), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(lat_min,lat_max,5), crs=ccrs.PlateCarree())

	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(6)

	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.add_feature(cfeat.BORDERS, linewidth=0.5)
	ax.coastlines(linewidth=0.5)


dict_plot = {'pr': ['Precipitation trend (mm yr$^-$$^1$)', np.arange(-20, 20.5, 0.5), cm.coolwarm],
'tas': ['Air temperature trend (°C yr$^-$$^1$)', np.arange(-10, 10.5, 0.5), cm.bwr],
'tasmax': ['Maximum air temperature trend (°C yr$^-$$^1$)', np.arange(-10, 10.5, 0.5), cm.bwr],
'tasmin': ['Minimum air temperature trend (°C yr$^-$$^1$)', np.arange(-10, 10.5, 0.5), cm.bwr]}

fig, axes = plt.subplots(1, 3, figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree()})

ax1 = axes[0]
plt_map = ax1.contourf(lon, lat, cru_trend, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
ax1.contourf(lon, lat, cru_sig, hatches=['...'], colors='none', transform=ccrs.PlateCarree())
ax1.set_title(u'(a) CRU', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

ax2 = axes[1]
plt_map = ax2.contourf(lon, lat, cru_trend, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
ax2.contourf(lon, lat, cru_sig, hatches=['...'], colors='none', transform=ccrs.PlateCarree())
ax2.set_title(u'(b) ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

ax3 = axes[2]
plt_map = ax3.contourf(lon, lat, cru_trend, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
ax3.contourf(lon, lat, cru_sig, hatches=['...'], colors='none', transform=ccrs.PlateCarree())
ax3.set_title(u'(c) RegCM5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

# Set colobar
cbar = fig.colorbar(plt_map, ax=fig.axes, orientation='horizontal', pad=0.1, aspect=50)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/figs/evaluate/rcm'.format(path)
name_out = 'pyplt_maps_trend_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


