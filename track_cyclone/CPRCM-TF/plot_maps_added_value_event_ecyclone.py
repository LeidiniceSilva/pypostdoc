# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Sep 22, 2025"
__description__ = "This script plots composite around cyclone center"

import os
import netCDF4
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeat

from matplotlib.patches import Circle
from scipy.interpolate import griddata
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/TRACK-CYCLONE/CPRCM-TF'
font_size = 10
day=24

def import_data(param):
	
	arq = '{0}/data/{1}_Mar{2}_lonlat.nc'.format(path, param, day)
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:]    
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	mean = var[:, :, :]
	
	if param == 'psl':
		#mean = var[-1, :, :] / 100.0
		mean = np.nanmean(var[:, :, :], axis=0) / 100.0
	elif param == 'pr':
		#mean = var[-1, :, :] * 3600
		mean = np.nansum(var[:, :, :], axis=0) * 3600
	else:
		#mean = var[-1, :, :]
		mean = np.nanmean(var[:, :, :], axis=0)
						
	return lat, lon, mean


# Import data
lat, lon, pr = import_data('pr')
lat, lon, psl = import_data('psl')
lat, lon, ta925 = import_data('ta925')
lat, lon, ua925 = import_data('ua925')
lat, lon, va925 = import_data('va925')

def configure_subplot(ax):

	ax.set_extent([-66, -38, -36, -14], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(-66,-38,4), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(-36,-14,4), crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.set_xlabel('Longitude', fontsize=font_size, fontweight='bold')
	ax.set_ylabel('Latitude',  fontsize=font_size, fontweight='bold')
	ax.grid(c='k', ls='--', alpha=0.5)
	ax.add_feature(cfeat.BORDERS, linewidth=2, color='gray')
	ax.coastlines(linewidth=2, color='gray')	


# Plot figure
fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()})
step = 4
levels_pr = np.arange(0, 31.5, 0.5)
levels_psl = np.arange(990, 1020.5, 0.5)
levels_ta925 = np.arange(286, 301, 1)

cf = ax.contourf(lon, lat, pr, levels=levels_pr, cmap='Blues', extend='max')
cb = plt.colorbar(cf)
cb.set_label('Precipitation (mm)', fontsize=font_size, fontweight='bold')

q = ax.quiver(lon[::step], lat[::step], ua925[::step, ::step], va925[::step, ::step], color='green', scale=200, width=0.002) 
ax.quiverkey(q, X=0.8, Y=1.03, U=10.5, label='10 m s⁻¹', labelpos='E', coordinates='axes', fontproperties={'size': 10})
    
ct1 = ax.contour(lon, lat, psl, levels=levels_psl, colors='black', linewidths=0.75)
ax.clabel(ct1, inline=1, fontsize=font_size-3)
configure_subplot(ax)

ct2 = ax.contour(lon, lat, ta925, levels=levels_ta925, colors='red', linewidths=0.75, linestyles='--')
ax.clabel(ct2, inline=1, fontsize=font_size-3)
configure_subplot(ax)

# Save figura
path_out = '{0}/figs/'.format(path)
name_out = 'pyplt_maps_Hurricane_Catarina_{0}Mar2004.png'.format(day)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


