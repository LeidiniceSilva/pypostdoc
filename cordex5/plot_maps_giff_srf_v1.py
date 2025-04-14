# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "April 14, 2025"
__description__ = "This script plot giff"

import os
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat

from cartopy import config
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5'


def import_rcm(param, dataset):

	arq   = '{0}/postproc/rcm/{1}_CSAM-3_RegCM5_1hr_200001_lonlat.nc'.format(path, param, dataset)    
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean


def configure_subplot(ax):

	ax.set_extent([-85, -30, -60, 15], crs=ccrs.PlateCarree())
	ax.coastlines()
	ax.add_feature(cfeat.LAND, facecolor='beige')
	ax.add_feature(cfeat.OCEAN, facecolor='white')
	ax.add_feature(cfeat.BORDERS, linestyle=':')
	gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray', alpha=0.5)
	gl.xlabel_style = {'color': 'black'}
	gl.ylabel_style = {'color': 'black'}
	gl.xformatter = LongitudeFormatter()
	gl.yformatter = LatitudeFormatter()


# Import model and obs dataset
lat, lon, rcm = import_rcm('pr', 'RegCM5')

# Creating mask of the border
rcm_ = rcm[30:-30, 30:-30]
border_mask = np.full((rcm.shape[1], rcm.shape[2]), np.nan)
border_mask[:29, :] = 1
border_mask[-29:, :] = 1
border_mask[:, :29] = 1
border_mask[:, -29:] = 1
lon_ = lon[30:-30,30:-30]
lat_ = lat[30:-30,30:-30]

t1 = pd.to_datetime('2000-01-01 01:00:00')
t2 = pd.to_datetime('2000-02-01 00:00:00')
series = pd.date_range(t1,t2,freq='60min')
iso8601_t = [t.strftime('%Y/%m/%d %H') for t in series]
iso8601_s = [t.strftime('%Y%m%dT%H') for t in series]

for i in range(0, rcm.shape[0]):
	print(iso8601_t[i])

	# Plot figure
	fig, ax = plt.subplots(1, 1, figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})
	contour = ax.contourf(lon, lat, rcm[i], levels=np.arange(0, 5.25, 0.25), cmap=cm.BuPu, extend='max', transform=ccrs.PlateCarree())
	cbar = fig.colorbar(contour, ax=ax, orientation='vertical', shrink=0.75, pad=0.07)

	ax.contourf(lon, lat, border_mask, levels=[0, 1], cmap='gray')
	configure_subplot(ax)

	# Path out to save figure
	path_out = '{0}/figs/'.format(path)
	name_out = 'pyplt_maps_giff_pr_RegCM5_{0}_ttt.png'.format(iso8601_s[i])
	plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
	plt.close('all')
	plt.cla()

exit()
