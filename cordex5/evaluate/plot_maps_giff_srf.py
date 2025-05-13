# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "April 14, 2025"
__description__ = "This script plot giff"

import os
import sys
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat

from cartopy import config
from netCDF4 import Dataset as nc
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5'


def import_dataset(param, domain, dataset):

	arq   = '{0}/postproc/obs/{1}_{2}_{3}_1hr_20000105-20000108.nc'.format(path, param, domain, dataset)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['latitude'][:]
	lon   = data.variables['longitude'][:]
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


# RegCM file
if len(sys.argv) > 1:
    RCMf = nc(sys.argv[1], mode='r')
else:
    RCMf = nc(os.path.join('{0}/postproc/rcm'.format(path),'pr_CSAM-3_RegCM5_1hr_20000105-20000108.nc'), mode='r')
    
lat  = RCMf.variables['xlat'][:,:]
lon  = RCMf.variables['xlon'][:,:]
rcm = RCMf.variables['pr'][:,:]
lonc = RCMf.longitude_of_projection_origin
latc = RCMf.latitude_of_projection_origin
RCMf.close()

# Import model and obs dataset
lat_obs, lon_obs, obs_new = import_dataset('tp', 'CSAM-3', 'ERA5')

print(obs_new.shape)

# Creating mask of the border
border_mask = np.full((rcm.shape[1], rcm.shape[2]), np.nan)
border_mask[:29, :] = 1
border_mask[-29:, :] = 1
border_mask[:, :29] = 1
border_mask[:, -29:] = 1

# Creating datetime
t1 = pd.to_datetime('2000-01-05 01:00:00')
t2 = pd.to_datetime('2000-01-08 00:00:00')
series = pd.date_range(t1,t2,freq='60min')
iso8601_t = [t.strftime('%Y/%m/%d %H') for t in series]
iso8601_s = [t.strftime('%Y%m%dT%H') for t in series]

for i in range(0, rcm.shape[0]):
	print(iso8601_t[i])

	# Plot figure
	fig, ax = plt.subplots(1, 1, figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})
	contour = ax.contourf(lon_obs, lat_obs, obs_new[i], levels=np.arange(0, 5.25, 0.25), cmap=cm.BuPu, extend='max', transform=ccrs.PlateCarree())
	contour = ax.contourf(lon, lat, rcm[i], levels=np.arange(0, 6.5, 0.5), cmap=cm.BuPu, extend='max', transform=ccrs.PlateCarree())
	cbar = fig.colorbar(contour, ax=ax, orientation='vertical', shrink=0.75, pad=0.07)

	ax.contourf(lon, lat, border_mask, levels=[0, 1], cmap='gray')

	ax.set_title(u'Total precipitation {0} UTC'.format(iso8601_t[i]))
	configure_subplot(ax)

	# Path out to save figure
	path_out = '{0}/figs/giff'.format(path)
	name_out = 'pyplt_maps_pr_RegCM5_ERA5_{0}.png'.format(iso8601_s[i])
	plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')
	plt.close('all')
	plt.cla()

exit()

