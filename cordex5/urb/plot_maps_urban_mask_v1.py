# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 16, 2023"
__description__ = "This script plot urban mask"

import os
import numpy as np
import netCDF4 as nc
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker

from matplotlib.colors import ListedColormap

var ='hfls' # sfcWind hfss hfls 
city = 'Buenos_Aires' # São_Paulo Buenos_Aires
city_ = 'BA'
freq = 'mon'

if freq == 'mon':
	dt = '200001'
	title = 'Jan 2000'
else:
	dt = '20000101'
	title = '01 Jan 2000'

inputs = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/postproc/urb'
mask_ctrl = 'urmask_{0}-CSAM-03_ECMWF-ERA5_evaluation_r1i1p1f1_RegCM5_fx.nc'.format(city)

def import_dataset(fname):

	with nc.Dataset(fname, 'r') as f:
		var_data = f.variables[var][:]
		lat = f.variables['lat'][:]
		lon = f.variables['lon'][:]

	return var_data, lat, lon


def import_urban_mask(fname):

	with nc.Dataset(fname, 'r') as f:
		urmask = f.variables['urmask'][:]  

	return urmask


# Import dataset
rcm_ctrl, lat_ctrl, lon_ctrl = import_dataset('{0}/{1}_{2}_CSAM-3_RegCM5_{3}_2000-2000.nc'.format(inputs, var, city_, freq))
rcm_urb, lat_urb, lon_urb = import_dataset('{0}/{1}_{2}_CSAM-3_RegCM5-urb_{3}_2000-2000.nc'.format(inputs, var, city_, freq))

mask_ctrl = import_urban_mask(os.path.join(inputs, mask_ctrl))
rcm_ctrl_mask = np.ma.masked_where(mask_ctrl < 1, rcm_ctrl[0, :, :])
rcm_urb_mask = np.ma.masked_where(mask_ctrl < 1, rcm_urb[0, :, :])

print('Ctrl shape:', rcm_ctrl.shape)
print('Urb shape:', rcm_urb.shape)
print('Urban mask shape:', mask_ctrl.shape)

# Plot figure
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})
font_size = 10

if city_ == 'SP':
	lonlat = [-47.2, -46, -24, -23.2]
else:
	lonlat = [-59.3, -57.7, -35.4, -34]

if var == 'sfcWind':
	vimn = 0
	vimx = 5.05
	int = 0.05
elif var == 'hfss':
	vimn = 0
	vimx = 150
	int = 1.5
else:
	vimn = 0
	vimx = 101
	int = 1

# Urban Fraction Plot
ax = axes[0]
ax.set_extent([lonlat[0], lonlat[1], lonlat[2], lonlat[3]])
ax.set_xticks(np.arange(lonlat[0],lonlat[1],0.2), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(lonlat[2],lonlat[3],0.2), crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE, edgecolor='black')
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.set_title('(a) RegCM5 {0} {1}'.format(var, title), loc='left', fontsize=font_size, fontweight='bold')
ax.text(-38, -38, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
im = ax.contourf(lon_ctrl, lat_ctrl, rcm_ctrl_mask, cmap='jet', levels=np.arange(vimn, vimx, int), transform=ccrs.PlateCarree())
cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05)
ax.gridlines(draw_labels=False, linewidth=0.5, color='black', alpha=0.5, linestyle='--')

# Urban Fraction Plot
ax = axes[1]
ax.set_extent([lonlat[0], lonlat[1], lonlat[2], lonlat[3]])
ax.set_xticks(np.arange(lonlat[0],lonlat[1],0.2), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(lonlat[2],lonlat[3],0.2), crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE, edgecolor='black')
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.set_title('(b) RegCM5 URB {0} {1}'.format(var, title), loc='left', fontsize=font_size, fontweight='bold')
ax.text(-38, -38, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
im = ax.contourf(lon_urb, lat_urb, rcm_urb_mask, cmap='jet', levels=np.arange(vimn, vimx, int), transform=ccrs.PlateCarree())
cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05)
ax.gridlines(draw_labels=False, linewidth=0.5, color='black', alpha=0.5, linestyle='--')

# Save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/figs/urb'
name_out = 'pyplt_maps_{0}_{1}_RegCM5_CSAM-3_{2}.png'.format(var, city_, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400)
plt.show()
exit()


