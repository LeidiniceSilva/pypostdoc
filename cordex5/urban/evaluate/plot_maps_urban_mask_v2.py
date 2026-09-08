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

city = 'Buenos_Aires' # São_Paulo Buenos_Aires
city_ = 'BA'
freq = 'day'

if freq == 'mon':
	dt = '200001'
	title = 'Jan 2000'
else:
	dt = '20000101'
	title = '01 Jan 2000'

inputs = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/postproc/urb'
mask_ctrl = 'urmask_{0}-CSAM-03_ECMWF-ERA5_evaluation_r1i1p1f1_RegCM5_fx.nc'.format(city)

def import_dataset(fname, var):

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
rcm_ctrl_hfls, lat_ctrl, lon_ctrl = import_dataset('{0}/{1}_{2}_CSAM-3_RegCM5_{3}_2000-2000.nc'.format(inputs, 'hfls', city_, freq), 'hfls')
rcm_ctrl_hfss, lat_ctrl, lon_ctrl = import_dataset('{0}/{1}_{2}_CSAM-3_RegCM5_{3}_2000-2000.nc'.format(inputs, 'hfss', city_, freq), 'hfss')
rcm_ctrl_tas, lat_ctrl, lon_ctrl = import_dataset('{0}/{1}_{2}_CSAM-3_RegCM5_{3}_2000-2000.nc'.format(inputs, 'tas', city_, freq), 'tas')
rcm_ctrl_sfcWind, lat_ctrl, lon_ctrl = import_dataset('{0}/{1}_{2}_CSAM-3_RegCM5_{3}_2000-2000.nc'.format(inputs, 'sfcWind', city_, freq), 'sfcWind')

rcm_urb_hfls, lat_urb, lon_urb = import_dataset('{0}/{1}_{2}_CSAM-3_RegCM5-urb_{3}_2000-2000.nc'.format(inputs, 'hfls', city_, freq), 'hfls')
rcm_urb_hfss, lat_urb, lon_urb = import_dataset('{0}/{1}_{2}_CSAM-3_RegCM5-urb_{3}_2000-2000.nc'.format(inputs, 'hfss', city_, freq), 'hfss')
rcm_urb_tas, lat_urb, lon_urb = import_dataset('{0}/{1}_{2}_CSAM-3_RegCM5-urb_{3}_2000-2000.nc'.format(inputs, 'tas', city_, freq), 'tas')
rcm_urb_sfcWind, lat_urb, lon_urb = import_dataset('{0}/{1}_{2}_CSAM-3_RegCM5-urb_{3}_2000-2000.nc'.format(inputs, 'sfcWind', city_, freq), 'sfcWind')

mask_ctrl = import_urban_mask(os.path.join(inputs, mask_ctrl))

rcm_ctrl_mask_hfls = np.ma.masked_where(mask_ctrl < 1, rcm_ctrl_hfls[0, :, :])
rcm_ctrl_mask_hfss = np.ma.masked_where(mask_ctrl < 1, rcm_ctrl_hfss[0, :, :])
rcm_ctrl_mask_tas = np.ma.masked_where(mask_ctrl < 1, rcm_ctrl_tas[0, :, :])
rcm_ctrl_mask_sfcWind = np.ma.masked_where(mask_ctrl < 1, rcm_ctrl_sfcWind[0, :, :])

rcm_urb_mask_hfls = np.ma.masked_where(mask_ctrl < 1, rcm_urb_hfls[0, :, :])
rcm_urb_mask_hfss = np.ma.masked_where(mask_ctrl < 1, rcm_urb_hfss[0, :, :])
rcm_urb_mask_tas = np.ma.masked_where(mask_ctrl < 1, rcm_urb_tas[0, :, :])
rcm_urb_mask_sfcWind = np.ma.masked_where(mask_ctrl < 1, rcm_urb_sfcWind[0, :, :])

print('Ctrl shape hfls:', rcm_ctrl_mask_hfls.shape)
print('Urb shape hfls:', rcm_urb_mask_hfls.shape)
print('Urban mask shape:', mask_ctrl.shape)

# Plot figure
def format_map(ax, lonlat, dx=0.4):

	ax.set_extent(lonlat)
	ax.set_xticks(np.arange(lonlat[0], lonlat[1], dx), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(lonlat[2], lonlat[3], dx), crs=ccrs.PlateCarree())
	ax.add_feature(cfeature.COASTLINE, edgecolor='black')
	ax.add_feature(cfeature.BORDERS, linestyle=':')
	ax.gridlines(draw_labels=False, linewidth=0.5, color='black', alpha=0.5, linestyle='--')

fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree()})
font_size = 10

if city_ == 'SP':
	lonlat = [-47.2, -46, -24, -23.2]
else:
	lonlat = [-59.3, -57.7, -35.4, -34]

# Urban Fraction Plot
ax1 = axes[0, 0]
im = ax1.contourf(lon_ctrl, lat_ctrl, rcm_ctrl_mask_hfls, cmap='jet', levels=np.arange(0, 100, 1), transform=ccrs.PlateCarree())
ax1.set_title('(a) CTRL hfls {0}'.format(title), loc='left', fontsize=font_size, fontweight='bold')
ax1.text(-38, -38, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(im, ax=ax1, orientation='horizontal', pad=0.12)
cbar.ax.tick_params(labelsize=6)
format_map(ax1, lonlat)

ax2 = axes[0, 1]
im = ax2.contourf(lon_ctrl, lat_ctrl, rcm_ctrl_mask_hfss, cmap='jet', levels=np.arange(0, 150, 1), transform=ccrs.PlateCarree())
ax2.set_title('(b) CTRL hfss {0}'.format(title), loc='left', fontsize=font_size, fontweight='bold')
ax2.text(-38, -38, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(im, ax=ax2, orientation='horizontal', pad=0.12)
cbar.ax.tick_params(labelsize=6)
format_map(ax2, lonlat)

ax3 = axes[0, 2]
im = ax3.contourf(lon_ctrl, lat_ctrl, rcm_ctrl_mask_tas-273.15, cmap='jet', levels=np.arange(10, 30, 0.5), transform=ccrs.PlateCarree())
ax3.set_title('(c) CTRL tas {0}'.format(title), loc='left', fontsize=font_size, fontweight='bold')
ax3.text(-38, -38, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(im, ax=ax3, orientation='horizontal', pad=0.12)
cbar.ax.tick_params(labelsize=6)
format_map(ax3, lonlat)

ax4 = axes[0, 3]
im = ax4.contourf(lon_ctrl, lat_ctrl, rcm_ctrl_mask_sfcWind, cmap='jet', levels=np.arange(0, 6, 0.1), transform=ccrs.PlateCarree())
ax4.set_title('(d) CTRL sfcWind {0}'.format(title), loc='left', fontsize=font_size, fontweight='bold')
ax4.text(-38, -38, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(im, ax=ax4, orientation='horizontal', pad=0.12)
cbar.ax.tick_params(labelsize=6)
format_map(ax4, lonlat)

# Urban Fraction Plot
ax5 = axes[1, 0]
ax5.set_title('(e) URB hfls {0}'.format(title), loc='left', fontsize=font_size, fontweight='bold')
ax5.text(-38, -38, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
im = ax5.contourf(lon_urb, lat_urb, rcm_urb_mask_hfls, cmap='jet', levels=np.arange(0, 100, 1), transform=ccrs.PlateCarree())
cbar = plt.colorbar(im, ax=ax5, orientation='horizontal', pad=0.12)
cbar.ax.tick_params(labelsize=6)
format_map(ax5, lonlat)

ax6 = axes[1, 1]
im = ax6.contourf(lon_urb, lat_urb, rcm_urb_mask_hfss, cmap='jet', levels=np.arange(0, 150, 1), transform=ccrs.PlateCarree())
ax6.set_title('(f) URB hfss {0}'.format(title), loc='left', fontsize=font_size, fontweight='bold')
ax6.text(-38, -38, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(im, ax=ax6, orientation='horizontal', pad=0.12)
cbar.ax.tick_params(labelsize=6)
format_map(ax6, lonlat)

ax7 = axes[1, 2]
im = ax7.contourf(lon_urb, lat_urb, rcm_urb_mask_tas-273.15, cmap='jet', levels=np.arange(10, 30, 0.5), transform=ccrs.PlateCarree())
ax7.set_title('(g) URB tas {0}'.format(title), loc='left', fontsize=font_size, fontweight='bold')
ax7.text(-38, -38, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(im, ax=ax7, orientation='horizontal', pad=0.12)
cbar.ax.tick_params(labelsize=6)
format_map(ax7, lonlat)

ax8 = axes[1, 3]
im = ax8.contourf(lon_urb, lat_urb, rcm_urb_mask_sfcWind, cmap='jet', levels=np.arange(0, 6, 0.1), transform=ccrs.PlateCarree())
ax8.set_title('(h) URB sfcWind {0}'.format(title), loc='left', fontsize=font_size, fontweight='bold')
ax8.text(-38, -38, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(im, ax=ax8, orientation='horizontal', pad=0.12)
cbar.ax.tick_params(labelsize=6)
format_map(ax8, lonlat)

# Save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/figs/urb'
name_out = 'pyplt_maps_{0}_RegCM5_CSAM-3_{1}.png'.format(city_, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400)
plt.show()
exit()


