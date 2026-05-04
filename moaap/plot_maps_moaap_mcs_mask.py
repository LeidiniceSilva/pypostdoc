# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 03, 2026"
__description__ = "This script plot MCSs"

import os
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat

from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


def open_era5(varname, domain, dt, time_sel=None):

    path = f'/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/ERA5/{domain}/preproc/'
    file = f'{varname}_{domain}_ERA5_reanalysis_1hr_{dt}.nc'

    ds = xr.open_dataset(path + file)

    if ds.longitude.max() > 180:
        ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180)).sortby('longitude')

    if time_sel is not None:

        if varname == 'Tb':
            ds = ds.sel(valid_time=time_sel)
        elif varname == 'tp':
            ds = ds.sel(time=time_sel)
        else:
            raise ValueError("Unknown time coordinate for variable")

    return ds[varname], ds['longitude'], ds['latitude']


# Import vars
Tb_car, lon_car, lat_car = open_era5('Tb', 'CAR-4', '2008070100', '2008-07-06T06:00')
tp_car, _, _ = open_era5('tp', 'CAR-4', '2008070100', '2008-07-06T06:00')

Tb_csam, lon_csam, lat_csam = open_era5('Tb', 'CSAM-3', '2008010100', '2008-01-01T00:00')
tp_csam, _, _ = open_era5('tp', 'CSAM-3', '2008010100', '2008-01-01T00:00')

Tb_eurr, lon_eurr, lat_eurr = open_era5('Tb', 'EURR-3', '2009090100', '2009-09-16T06:00')
tp_eurr, _, _ = open_era5('tp', 'EURR-3', '2009090100', '2009-09-16T06:00')

tp_car_masked = np.ma.masked_less(tp_car, 1)
tp_csam_masked = np.ma.masked_less(tp_csam, 1)
tp_eurr_masked = np.ma.masked_less(tp_eurr, 1)

# Plot figure
plt.figure(figsize=(14, 6))
font_size = 10

def configure_subplot(ax, lon, lat):

	ax.set_extent([float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())], crs=ccrs.PlateCarree())
	xticks = np.linspace(float(lon.min()), float(lon.max()), 5)
	yticks = np.linspace(float(lat.min()), float(lat.max()), 5)

	ax.set_xticks(xticks, crs=ccrs.PlateCarree())
	ax.set_yticks(yticks, crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())

	ax.grid(color='gray', ls='--', alpha=0.75)
	ax.coastlines(linewidth=0.5)
	ax.add_feature(cfeat.BORDERS, linewidth=0.5)

	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(6)


tb_levels = np.arange(180, 301, 1)       # K
tp_levels = np.arange(0, 20.5, 0.5)      # mm/h

# CAR-4
plt.subplot(1, 3, 1, projection=ccrs.PlateCarree())
ax = plt.gca()
cf = plt.contourf(lon_car, lat_car, Tb_car, levels=tb_levels, cmap='gist_gray_r', extend='both', transform=ccrs.PlateCarree())
cf2 = plt.contourf(lon_car, lat_car, tp_car_masked, levels=tp_levels, cmap='Pastel1', extend='max', transform=ccrs.PlateCarree())
plt.title('(a)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('CAR-4', fontsize=font_size, fontweight='bold')
configure_subplot(ax, lon_car, lat_car)

# CSAM-3
plt.subplot(1, 3, 2, projection=ccrs.PlateCarree())
ax = plt.gca()
cf1 = plt.contourf(lon_csam, lat_csam, Tb_csam, levels=tb_levels, cmap='gist_gray_r', extend='both', transform=ccrs.PlateCarree())
cf2 = plt.contourf(lon_csam, lat_csam, tp_csam_masked, levels=tp_levels, cmap='Pastel1', extend='max', transform=ccrs.PlateCarree())
plt.title('(b)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('CSAM-3', fontsize=font_size, fontweight='bold')
configure_subplot(ax, lon_csam, lat_csam)

cbar = plt.colorbar(cf1, ax=ax, fraction=0.06, pad=0.08, aspect=40, orientation='horizontal')
cbar.set_label('Tb (K)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

cbar = plt.colorbar(cf2, ax=ax, fraction=0.06, pad=0.1, aspect=40, orientation='horizontal')
cbar.set_label('Precipitation (mm/h)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# EURR-3 
plt.subplot(1, 3, 3, projection=ccrs.PlateCarree())
ax = plt.gca()
cf1 = plt.contourf(lon_eurr, lat_eurr, Tb_eurr, levels=tb_levels, cmap='gist_gray_r', extend='both', transform=ccrs.PlateCarree())
cf2 = plt.contourf(lon_eurr, lat_eurr, tp_eurr_masked, levels=tp_levels, cmap='Pastel1', extend='max', transform=ccrs.PlateCarree())
plt.title('(c)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('EURR-3', fontsize=font_size, fontweight='bold')
configure_subplot(ax, lon_eurr, lat_eurr)

# save
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/figs'
name_out = f'pyplt_maps_moaap_mcs_mask_domains_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()



