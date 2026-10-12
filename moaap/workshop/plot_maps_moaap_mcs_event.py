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


def open_obs(varname, domain, dt, time_sel=None):

    path = f'/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/dataset/GPM/{domain}/preproc/'

    if  varname == 'pr':
        file = f'{varname}_GPM_IMERG_{domain}_1hr_{dt}.nc'
    else:
        file = f'{varname}_GPM_MERGIR_{domain}_1hr_{dt}.nc'

    ds = xr.open_dataset(path + file)

    if ds.lon.max() > 180:
        ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180)).sortby('lon')

    if time_sel is not None:

        if varname == 'Tb':
            ds = ds.sel(time=time_sel, method='nearest')
        elif varname == 'pr':
            ds = ds.sel(time=time_sel, method='nearest')
        else:
            raise ValueError("Unknown time coordinate for variable")

    return ds[varname], ds['lon'], ds['lat']


def open_moaap_mcs(domain, dt, time_sel):

    path = f'/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/dataset/GPM/{domain}/output/'
    
    # Filename
    if domain == 'CAR-4':
        file = f'{dt}_GPM_ObjectMasks__dt-1h_MOAAP-masks.nc'
    elif domain == 'CSAM-3':
        file = f'{dt}_GPM_ObjectMasks__dt-1h_MOAAP-masks.nc'
    elif domain == 'EURR-3':
        file = f'{dt}_GPM_ObjectMasks__dt-1h_MOAAP-masks.nc'
    else:
        raise ValueError("Unknown domain")
    
    ds = xr.open_dataset(path + file)
    
    print(ds['time'])

    # Select time
    ds = ds.sel(time=time_sel, method='nearest')
    
    return ds['MCS_Tb_Objects'], ds['lon'], ds['lat']


# Import vars
Tb_car, lon_car, lat_car = open_obs('Tb', 'CAR-4', '200807', '2008-07-06T06:00')
tp_car, _, _ = open_obs('pr', 'CAR-4', '200807', '2008-07-06T06:00')
mcs_car, lon_mcs_car, lat_mcs_car = open_moaap_mcs('CAR-4', '20080701', '2008-07-06T06:00')

Tb_csam, lon_csam, lat_csam = open_obs('Tb', 'CSAM-3', '200712', '2007-12-05T18:00')
tp_csam, _, _ = open_obs('pr', 'CSAM-3', '200712', '2007-12-05T18:00')
mcs_csam, lon_mcs_csam, lat_mcs_csam = open_moaap_mcs('CSAM-3', '20071201', '2007-12-05T18:00')

Tb_eurr, lon_eurr, lat_eurr = open_obs('Tb', 'EURR-3', '200909', '2009-09-16T06:00')
tp_eurr, _, _ = open_obs('pr', 'EURR-3', '200909', '2009-09-16T06:00')
mcs_eurr, lon_mcs_eurr, lat_mcs_eurr = open_moaap_mcs('EURR-3', '20090901', '2009-09-16T06:00')

tp_car_masked = np.ma.masked_less(tp_car, 1)
tp_csam_masked = np.ma.masked_less(tp_csam, 1)
tp_eurr_masked = np.ma.masked_less(tp_eurr, 1)

print(Tb_car.shape)
print(tp_car.shape)
print(mcs_car.shape)

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


tb_levels = np.arange(220, 321, 1)       # K
tp_levels = np.arange(0, 20.5, 0.5)      # mm/h
mcs_levels = np.arange(0, 1.1, 1)        # MCS objects

# CAR-4
plt.subplot(1, 3, 1, projection=ccrs.PlateCarree())
ax = plt.gca()
cf = plt.contourf(lon_car, lat_car, Tb_car, levels=tb_levels, cmap='gist_gray_r', extend='both', transform=ccrs.PlateCarree())
cf2 = plt.contourf(lon_car, lat_car, tp_car_masked, levels=tp_levels, cmap='Pastel1', extend='max', transform=ccrs.PlateCarree())
mcs_contour = plt.contour(lon_mcs_car, lat_mcs_car, mcs_car, levels=[0.5], colors='#33a02c', linewidths=2, transform=ccrs.PlateCarree())
plt.title('(a)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('CAR-4', fontsize=font_size, fontweight='bold')
configure_subplot(ax, lon_car, lat_car)

# CSAM-3
plt.subplot(1, 3, 2, projection=ccrs.PlateCarree())
ax = plt.gca()
cf1 = plt.contourf(lon_csam, lat_csam, Tb_csam, levels=tb_levels, cmap='gist_gray_r', extend='both', transform=ccrs.PlateCarree())
cf2 = plt.contourf(lon_csam, lat_csam, tp_csam_masked, levels=tp_levels, cmap='Pastel1', extend='max', transform=ccrs.PlateCarree())
mcs_contour = plt.contour(lon_mcs_csam, lat_mcs_csam, mcs_csam, levels=[0.5], colors='#33a02c', linewidths=2, transform=ccrs.PlateCarree())
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
mcs_contour = plt.contour(lon_mcs_eurr, lat_mcs_eurr, mcs_eurr, levels=[0.5], colors='#33a02c', linewidths=2, transform=ccrs.PlateCarree())
plt.title('(c)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('EURR-3', fontsize=font_size, fontweight='bold')
configure_subplot(ax, lon_eurr, lat_eurr)

# Add MCS legend
plt.plot([], [], color='#33a02c', linewidth=2, label='MCS')
ax.legend(bbox_to_anchor=(0.75, 1.25), loc='upper left', frameon=True, fontsize=8)

# save
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/figs/v2'
name_out = f'pyplt_maps_moaap_mcs_event_domains_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
