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


def open_obs(ref, varname, domain, dt, time_sel=None):

    path = f'/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/dataset/{ref}/{domain}/preproc/'

    if  ref == 'GPM':
        ref_i = 'GPM_IMERG'
        ref_ii = 'GPM_MERGIR'

        if  varname == 'pr':
            file = f'{varname}_{ref_i}_{domain}_1hr_{dt}.nc'
        else:
            file = f'{varname}_{ref_ii}_{domain}_1hr_{dt}.nc'

    else:
        ref_i = 'ERA5_reanalysis'
        ref_ii = 'ERA5_reanalysis'

        if  varname == 'tp':
            file = f'{varname}_{domain}_{ref_i}_1hr_{dt}.nc'
        else:
            file = f'{varname}_{domain}_{ref_ii}_1hr_{dt}.nc'

    ds = xr.open_dataset(path + file)

    rename = {}
    if 'longitude' in ds.coords: rename['longitude'] = 'lon'
    if 'latitude' in ds.coords:  rename['latitude'] = 'lat'
    if 'valid_time' in ds.coords: rename['valid_time'] = 'time'

    if rename:
        ds = ds.rename(rename)

    if ds.lon.max() > 180:
        ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180)).sortby('lon')

    if time_sel is not None:

        if varname == 'Tb':
            ds = ds.sel(time=time_sel, method='nearest')
        elif varname == 'tp':
            ds = ds.sel(time=time_sel, method='nearest')
        elif varname == 'pr':
            ds = ds.sel(time=time_sel, method='nearest')
        else:
            raise ValueError("Unknown time coordinate for variable")

    return ds[varname], ds['lon'], ds['lat']


def open_moaap_mcs(ref, domain, dt, time_sel):

    path = f'/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/dataset/{ref}/{domain}/output/'
    
    # Filename

    if  ref == 'GPM':
        ref_ = 'GPM'
    else:
        ref_ = 'ERA5_evaluation'

    file = f'{dt}_{ref_}_ObjectMasks__dt-1h_MOAAP-masks.nc'    
    ds = xr.open_dataset(path + file)
    
    print(ds['time'])

    # Select time
    ds = ds.sel(time=time_sel, method='nearest')
    
    return ds['MCS_Tb_Objects'], ds['lon'], ds['lat']


# Import vars
Tb_eurr_era5, lon_eurr_era5, lat_eurr_era5 = open_obs('ERA5', 'Tb', 'EURR-3', '2009100100', '2009-10-02T02:00')
tp_eurr_era5, lon_eurr_era5, lat_eurr_era5 = open_obs('ERA5', 'tp', 'EURR-3', '2009100100', '2009-10-02T02:00')
mcs_eurr_era5, lon_mcs_eurr_era5, lat_mcs_eurr_era5 = open_moaap_mcs('ERA5', 'EURR-3', '20091001', '2009-10-02T02:00')
tp_eurr_era5_masked = np.ma.masked_less(tp_eurr_era5, 1)

Tb_eurr_gpm, lon_eurr_gpm, lat_eurr_gpm = open_obs('GPM', 'Tb', 'EURR-3', '200910', '2009-10-02T02:00')
tp_eurr_gpm, lon_eurr_gpm, lat_eurr_gpm = open_obs('GPM', 'pr', 'EURR-3', '200910', '2009-10-02T02:00')
mcs_eurr_gpm, lon_mcs_eurr_gpm, lat_mcs_eurr_gpm = open_moaap_mcs('GPM', 'EURR-3', '20091001', '2009-10-02T02:00')
tp_eurr_gpm_masked = np.ma.masked_less(tp_eurr_gpm, 1)

print(Tb_eurr_era5.shape)
print(tp_eurr_era5.shape)
print(mcs_eurr_era5.shape)
print()
print(Tb_eurr_gpm.shape)
print(tp_eurr_gpm.shape)
print(mcs_eurr_gpm.shape)

# Plot figure
fig = plt.figure(figsize=(10, 6))
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
		label.set_fontsize(8)

tb_levels = np.arange(220, 321, 1)       # K
tp_levels = np.arange(0, 20.5, 0.5)      # mm/h
mcs_levels = np.arange(0, 1.1, 1)        # MCS objects

# ERA5
plt.subplot(1, 2, 1, projection=ccrs.PlateCarree())
ax = plt.gca()
cf1 = plt.contourf(lon_eurr_era5, lat_eurr_era5, Tb_eurr_era5, levels=tb_levels, cmap='gist_gray_r', extend='both', transform=ccrs.PlateCarree())
cf2 = plt.contourf(lon_eurr_era5, lat_eurr_era5, tp_eurr_era5_masked, levels=tp_levels, cmap='Pastel1', extend='max', transform=ccrs.PlateCarree())
mcs_contour = plt.contour(lon_mcs_eurr_era5, lat_mcs_eurr_era5, mcs_eurr_era5, levels=[0.5], colors='#33a02c', linewidths=2, transform=ccrs.PlateCarree())
plt.title('(a) ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax, lon_eurr_gpm, lat_eurr_gpm)

# GPM
plt.subplot(1, 2, 2, projection=ccrs.PlateCarree())
ax = plt.gca()
cf1 = plt.contourf(lon_eurr_gpm, lat_eurr_gpm, Tb_eurr_gpm, levels=tb_levels, cmap='gist_gray_r', extend='both', transform=ccrs.PlateCarree())
cf2 = plt.contourf(lon_eurr_gpm, lat_eurr_gpm, tp_eurr_gpm_masked, levels=tp_levels, cmap='Pastel1', extend='max', transform=ccrs.PlateCarree())
mcs_contour = plt.contour(lon_mcs_eurr_gpm, lat_mcs_eurr_gpm, mcs_eurr_gpm, levels=[0.5], colors='#33a02c', linewidths=2, transform=ccrs.PlateCarree())
plt.title('(b) GPM (IMERG+MERGIR)', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax, lon_eurr_gpm, lat_eurr_gpm)

# Add MCS legend
plt.plot([], [], color='#33a02c', linewidth=2, label='MCS')
ax.legend(bbox_to_anchor=(0.75, 1.25), loc='upper left', frameon=True, fontsize=8)

cax1 = fig.add_axes([0.30, 0.3, 0.40, 0.025])
cbar1 = fig.colorbar(cf1, cax=cax1, orientation='horizontal')
cbar1.set_label('Tb (K)', fontsize=font_size, fontweight='bold')
cbar1.ax.tick_params(labelsize=font_size)

cax2 = fig.add_axes([0.30, 0.2, 0.40, 0.025])
cbar2 = fig.colorbar(cf2, cax=cax2, orientation='horizontal')
cbar2.set_label('Precipitation (mm/h)', fontsize=font_size, fontweight='bold')
cbar2.ax.tick_params(labelsize=font_size)

# save
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/figs/workshop'
name_out = f'pyplt_maps_moaap_mcs_event_EUR_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
