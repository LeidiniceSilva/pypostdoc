# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 03, 2026"
__description__ = "This script plot MCSs"

import os
import glob
import numpy as np
import xarray as xr
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat

from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


def open_mcs_fraction(domain):

    path = f'/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/ERA5/{domain}/output/'
    pattern = '*_ERA5_evaluation_ObjectMasks__dt-1h_MOAAP-masks.nc'

    files = sorted(glob.glob(os.path.join(path, pattern)))

    pr_tot = None
    pr_mcs = None

    for f in files:
        print(f)

        ds = xr.open_dataset(f)

        pr = ds['PR']
        mcs_mask = (ds['MCS_Tb_Objects'] > 0)

        pr_mcs_tmp = pr * mcs_mask

        if pr_tot is None:
            pr_tot = pr.sum(dim='time')
            pr_mcs = pr_mcs_tmp.sum(dim='time')
            lon = ds['lon']
            lat = ds['lat']
        else:
            pr_tot += pr.sum(dim='time')
            pr_mcs += pr_mcs_tmp.sum(dim='time')

    frac = (pr_mcs / pr_tot.where(pr_tot > 0)) * 100.

    return frac, lon, lat


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

# Import vars
tp_car, lon_car, lat_car = open_mcs_fraction('CAR-4')
lon_car = ((lon_car + 180) % 360) - 180

tp_csam, lon_csam, lat_csam = open_mcs_fraction('CSAM-3')
lon_csam = ((lon_csam + 180) % 360) - 180

tp_eurr, lon_eurr, lat_eurr = open_mcs_fraction('EURR-3')

# Plot figure
fig = plt.figure(figsize=(6, 12))
font_size = 10

color = ['#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff',
        '#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']
cmap = matplotlib.colors.ListedColormap(color)
cmap.set_under('white') 
tp_levels = np.arange(1,100,5)

# CAR-4
ax1 = fig.add_subplot(3, 1, 1, projection=ccrs.PlateCarree())
cf = ax1.contourf(lon_car, lat_car, tp_car, levels=tp_levels, cmap=cmap, extend='both', transform=ccrs.PlateCarree())
plt.title('(a)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('CAR-4', fontsize=font_size, fontweight='bold')
configure_subplot(ax1, lon_car, lat_car)

# CSAM-3
ax2 = fig.add_subplot(3, 1, 2, projection=ccrs.PlateCarree())
cf = ax2.contourf(lon_csam, lat_csam, tp_csam, levels=tp_levels, cmap=cmap, extend='both', transform=ccrs.PlateCarree())
plt.title('(b)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('CSAM-3', fontsize=font_size, fontweight='bold')
configure_subplot(ax2, lon_csam, lat_csam)

# EURR-3 
ax3 = fig.add_subplot(3, 1, 3, projection=ccrs.PlateCarree())
cf = ax3.contourf(lon_eurr, lat_eurr, tp_eurr, levels=tp_levels, cmap=cmap, extend='both', transform=ccrs.PlateCarree())
plt.title('(c)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('EURR-3', fontsize=font_size, fontweight='bold')
configure_subplot(ax3, lon_eurr, lat_eurr)

cbar = fig.colorbar(cf, ax=[ax1, ax2, ax3], orientation='horizontal', fraction=0.05, pad=0.04, aspect=40)
cbar.set_label('Precipitation fraction (%)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/figs'
name_out = f'pyplt_maps_moaap_mcs_frac_domains_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()





