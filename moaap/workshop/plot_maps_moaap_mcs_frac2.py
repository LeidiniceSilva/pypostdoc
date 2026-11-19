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


def load_dataset(path_):

    path = f'/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/dataset/{path_}'
    pattern = '*_MOAAP-masks.nc'

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

    return frac, lat, lon


def configure_subplot(ax, lon, lat):

	ax.set_extent([-10, 20, 38, 53], crs=ccrs.PlateCarree())
	xticks = np.linspace(-10, 20, 4)
	yticks = np.linspace(38, 53, 4)

	ax.set_xticks(xticks, crs=ccrs.PlateCarree())
	ax.set_yticks(yticks, crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())

	ax.grid(color='gray', ls='--', alpha=0.75)
	ax.coastlines(linewidth=0.5)
	ax.add_feature(cfeat.BORDERS, linewidth=0.5)

	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(8)


# domain
domain = 'EUR'

# Import datasets
mcs_era5, lat_era5, lon_era5 = load_dataset('/ERA5/EURR-3/output')
mcs_gpm,  lat_gpm,  lon_gpm  = load_dataset('/GPM/EURR-3/output')

mcs_cnrm_cpm_eval, lat_cnrm_cpm_eval, lon_cnrm_cpm_eval = load_dataset('/CPMs/CNRM/ALPX-3/evaluation/ERA5/output')
mcs_cnrm_cpm_hist, lat_cnrm_cpm_hist, lon_cnrm_cpm_hist = load_dataset('/CPMs/CNRM/ALPX-3/historical/CNRM-ESM2-1/output')
mcs_cnrm_rcm_eval, lat_cnrm_rcm_eval, lon_cnrm_rcm_eval = load_dataset('/RCMs/CNRM/EUR-12/evaluation/ERA5/output')
mcs_cnrm_rcm_hist, lat_cnrm_rcm_hist, lon_cnrm_rcm_hist = load_dataset('/RCMs/CNRM/EUR-12/historical/CNRM-ESM2-1/output')

mcs_ictp_rcm_eval, lat_ictp_rcm_eval, lon_ictp_rcm_eval = load_dataset('/RCMs/ICTP/EUR-12/evaluation/ERA5/output')
mcs_ictp_rcm_hist, lat_ictp_rcm_hist, lon_ictp_rcm_hist = load_dataset('/RCMs/ICTP/EUR-12/historical/ECEarth/output')

print(mcs_era5)
print()
print(mcs_gpm)
print()
print(mcs_cnrm_cpm_eval)
print()
print(mcs_cnrm_cpm_hist)
print()
print(mcs_cnrm_rcm_eval)
print()
print(mcs_cnrm_rcm_eval)
print()
print(mcs_ictp_rcm_eval)
print()
print(mcs_ictp_rcm_eval)

# Plot figure
fig = plt.figure(figsize=(10, 12))
font_size = 10

color = ['#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff',
        '#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']
cmap = matplotlib.colors.ListedColormap(color)
cmap.set_under('white') 
mcs_levels = np.arange(1,51,1)

# ERA5
ax1 = fig.add_subplot(2, 4, 1, projection=ccrs.PlateCarree())
cf = ax1.contourf(lon_era5, lat_era5, mcs_era5, levels=mcs_levels, cmap=cmap, extend='min', transform=ccrs.PlateCarree())
plt.title('(a) ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax1, lon_cnrm_cpm_eval, lat_cnrm_cpm_eval)

# GPM
ax2 = fig.add_subplot(2, 4, 2, projection=ccrs.PlateCarree())
cf = ax2.contourf(lon_gpm, lat_gpm, mcs_gpm, levels=mcs_levels, cmap=cmap, extend='min', transform=ccrs.PlateCarree())
plt.title('(b) GPM', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2, lon_cnrm_cpm_eval, lat_cnrm_cpm_eval)

# CNRM ALPX-3 eval
ax3 = fig.add_subplot(2, 4, 3, projection=ccrs.PlateCarree())
cf = ax3.contourf(lon_cnrm_cpm_eval, lat_cnrm_cpm_eval, mcs_cnrm_cpm_eval, levels=mcs_levels, cmap=cmap, extend='min', transform=ccrs.PlateCarree())
plt.title('(c) CNRM ALPX-3 Eval', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3, lon_cnrm_cpm_eval, lat_cnrm_cpm_eval)

# CNRM ALPX-3 hist
ax4 = fig.add_subplot(2, 4, 4, projection=ccrs.PlateCarree())
cf = ax4.contourf(lon_cnrm_cpm_hist, lat_cnrm_cpm_hist, mcs_cnrm_cpm_hist, levels=mcs_levels, cmap=cmap, extend='min', transform=ccrs.PlateCarree())
plt.title('(d) CNRM ALPX-3 Hist', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax4, lon_cnrm_cpm_eval, lat_cnrm_cpm_eval)

# CNRM EUR-12 eval
ax5 = fig.add_subplot(2, 4, 5, projection=ccrs.PlateCarree())
cf = ax5.contourf(lon_cnrm_rcm_eval, lat_cnrm_rcm_eval, mcs_cnrm_rcm_eval, levels=mcs_levels, cmap=cmap, extend='min', transform=ccrs.PlateCarree())
plt.title('(e) CNRM EUR-12 Eval', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax5, lon_cnrm_cpm_eval, lat_cnrm_cpm_eval)

# CNRM EUR-12 hist
ax6 = fig.add_subplot(2, 4, 6, projection=ccrs.PlateCarree())
cf = ax6.contourf(lon_cnrm_rcm_hist, lat_cnrm_rcm_hist, mcs_cnrm_rcm_hist, levels=mcs_levels, cmap=cmap, extend='min', transform=ccrs.PlateCarree())
plt.title('(f) CNRM EUR-12 Hist', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax6, lon_cnrm_cpm_eval, lat_cnrm_cpm_eval)

# ICTP EUR-12 eval
ax7 = fig.add_subplot(2, 4, 7, projection=ccrs.PlateCarree())
cf = ax7.contourf(lon_ictp_rcm_eval, lat_ictp_rcm_eval, mcs_ictp_rcm_eval, levels=mcs_levels, cmap=cmap, extend='min', transform=ccrs.PlateCarree())
plt.title('(e) ICTP EUR-12 Eval', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax7, lon_cnrm_cpm_eval, lat_cnrm_cpm_eval)

# ICTP EUR-12 hist
ax8 = fig.add_subplot(2, 4, 8, projection=ccrs.PlateCarree())
cf = ax8.contourf(lon_ictp_rcm_hist, lat_ictp_rcm_hist, mcs_ictp_rcm_hist, levels=mcs_levels, cmap=cmap, extend='min', transform=ccrs.PlateCarree())
plt.title('(f) ICTP EUR-12 Hist', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax8, lon_cnrm_cpm_eval, lat_cnrm_cpm_eval)

cbar_ax = fig.add_axes([0.25, 0.05, 0.5, 0.02])  
cbar = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal')
cbar.set_label('Precipitation fraction (%)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/figs/workshop'
name_out = f'pyplt_maps_moaap_mcs_frac_{domain}_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()





