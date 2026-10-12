
# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 03, 2026"
__description__ = "This script plot MCSs"

import os
import glob
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from tqdm import tqdm
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


def load_dataset(path_, pattern="*_MOAAP-masks.nc"):

    data_path = f"/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/dataset/{path_}"
    file_list = sorted(glob.glob(os.path.join(data_path, pattern)))

    mcs_sum, lat, lon = None, None, None
    for f in tqdm(file_list):
        ds = xr.open_dataset(f)
        mcs = ds["MCS_Tb_Objects"].values

        if mcs_sum is None:
            mcs_sum = np.zeros((mcs.shape[1], mcs.shape[2]))
            lat = ds["lat"].values
            lon = ds["lon"].values 

        mcs_sum += np.nansum(mcs, axis=0) / 100
        ds.close()

    return mcs_sum, lat, lon


def configure_subplot(ax, lon, lat):

	ax.set_extent([-12, 26, 36, 58], crs=ccrs.PlateCarree())
	xticks = np.linspace(-12, 26, 2)
	yticks = np.linspace(36, 58, 2)

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
mcs_eur_gpm, lat_eur_gpm, lon_eur_gpm = load_dataset('/GPM/EURR-3/output')
mcs_eur_cpm_eval, lat_eur_cpm_eval, lon_eur_cpm_eval = load_dataset('/CPMs/ICTP/EURR-3/evaluation/ERA5/output')
mcs_eur_cpm_hist, lat_eur_cpm_hist, lon_eur_cpm_hist = load_dataset('/CPMs/ICTP/EURR-3/historical/ECEarth/output')
mcs_eur_rcm_eval, lat_eur_rcm_eval, lon_eur_rcm_eval = load_dataset('/RCMs/ICTP/EUR-12/evaluation/ERA5/output')
mcs_eur_rcm_hist, lat_eur_rcm_hist, lon_eur_rcm_hist = load_dataset('/RCMs/ICTP/EUR-12/historical/ECEarth/output')

print(mcs_eur_gpm)
print()
print(mcs_eur_cpm_eval)

# Plot figure
fig = plt.figure(figsize=(12, 5))
font_size = 10

orig_cmap = plt.cm.get_cmap('RdYlBu_r', 256)
colors = orig_cmap(np.linspace(0, 1, 256))
colors[0] = [1, 1, 1, 1]
cmap = mcolors.ListedColormap(colors)
mcs_levels = np.arange(1, 51, 1)

# GPM
ax1 = fig.add_subplot(2, 3, 1, projection=ccrs.PlateCarree())
cf = ax1.contourf(lon_eur_gpm, lat_eur_gpm, mcs_eur_gpm/15, levels=mcs_levels, cmap=cmap, extend='min', transform=ccrs.PlateCarree())
plt.title('(a) GPM', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax1, lon_eur_gpm, lat_eur_gpm)

# CPM eval
ax2 = fig.add_subplot(2, 3, 2, projection=ccrs.PlateCarree())
cf = ax2.contourf(lon_eur_cpm_eval, lat_eur_cpm_eval, mcs_eur_cpm_eval, levels=mcs_levels, cmap=cmap, extend='min', transform=ccrs.PlateCarree())
plt.title('(b) EURR-3 Eval', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2, lon_eur_gpm, lat_eur_gpm)

# CPM hist
ax3 = fig.add_subplot(2, 3, 3, projection=ccrs.PlateCarree())
cf = ax3.contourf(lon_eur_cpm_hist, lat_eur_cpm_hist, mcs_eur_cpm_hist, levels=mcs_levels, cmap=cmap, extend='min', transform=ccrs.PlateCarree())
plt.title('(c) EURR-3 Hist', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3, lon_eur_gpm, lat_eur_gpm)

# RCM eval
ax4 = fig.add_subplot(2, 3, 5, projection=ccrs.PlateCarree())
cf = ax4.contourf(lon_eur_rcm_eval, lat_eur_rcm_eval, mcs_eur_rcm_eval, levels=mcs_levels, cmap=cmap, extend='min', transform=ccrs.PlateCarree())
plt.title('(d) EUR-12 Eval', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax4, lon_eur_gpm, lat_eur_gpm)

# RCM hist
ax5 = fig.add_subplot(2, 3, 6, projection=ccrs.PlateCarree())
cf = ax5.contourf(lon_eur_rcm_hist, lat_eur_rcm_hist, mcs_eur_rcm_hist, levels=mcs_levels, cmap=cmap, extend='min', transform=ccrs.PlateCarree())
plt.title('(e) EUR-12 Hist', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax5, lon_eur_gpm, lat_eur_gpm)

cbar_ax = fig.add_axes([0.25, 0.05, 0.5, 0.02])  # [left, bottom, width, height]
cbar = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal')
cbar.set_label('MCS frequency (days/year)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/figs/v2'
name_out = f'pyplt_maps_moaap_mcs_clim_{domain}_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


