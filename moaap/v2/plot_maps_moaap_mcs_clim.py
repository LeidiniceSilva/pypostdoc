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

def load_obs(domain, pattern="*_MOAAP-masks.nc"):

    data_path = f"/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/dataset/GPM/{domain}/output/"
    file_list = sorted(glob.glob(os.path.join(data_path, pattern)))

    mcs_sum, lat, lon = None, None, None
    for f in tqdm(file_list):
        ds = xr.open_dataset(f)
        mcs = ds["MCS_Tb_Objects"].values  # (time, lat, lon)

        if mcs_sum is None:
            mcs_sum = np.zeros((mcs.shape[1], mcs.shape[2]))
            lat = ds["lat"].values
            lon = ds["lon"].values 

        # climatological occurrence
        mcs_sum += np.nansum(mcs, axis=0) / 100
        ds.close()

    return mcs_sum, lat, lon


def load_cpm(domain, exp, dataset, pattern="*_MOAAP-masks.nc"):

    data_path = f"/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/dataset/CPMs/ICTP/{domain}/{exp}/{dataset}/output/"
    file_list = sorted(glob.glob(os.path.join(data_path, pattern)))

    mcs_sum, lat, lon = None, None, None
    for f in tqdm(file_list):
        ds = xr.open_dataset(f)
        mcs = ds["MCS_Tb_Objects"].values  # (time, lat, lon)

        if mcs_sum is None:
            mcs_sum = np.zeros((mcs.shape[1], mcs.shape[2]))
            lat = ds["lat"].values
            lon = ds["lon"].values 

        # climatological occurrence
        mcs_sum += np.nansum(mcs, axis=0) / 100
        ds.close()

    return mcs_sum, lat, lon


def load_rcm(domain, exp, dataset, pattern="*_MOAAP-masks.nc"):

    data_path = f"/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/dataset/CRMs/ICTP/{domain}/{exp}/{dataset}/output/"
    file_list = sorted(glob.glob(os.path.join(data_path, pattern)))

    mcs_sum, lat, lon = None, None, None
    for f in tqdm(file_list):
        ds = xr.open_dataset(f)
        mcs = ds["MCS_Tb_Objects"].values  # (time, lat, lon)

        if mcs_sum is None:
            mcs_sum = np.zeros((mcs.shape[1], mcs.shape[2]))
            lat = ds["lat"].values
            lon = ds["lon"].values 

        # climatological occurrence
        mcs_sum += np.nansum(mcs, axis=0) / 100
        ds.close()

    return mcs_sum, lat, lon


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
mcs_eur_gpm, lat_eur_gpm, lon_eur_gpm = load_obs('EURR-3')

# Plot figure
fig = plt.figure(figsize=(12, 6))
font_size = 10

orig_cmap = plt.cm.get_cmap('RdYlBu_r', 256)
colors = orig_cmap(np.linspace(0, 1, 256))
colors[0] = [1, 1, 1, 1]
cmap = mcolors.ListedColormap(colors)
mcs_levels = 20 # np.arange(1, 105, 5)

# GPM
ax1 = fig.add_subplot(2, 3, 1, projection=ccrs.PlateCarree())
cf = ax1.contourf(lon_eur_gpm, lat_eur_gpm, mcs_eur_gpm, levels=mcs_levels, cmap=cmap, extend='min', transform=ccrs.PlateCarree())
plt.title('(a) GPM', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax1, lon_car_era5, lat_car_era5)

# CPM eval
ax2 = fig.add_subplot(2, 3, 2, projection=ccrs.PlateCarree())
cf = ax2.contourf(lon_eur_gpm, lat_eur_gpm, mcs_eur_gpm, levels=mcs_levels, cmap=cmap, extend='min', transform=ccrs.PlateCarree())
plt.title('(b) ICTP EURR-3 Evaluation', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2, lon_eur_gpm, lat_eur_gpm)

# CPM hist
ax3 = fig.add_subplot(2, 3, 3, projection=ccrs.PlateCarree())
cf = ax3.contourf(lon_eur_gpm, lat_eur_gpm, mcs_eur_gpm, levels=mcs_levels, cmap=cmap, extend='min', transform=ccrs.PlateCarree())
plt.title('(c) ICTP EURR-3 Historical', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3, lon_eur_gpm, lat_eur_gpm)

# RCM eval
ax4 = fig.add_subplot(2, 3, 5, projection=ccrs.PlateCarree())
cf = ax4.contourf(lon_eur_gpm, lat_eur_gpm, mcs_eur_gpm, levels=mcs_levels, cmap=cmap, extend='min', transform=ccrs.PlateCarree())
plt.title('(d) ICTP EUR-12 Evaluation', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax4, lon_eur_gpm, lat_eur_gpm)

# RCM hist
ax5 = fig.add_subplot(2, 3, 6, projection=ccrs.PlateCarree())
cf = ax5.contourf(lon_eur_gpm, lat_eur_gpm, mcs_eur_gpm, levels=mcs_levels, cmap=cmap, extend='min', transform=ccrs.PlateCarree())
plt.title('(e) ICTP EUR-12 Historical', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax5, lon_eur_gpm, lat_eur_gpm)

cbar_ax = fig.add_axes([0.25, 0.05, 0.5, 0.02])  # [left, bottom, width, height]
cbar = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal')
cbar.set_label('MCS frequency (days/year)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/figs/v2'
name_out = f'pyplt_maps_moaap_mcs_clim_domains_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


