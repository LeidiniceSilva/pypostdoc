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

from tqdm import tqdm
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

def load_era5(domain, pattern="*_MOAAP-masks.nc"):

    data_path = f"/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/ERA5/{domain}/output/"
    file_list = sorted(glob.glob(os.path.join(data_path, pattern)))

    mcs_sum = None
    lat, lon = None, None

    for f in tqdm(file_list):

        ds = xr.open_dataset(f)
        mcs = ds["MCS_Tb_Objects"].values  # (time, lat, lon)

        if mcs_sum is None:
            mcs_sum = np.zeros((mcs.shape[1], mcs.shape[2]))
            lat = ds["lat"].values
            lon = ds["lon"].values 

        # climatological occurrence
        mcs_sum += np.nansum(mcs, axis=0) / 10

        ds.close()

    return mcs_sum, lat, lon


def load_cpm(domain, pattern="*_MOAAP-masks.nc"):

    data_path = f"/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/ERA5/{domain}/output/"
    file_list = sorted(glob.glob(os.path.join(data_path, pattern)))

    mcs_sum = None
    lat, lon = None, None

    for f in tqdm(file_list):

        ds = xr.open_dataset(f)
        mcs = ds["MCS_Tb_Objects"].values  # (time, lat, lon)

        if mcs_sum is None:
            mcs_sum = np.zeros((mcs.shape[1], mcs.shape[2]))
            lat = ds["lat"].values
            lon = ds["lon"].values 

        # climatological occurrence
        mcs_sum += np.nansum(mcs, axis=0) / 10

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
mcs_car_era5, lat_car_era5, lon_car_era5 = load_era5('CAR-4')
lon_car_era5 = ((lon_car_era5 + 180) % 360) - 180

mcs_csam_era5, lat_csam_era5, lon_csam_era5 = load_era5('CSAM-3')
lon_csam_era5 = ((lon_csam_era5 + 180) % 360) - 180

mcs_eurr_era5, lat_eurr_era5, lon_eurr_era5 = load_era5('EURR-3')

mcs_car_cpm, lat_car_cpm, lon_car_cpm = load_cpm('CAR-4')
lon_car_cpm = ((lon_car_cpm + 180) % 360) - 180

mcs_csam_cpm, lat_csam_cpm, lon_csam_cpm = load_cpm('CSAM-3')
lon_csam_cpm = ((lon_csam_cpm + 180) % 360) - 180

mcs_eurr_cpm, lat_eurr_cpm, lon_eurr_cpm = load_cpm('EURR-3')

# Plot figure
fig = plt.figure(figsize=(10, 8))
font_size = 10

# CAR-4
ax1 = fig.add_subplot(3, 2, 1, projection=ccrs.PlateCarree())
cf = ax1.contourf(lon_car_era5, lat_car_era5, mcs_car_era5, levels=20, cmap="RdYlBu_r", transform=ccrs.PlateCarree())
plt.title('(a)', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax1, lon_car_era5, lat_car_era5)

ax2 = fig.add_subplot(3, 2, 2, projection=ccrs.PlateCarree())
cf = ax2.contourf(lon_car_cpm, lat_car_cpm, mcs_car_cpm, levels=20, cmap="RdYlBu_r", transform=ccrs.PlateCarree())
plt.title('(b)', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2, lon_car_cpm, lat_car_cpm)

# CSAM-3
ax3 = fig.add_subplot(3, 2, 3, projection=ccrs.PlateCarree())
cf = ax3.contourf(lon_csam_era5, lat_csam_era5, mcs_csam_era5, levels=20, cmap="RdYlBu_r", transform=ccrs.PlateCarree())
plt.title('(c)', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3, lon_csam_era5, lat_csam_era5)

ax4 = fig.add_subplot(3, 2, 4, projection=ccrs.PlateCarree())
cf = ax4.contourf(lon_csam_cpm, lat_csam_cpm, mcs_csam_cpm, levels=20, cmap="RdYlBu_r", transform=ccrs.PlateCarree())
plt.title('(d)', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax4, lon_csam_cpm, lat_csam_cpm)

# EURR-3 
ax5 = fig.add_subplot(3, 2, 5, projection=ccrs.PlateCarree())
cf = ax5.contourf(lon_eurr_era5, lat_eurr_era5, mcs_eurr_era5, levels=20, cmap="RdYlBu_r", transform=ccrs.PlateCarree())
plt.title('(e)', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax5, lon_eurr_era5, lat_eurr_era5)

ax6 = fig.add_subplot(3, 2, 6, projection=ccrs.PlateCarree())
cf = ax6.contourf(lon_eurr_cpm, lat_eurr_cpm, mcs_eurr_cpm, levels=20, cmap="RdYlBu_r", transform=ccrs.PlateCarree())
plt.title('(f)', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax6, lon_eurr_cpm, lat_eurr_cpm)

cbar_ax = fig.add_axes([0.25, 0.08, 0.5, 0.02])  # [left, bottom, width, height]
cbar = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal')
cbar.set_label('MCS frequency (days/year)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/figs'
name_out = f'pyplt_maps_moaap_mcs_clim_domains_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


