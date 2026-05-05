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
mcs_car, lat_car, lon_car = load_era5('CAR-4')
lon_car = ((lon_car + 180) % 360) - 180

mcs_csam, lat_csam, lon_csam = load_era5('CSAM-3')
lon_csam = ((lon_csam + 180) % 360) - 180

mcs_eurr, lat_eurr, lon_eurr = load_era5('EURR-3')

# Plot figure
fig = plt.figure(figsize=(6, 12))
font_size = 10

# CAR-4
ax1 = fig.add_subplot(3, 1, 1, projection=ccrs.PlateCarree())
cf = ax1.contourf(lon_car, lat_car, mcs_car, levels=20, cmap="RdYlBu_r", transform=ccrs.PlateCarree())
plt.title('(a)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('CAR-4', fontsize=font_size, fontweight='bold')
configure_subplot(ax1, lon_car, lat_car)

# CSAM-3
ax2 = fig.add_subplot(3, 1, 2, projection=ccrs.PlateCarree())
cf = ax2.contourf(lon_csam, lat_csam, mcs_csam, levels=20, cmap="RdYlBu_r", transform=ccrs.PlateCarree())
plt.title('(b)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('CSAM-3', fontsize=font_size, fontweight='bold')
configure_subplot(ax2, lon_csam, lat_csam)

# EURR-3 
ax3 = fig.add_subplot(3, 1, 3, projection=ccrs.PlateCarree())
cf = ax3.contourf(lon_eurr, lat_eurr, mcs_eurr, levels=20, cmap="RdYlBu_r", transform=ccrs.PlateCarree())
plt.title('(c)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('EURR-3', fontsize=font_size, fontweight='bold')
configure_subplot(ax3, lon_eurr, lat_eurr)

cbar = fig.colorbar(cf, ax=[ax1, ax2, ax3], orientation='horizontal', fraction=0.05, pad=0.04, aspect=40)
cbar.set_label('MCS frequency (days/year)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/figs'
name_out = f'pyplt_maps_moaap_mcs_clim_domains_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


