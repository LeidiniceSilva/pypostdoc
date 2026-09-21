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


def load_dataset_annual_mean(path_, start_year, end_year, pattern="*_MOAAP-masks.nc"):
    """Loads datasets for a specific year range and computes the annual mean frequency."""
    data_path = f"/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/dataset/{path_}"
    all_files = sorted(glob.glob(os.path.join(data_path, pattern)))
    
    # Filter files strictly within the year range
    file_list = []
    for f in all_files:
        filename = os.path.basename(f)
        for y in range(start_year, end_year + 1):
            if str(y) in filename:
                file_list.append(f)
                break

    mcs_sum, lat, lon = None, None, None
    for f in tqdm(file_list, desc=f"Loading {path_} ({start_year}-{end_year})"):
        ds = xr.open_dataset(f)
        mcs = ds["MCS_Tb_Objects"].values

        if mcs_sum is None:
            mcs_sum = np.zeros((mcs.shape[1], mcs.shape[2]))
            lat = ds["lat"].values
            lon = ds["lon"].values 

        mcs_sum += np.nansum(mcs, axis=0) / 100
        ds.close()

    # Compute annual mean frequency across the period
    n_years = (end_year - start_year) + 1
    mcs_annual_mean = mcs_sum / n_years

    return mcs_annual_mean, lat, lon


def configure_subplot(ax):
    """Configures geographic projection, ticks, grid, and borders."""
    ax.set_extent([-12, 26, 36, 58], crs=ccrs.PlateCarree())
    xticks = np.linspace(-12, 26, 7)
    yticks = np.linspace(36, 58, 6)

    ax.set_xticks(xticks, crs=ccrs.PlateCarree())
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())

    ax.grid(color='gray', ls='--', alpha=0.5)
    ax.coastlines(linewidth=0.5)
    ax.add_feature(cfeat.BORDERS, linewidth=0.5)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(8)


# Domain
domain = 'EUR'

# Paths
rcm_hist_path = 'RCMs/ICTP/EUR-12/historical/ECEarth/output'
rcm_ssp_path  = 'RCMs/ICTP/EUR-12/ssp370/ECEarth/output'

# 1. Historical (2000-2009)
mcs_hist, lat_hist, lon_hist = load_dataset_annual_mean(rcm_hist_path, 2000, 2009)

# 2. GWL 1.5°C (Combining Historical 2000-2009 and SSP370 2015-2021)
mcs_gwl15_part1, _, _ = load_dataset_annual_mean(rcm_hist_path, 2002, 2009)
mcs_gwl15_part2, _, _ = load_dataset_annual_mean(rcm_ssp_path, 2015, 2021)

# Combined weighted average for GWL 1.5°C (10 yrs + 7 yrs = 17 yrs total)
mcs_gwl15 = (mcs_gwl15_part1 * 10 + mcs_gwl15_part2 * 7) / 17

# 3. GWL 2.0°C (SSP370 2023-2042)
mcs_gwl20, _, _ = load_dataset_annual_mean(rcm_ssp_path, 2023, 2042)

# Compute Percentage Change relative to Historical
# Change (%) = ((GWL - Historical) / Historical) * 100
# Handle zero/small division safely
with np.errstate(divide='ignore', invalid='ignore'):
    change_gwl15 = np.where(mcs_hist > 0.01, ((mcs_gwl15 - mcs_hist) / mcs_hist) * 100, np.nan)
    change_gwl20 = np.where(mcs_hist > 0.01, ((mcs_gwl20 - mcs_hist) / mcs_hist) * 100, np.nan)

# Plot Figure (1 row, 2 columns)
fig = plt.figure(figsize=(12, 4))
font_size = 10

# Diverging color map centered at zero
cmap = plt.cm.get_cmap('BrBG', 40)
diff_levels = np.arange(-100, 105, 5)

# Panel (a): GWL 1.5°C - Historical
ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
cf1 = ax1.contourf(lon_hist, lat_hist, change_gwl15, levels=diff_levels, cmap=cmap, extend='both', transform=ccrs.PlateCarree())
ax1.set_title('(a) GWL 1.5°C- Historical', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

# Panel (b): GWL 2.0°C - Historical
ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())
cf2 = ax2.contourf(lon_hist, lat_hist, change_gwl20, levels=diff_levels, cmap=cmap, extend='both', transform=ccrs.PlateCarree())
ax2.set_title('(b) GWL 2.0° - Historical', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

# Colorbar configuration
cbar_ax = fig.add_axes([0.25, 0.08, 0.5, 0.03])  # [left, bottom, width, height]
cbar = fig.colorbar(cf1, cax=cbar_ax, orientation='horizontal')
cbar.set_label('MCS frequency Relative change (%)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/figs/v2'
name_out = f'pyplt_maps_moaap_mcs_clim_{domain}_GWL_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
