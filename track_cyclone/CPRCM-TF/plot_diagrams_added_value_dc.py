# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Sep 22, 2025"
__description__ = "This script plots Hovmöller diagrams"

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

path = "/leonardo/home/userexternal/mdasilva/leonardo_work/TRACK-CYCLONE/CPRCM-TF"
font_size = 9
seas = 'DJF_2000-2009'


def fmt_lon(x):
    if x < 0:
        return f"{abs(x):.0f}°W"
    elif x > 0:
        return f"{x:.0f}°E"
    else:
        return "0°"


def fmt_lat(x):
    if x < 0:
        return f"{abs(x):.0f}°S"
    elif x > 0:
        return f"{x:.0f}°N"
    else:
        return "0°"


def import_hovmoller_data(param, dataset, lat_band, lon_band):
    arq = f"{path}/data/Dataset/{param}_{dataset}_mmhr_{seas}_dc.nc"
    ds = xr.open_dataset(arq, chunks={"time": -1})

    lat_key = next((k for k in ['lat', 'latitude', 'XLAT', 'xlat', 'rlat'] if k in ds.variables), 'lat')
    lon_key = next((k for k in ['lon', 'longitude', 'XLONG', 'xlon', 'rlon'] if k in ds.variables), 'lon')

    lon_vals_all = ds[lon_key].values
    lat_vals_all = ds[lat_key].values

    # Automatically normalize 0-360 datasets to -180 to 180 so negative bands work universally
    if np.nanmax(lon_vals_all) > 180:
        if lon_vals_all.ndim == 1:
            ds = ds.assign_coords({lon_key: ((ds[lon_key] + 180) % 360) - 180})
            ds = ds.sortby(lon_key)
            lon_vals_all = ds[lon_key].values
            precip = ds[param].values
        else:
            # 2D curvilinear grid conversion
            lon_vals_all = np.where(lon_vals_all > 180, lon_vals_all - 360, lon_vals_all)
            precip = ds[param].values
    else:
        precip = ds[param].values

    lat_min, lat_max = lat_band
    lon_min, lon_max = lon_band

    # Files are already the diurnal cycle with 24 time steps
    unique_hours = np.arange(24)

    if lat_vals_all.ndim == 2:
        # 2D coordinate grids (e.g. WRF/CPRCM)
        lat_mask = (lat_vals_all >= lat_min) & (lat_vals_all <= lat_max)
        lon_mask = (lon_vals_all >= lon_min) & (lon_vals_all <= lon_max)

        lon_1d = np.nanmean(np.where(lat_mask, lon_vals_all, np.nan), axis=0)
        lat_1d = np.nanmean(np.where(lon_mask, lat_vals_all, np.nan), axis=1)

        lat_mask_t = lat_mask[np.newaxis, :, :]
        precip_lat_masked = np.where(lat_mask_t, precip, np.nan)
        lon_time_vals = np.nanmean(precip_lat_masked, axis=1)

        lon_mask_t = lon_mask[np.newaxis, :, :]
        precip_lon_masked = np.where(lon_mask_t, precip, np.nan)
        lat_time_vals = np.nanmean(precip_lon_masked, axis=2)

    else:
        # Standard 1D coordinate grids (handling descending latitude correctly)
        lat_vals_all = lat_vals_all.astype(float)
        lon_vals_all = lon_vals_all.astype(float)

        lat_mask = (lat_vals_all >= lat_min) & (lat_vals_all <= lat_max)
        lon_mask = (lon_vals_all >= lon_min) & (lon_vals_all <= lon_max)

        lon_1d = lon_vals_all[lon_mask]
        lat_1d = lat_vals_all[lat_mask]

        precip_region = precip[:, lat_mask, :]
        precip_region = precip_region[:, :, lon_mask]  # shape: (time, lat_sub, lon_sub)

        # Average over latitude band (axis 1) and longitude band (axis 2) directly across the 24 time steps
        lon_time_vals = np.nanmean(precip_region, axis=1)
        lat_time_vals = np.nanmean(precip_region, axis=2)

    ds.close()
    return unique_hours, lon_1d, lat_1d, lon_time_vals, lat_time_vals


# Define averaging bands for regional analysis using standard negative coordinates
lat_averaging_band = (-35.0, -13.0)   # Latitude bounds
lon_averaging_band = (-60.0, -40.0)   # Longitude bounds (-60° to -40°)

# Import data for all 4 configurations
hours, lons_era5, lats_era5, lon_t_era5, lat_t_era5 = import_hovmoller_data(
    'tp', 'CSAM-3_ERA5', lat_averaging_band, lon_averaging_band)

_,     lons_gpm,  lats_gpm,  lon_t_gpm,  lat_t_gpm  = import_hovmoller_data(
    'pr', 'CSAM-3_GPM_IMERG', lat_averaging_band, lon_averaging_band)

_,     lons_cpr,  lats_cpr,  lon_t_cpr,  lat_t_cpr  = import_hovmoller_data(
    'pr', 'CSAM-3_ERA5_evaluation_RegCM5', lat_averaging_band, lon_averaging_band)

_,     lons_rcm,  lats_rcm,  lon_t_rcm,  lat_t_rcm  = import_hovmoller_data(
    'pr', 'SAM-12_ERA5_evaluation_RegCM5', lat_averaging_band, lon_averaging_band)

# Set up the 2x4 Subplot Figure Layout
fig, axes = plt.subplots(2, 4, figsize=(16, 7))
(ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8) = axes

levels = np.linspace(0.0, 1.0, 22)
kw = dict(levels=levels, cmap="terrain_r", extend="max", rasterized=True)

# --- Top Row: Longitude-Time Sections (X: Longitude, Y: Hour) ---
cf1 = ax1.contourf(lons_era5, hours, lon_t_era5, **kw)
ax1.set_title("(a) Reanalysis (25 km)", loc="left", fontsize=font_size, fontweight="bold")

cf2 = ax2.contourf(lons_gpm, hours, lon_t_gpm, **kw)
ax2.set_title("(b) Satellite estimation (10 km)", loc="left", fontsize=font_size, fontweight="bold")

cf3 = ax3.contourf(lons_cpr, hours, lon_t_cpr, **kw)
ax3.set_title("(c) CPRCM (3 km)", loc="left", fontsize=font_size, fontweight="bold")

cf4 = ax4.contourf(lons_rcm, hours, lon_t_rcm, **kw)
ax4.set_title("(d) RCM (12 km)", loc="left", fontsize=font_size, fontweight="bold")

for ax in [ax1, ax2, ax3, ax4]:
    ax.set_ylabel("Hour (LST)", fontsize=font_size, fontweight="bold")
    ax.set_yticks([0, 6, 12, 18, 23])
    ax.set_xlim(lon_averaging_band[0], lon_averaging_band[1])   # -60 to -40
    ax.set_xticks(np.arange(lon_averaging_band[0], lon_averaging_band[1] + 1, 5))
    ax.set_xticklabels([fmt_lon(x) for x in np.arange(lon_averaging_band[0], lon_averaging_band[1] + 1, 5)])
    ax.grid(True, linestyle="--", alpha=0.5)

# --- Bottom Row: Latitude-Time Sections (X: Latitude, Y: Hour) ---
cf5 = ax5.contourf(lats_era5, hours, lat_t_era5, **kw)
ax5.set_title("(e) Reanalysis (25 km)", loc="left", fontsize=font_size, fontweight="bold")

cf6 = ax6.contourf(lats_gpm, hours, lat_t_gpm, **kw)
ax6.set_title("(f) Satellite estimation (10 km)", loc="left", fontsize=font_size, fontweight="bold")

cf7 = ax7.contourf(lats_cpr, hours, lat_t_cpr, **kw)
ax7.set_title("(g) CPRCM (3 km)", loc="left", fontsize=font_size, fontweight="bold")

cf8 = ax8.contourf(lats_rcm, hours, lat_t_rcm, **kw)
ax8.set_title("(h) RCM (12 km)", loc="left", fontsize=font_size, fontweight="bold")

for ax in [ax5, ax6, ax7, ax8]:
    ax.set_ylabel("Hour (LST)", fontsize=font_size, fontweight="bold")
    ax.set_yticks([0, 6, 12, 18, 23])
    ax.set_xlim(lat_averaging_band[0], lat_averaging_band[1])   # -35 to -13
    ax.set_xticks(np.arange(lat_averaging_band[0], lat_averaging_band[1] + 1, 5))
    ax.set_xticklabels([fmt_lat(x) for x in np.arange(lat_averaging_band[0], lat_averaging_band[1] + 1, 5)])
    ax.grid(True, linestyle="--", alpha=0.5)

# Shared Colorbar
cbar_ax = fig.add_axes([0.91, 0.15, 0.012, 0.7])
cbar = fig.colorbar(cf8, cax=cbar_ax)
cbar.set_label("Precipitation amount (mm hr$^{-1}$)", fontsize=font_size, fontweight="bold")
cbar.ax.tick_params(labelsize=font_size)

# Save figure
path_out = f'{path}/figs/'
os.makedirs(path_out, exist_ok=True)
name_out = f'pyplt_diagrams_diurnal_cycle_{seas}.png'
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
