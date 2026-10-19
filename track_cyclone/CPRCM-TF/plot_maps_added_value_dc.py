# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Sep 22, 2025"
__description__ = "This script plots added value"

import os
import netCDF4
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.colors
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat

from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

path = "/leonardo/home/userexternal/mdasilva/leonardo_work/TRACK-CYCLONE/CPRCM-TF"
font_size = 9
seas = 'JJA_2000-2009'


def import_diurnal_data(param, dataset):

    arq = "{0}/data/Dataset/{1}_{2}_mmhr_{3}.nc".format(path, param, dataset, seas)
    data = netCDF4.Dataset(arq)

    lat_key = next((k for k in ['lat', 'latitude', 'XLAT', 'xlat', 'rlat'] if k in data.variables), 'lat')
    lon_key = next((k for k in ['lon', 'longitude', 'XLONG', 'xlon', 'rlon'] if k in data.variables), 'lon')
    lat = data.variables[lat_key][:]
    lon = data.variables[lon_key][:]
    data.close()

    ds = xr.open_dataset(arq)
    precip = ds[param]

    # Mean diurnal cycle: shape (24, lat, lon)
    diurnal_cycle = precip.groupby("time.hour").mean(dim="time")
    hours = diurnal_cycle["hour"].values  # shape (N,)
    N = len(hours)

    omega = 2.0 * np.pi / 24.0

    # First-harmonic Fourier coefficients
    h = xr.DataArray(hours, dims=["hour"])
    a1 = (2.0 / N) * (diurnal_cycle * np.cos(omega * h)).sum(dim="hour")
    b1 = (2.0 / N) * (diurnal_cycle * np.sin(omega * h)).sum(dim="hour")

    # Amplitude of the first harmonic
    amplitude = np.sqrt(a1**2 + b1**2)

    # Phase -> peak hour in UTC
    phase_rad = np.arctan2(b1, a1)
    peak_hour_utc = (phase_rad / omega) % 24

    # --- FIX: handle both 1D and 2D coordinate grids ---
    spatial_dims = diurnal_cycle.dims[1:]  # e.g. ('lat', 'lon') or ('y', 'x')

    if lon.ndim == 2:
        # Curvilinear grid (e.g. WRF/RegCM5): lat/lon already 2D, use as-is
        lon_2d = xr.DataArray(np.array(lon), dims=spatial_dims)
        lat_2d = xr.DataArray(np.array(lat), dims=spatial_dims)
    else:
        # Regular grid: broadcast 1D lat/lon to 2D
        lat_1d = np.array(lat)
        lon_1d = np.array(lon)
        lon_2d = xr.DataArray(
            np.broadcast_to(lon_1d[np.newaxis, :], (len(lat_1d), len(lon_1d))).copy(),
            dims=spatial_dims
        )
        lat_2d = xr.DataArray(
            np.broadcast_to(lat_1d[:, np.newaxis], (len(lat_1d), len(lon_1d))).copy(),
            dims=spatial_dims
        )

    # Convert UTC peak hour to Local Solar Time (LST)
    time_offset = lon_2d / 15.0  # degrees -> hours
    peak_hour_lst = (peak_hour_utc + time_offset) % 24

    ds.close()
    return lat_2d.values, lon_2d.values, amplitude.values, peak_hour_lst.values


# --- Import data for all 4 configurations ---
lat_cprcm,  lon_cprcm,  amp_cprcm,  peak_hr_cprcm = import_diurnal_data('pr', 'CSAM-3_ERA5_evaluation_RegCM5')
lat_rcm,    lon_rcm,    amp_rcm,    peak_hr_rcm   = import_diurnal_data('pr', 'SAM-12_ERA5_evaluation_RegCM5')
lat_gpm,    lon_gpm,    amp_gpm,    peak_hr_gpm   = import_diurnal_data('pr', 'CSAM-3_GPM_IMERG')
lat_era5,   lon_era5,   amp_era5,   peak_hr_era5  = import_diurnal_data('tp', 'CSAM-3_ERA5')

print(amp_cprcm.shape)
print(np.min(amp_cprcm), np.max(amp_cprcm))
print()
print(peak_hr_cprcm.shape)
print(np.min(peak_hr_cprcm), np.max(peak_hr_cprcm))


def configure_subplot(ax):
    ax.set_extent([-66, -38, -36, -14], crs=ccrs.PlateCarree())
    ax.set_xticks(np.arange(-66, -38, 4), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(-36, -14, 4), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.grid(c='k', ls='--', alpha=0.5)
    ax.add_feature(cfeat.BORDERS, linewidth=0.75, color='black')
    ax.coastlines(linewidth=0.75, color='black')


# Plot figure
fig, axes = plt.subplots(2, 4, figsize=(16, 7), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8) = axes

levels_time = np.linspace(0, 24, 25)
levels_amp  = np.linspace(0, 1.5, 16)

cf1 = ax1.contourf(lon_era5,   lat_era5,   peak_hr_era5,   levels=levels_time, cmap="hsv",       extend="neither", transform=ccrs.PlateCarree())
ax1.set_title("(a) Reanalysis (25 km)",           loc="left", fontsize=font_size, fontweight="bold")
configure_subplot(ax1)

cf2 = ax2.contourf(lon_gpm,    lat_gpm,    peak_hr_gpm,    levels=levels_time, cmap="hsv",       extend="neither", transform=ccrs.PlateCarree())
ax2.set_title("(b) Satellite estimation (10 km)", loc="left", fontsize=font_size, fontweight="bold")
configure_subplot(ax2)

cf3 = ax3.contourf(lon_cprcm,  lat_cprcm,  peak_hr_cprcm,  levels=levels_time, cmap="hsv",       extend="neither", transform=ccrs.PlateCarree())
ax3.set_title("(c) CPRCM (3 km)",                loc="left", fontsize=font_size, fontweight="bold")
configure_subplot(ax3)

cf4 = ax4.contourf(lon_rcm,    lat_rcm,    peak_hr_rcm,    levels=levels_time, cmap="hsv",       extend="neither", transform=ccrs.PlateCarree())
ax4.set_title("(d) RCM (12 km)",                 loc="left", fontsize=font_size, fontweight="bold")
configure_subplot(ax4)

cbar1 = plt.colorbar(cf1, cax=fig.add_axes([0.91, 0.53, 0.012, 0.35]), pad=0.02, shrink=0.9, ticks=[0, 3, 6, 9, 12, 15, 18, 21, 24])
cbar1.set_label("Time of diurnal cycle maximum (hour)", fontsize=font_size, fontweight="bold")
cbar1.ax.tick_params(labelsize=font_size)

cf5 = ax5.contourf(lon_era5,   lat_era5,   amp_era5,   levels=levels_amp, cmap="terrain_r", extend="neither", transform=ccrs.PlateCarree())
ax5.set_title("(e) Reanalysis (25 km)",           loc="left", fontsize=font_size, fontweight="bold")
configure_subplot(ax5)

cf6 = ax6.contourf(lon_gpm,    lat_gpm,    amp_gpm,    levels=levels_amp, cmap="terrain_r", extend="neither", transform=ccrs.PlateCarree())
ax6.set_title("(f) Satellite estimation (10 km)", loc="left", fontsize=font_size, fontweight="bold")
configure_subplot(ax6)

cf7 = ax7.contourf(lon_cprcm,  lat_cprcm,  amp_cprcm,  levels=levels_amp, cmap="terrain_r", extend="neither", transform=ccrs.PlateCarree())
ax7.set_title("(g) CPRCM (3 km)",                loc="left", fontsize=font_size, fontweight="bold")
configure_subplot(ax7)

cf8 = ax8.contourf(lon_rcm,    lat_rcm,    amp_rcm,    levels=levels_amp, cmap="terrain_r", extend="neither", transform=ccrs.PlateCarree())
ax8.set_title("(h) RCM (12 km)",                 loc="left", fontsize=font_size, fontweight="bold")
configure_subplot(ax8)

cbar2 = plt.colorbar(cf8, cax=fig.add_axes([0.91, 0.1, 0.012, 0.35]), pad=0.02, shrink=0.9)
cbar2.set_label("Amplitude (mm hr$^{-1}$)", fontsize=font_size, fontweight="bold")
cbar2.ax.tick_params(labelsize=font_size)

# Save figure
path_out = f'{path}/figs/'
os.makedirs(path_out, exist_ok=True)
name_out = f'pyplt_maps_diurnal_cycle_CSAM-3_{seas}.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
