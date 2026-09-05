# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 28, 2026"
__description__ = "This script plot the cyclone track"

import pathlib
import glob
import os
import string
import numpy as np
import netCDF4
import matplotlib.colors
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def import_regcm5_track(path, exp, lon_min, lon_max, lat_min, lat_max):

    # Import psl for tracking position
    arq_psl = f"{path}/RegCM5/psl_{exp}_1hr_2006sep26.nc"
    files_psl = glob.glob(arq_psl)
    if len(files_psl) == 0:
        raise FileNotFoundError(f"No file found: {arq_psl}")

    data_psl = netCDF4.Dataset(files_psl[0])
    psl = np.squeeze(data_psl.variables["psl"][:])
    if np.ma.isMaskedArray(psl):
        psl = psl.filled(np.nan)
    lat = np.squeeze(data_psl.variables["xlat"][:])
    lon = np.squeeze(data_psl.variables["xlon"][:])
    data_psl.close()

    if np.nanmean(psl) > 10000:
        psl = psl / 100.0

    mask = (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)
    psl_domain = np.where(mask[np.newaxis, :, :], psl, np.nan)

    # Calculate track coordinates and minimum MSLP value at each hourly timestep
    track_lat = []
    track_lon = []
    min_mslp = []

    for t in range(psl_domain.shape[0]):
        field_t = psl_domain[t, :, :]
        if not np.all(np.isnan(field_t)):
            min_val = np.nanmin(field_t)
            min_idx = np.unravel_index(np.nanargmin(field_t), field_t.shape)
            track_lat.append(lat[min_idx])
            track_lon.append(lon[min_idx])
            min_mslp.append(min_val)
        else:
            track_lat.append(np.nan)
            track_lon.append(np.nan)
            min_mslp.append(np.nan)

    return np.array(track_lon), np.array(track_lat), np.array(min_mslp)


def import_obs_track(path, obs_name, lon_min, lon_max, lat_min, lat_max):

    # Import msl for tracking position
    arq_msl = f"{path}/ERA5/msl_{obs_name}_1hr_2006sep26.nc"
    files_msl = glob.glob(arq_msl)
    if len(files_msl) == 0:
        raise FileNotFoundError(f"No file found: {arq_msl}")

    data_msl = netCDF4.Dataset(files_msl[0])
    msl = np.squeeze(data_msl.variables["msl"][:])
    if np.ma.isMaskedArray(msl):
        msl = msl.filled(np.nan)
    lat = np.squeeze(data_msl.variables["latitude"][:])
    lon = np.squeeze(data_msl.variables["longitude"][:])
    data_msl.close()

    # Convert longitude from [0, 360] to [-180, 180] and re-sort array
    lon = np.where(lon > 180, lon - 360, lon)
    lon_idx = np.argsort(lon)
    lon = lon[lon_idx]
    msl = msl[:, :, lon_idx]

    if np.nanmean(msl) > 10000:
        msl = msl / 100.0

    lon_2d, lat_2d = np.meshgrid(lon, lat)
    mask = (lon_2d >= lon_min) & (lon_2d <= lon_max) & (lat_2d >= lat_min) & (lat_2d <= lat_max)
    msl_domain = np.where(mask[np.newaxis, :, :], msl, np.nan)

    # Calculate track coordinates and minimum MSLP value at each hourly timestep
    track_lat = []
    track_lon = []
    min_mslp = []

    for t in range(msl_domain.shape[0]):
        field_t = msl_domain[t, :, :]
        if not np.all(np.isnan(field_t)):
            min_val = np.nanmin(field_t)
            min_idx = np.unravel_index(np.nanargmin(field_t), field_t.shape)
            track_lat.append(lat_2d[min_idx])
            track_lon.append(lon_2d[min_idx])
            min_mslp.append(min_val)
        else:
            track_lat.append(np.nan)
            track_lon.append(np.nan)
            min_mslp.append(np.nan)

    return np.array(track_lon), np.array(track_lat), np.array(min_mslp)


def main():

    path = pathlib.Path("/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11/postproc/paper/cyc")

    obs_name = "ERA5"
    dt = "2006sep26"

    # Get spatial limits
    lon_min, lon_max = -10, 20
    lat_min, lat_max = 35, 60

    # Import pressure and cyclone track coordinates
    trk_lon_obs, trk_lat_obs, mslp_obs = import_obs_track(path, obs_name, lon_min, lon_max, lat_min, lat_max)
    trk_lon_exp_1, trk_lat_exp_1, mslp_exp_1 = import_regcm5_track(path, "NoTo-EUR", lon_min, lon_max, lat_min, lat_max)
    trk_lon_exp_2, trk_lat_exp_2, mslp_exp_2 = import_regcm5_track(path, "WSM5-EUR", lon_min, lon_max, lat_min, lat_max)
    trk_lon_exp_3, trk_lat_exp_3, mslp_exp_3 = import_regcm5_track(path, "WSM7-EUR", lon_min, lon_max, lat_min, lat_max)
    trk_lon_exp_4, trk_lat_exp_4, mslp_exp_4 = import_regcm5_track(path, "WDM7-EUR", lon_min, lon_max, lat_min, lat_max)

    hours = np.arange(len(mslp_obs))

    # Figure setup
    fig = plt.figure(figsize=(7, 6))
    font_size = 10

    # Main Plot (a)
    ax_map = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax_map.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # Add features: Light brown land, light blue ocean
    ax_map.add_feature(cfeature.OCEAN, facecolor='lightblue')
    ax_map.add_feature(cfeature.LAND, facecolor='tan', edgecolor='none')
    ax_map.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax_map.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle=':')

    # Plot tracks with slighter/slimmer line widths
    ax_map.plot(trk_lon_obs, trk_lat_obs, color='black', linewidth=1.0, marker='o', markersize=2, label='ERA5', transform=ccrs.PlateCarree())
    ax_map.plot(trk_lon_exp_1, trk_lat_exp_1, color='red', linewidth=1.0, marker='o', markersize=2, label='NoTo', transform=ccrs.PlateCarree())
    ax_map.plot(trk_lon_exp_2, trk_lat_exp_2, color='blue', linewidth=1.0, marker='o', markersize=2, label='WSM5', transform=ccrs.PlateCarree())
    ax_map.plot(trk_lon_exp_3, trk_lat_exp_3, color='green', linewidth=1.0, marker='o', markersize=2, label='WSM7', transform=ccrs.PlateCarree())
    ax_map.plot(trk_lon_exp_4, trk_lat_exp_4, color='purple', linewidth=1.0, marker='o', markersize=2, label='WDM7', transform=ccrs.PlateCarree())

    # Start point markers
    ax_map.plot(trk_lon_obs[0], trk_lat_obs[0], color='black', marker='*', markersize=4, transform=ccrs.PlateCarree())
    ax_map.plot(trk_lon_exp_1[0], trk_lat_exp_1[0], color='red', marker='*', markersize=4, transform=ccrs.PlateCarree())
    ax_map.plot(trk_lon_exp_2[0], trk_lat_exp_2[0], color='blue', marker='*', markersize=4, transform=ccrs.PlateCarree())
    ax_map.plot(trk_lon_exp_3[0], trk_lat_exp_3[0], color='green', marker='*', markersize=4, transform=ccrs.PlateCarree())
    ax_map.plot(trk_lon_exp_4[0], trk_lat_exp_4[0], color='purple', marker='*', markersize=4, transform=ccrs.PlateCarree())

    ax_map.set_title("(a)", loc='left', fontsize=font_size, fontweight='bold')
    ax_map.legend(loc='upper right', fontsize=font_size, framealpha=0.9)

    gl = ax_map.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    # Inset Plot (b)
    ax_inset = fig.add_axes([0.25, 0.5, 0.35, 0.28])

    # Plot MSLP time series with slimmer line widths
    ax_inset.plot(hours, mslp_obs, color='black', linewidth=1.0, marker='o', markersize=2)
    ax_inset.plot(hours, mslp_exp_1, color='red', linewidth=1.0, marker='o', markersize=2)
    ax_inset.plot(hours, mslp_exp_2, color='blue', linewidth=1.0, marker='o', markersize=2)
    ax_inset.plot(hours, mslp_exp_3, color='green', linewidth=1.0, marker='o', markersize=2)
    ax_inset.plot(hours, mslp_exp_4, color='orange', linewidth=1.0, marker='o', markersize=2)

    ax_inset.set_title("(b)", loc='left', fontsize=font_size, fontweight='bold')
    ax_inset.set_xlabel("Hours of the day", fontsize=font_size-1)
    ax_inset.set_ylabel("MSLP (hPa)", fontsize=font_size-1)
    ax_inset.tick_params(axis='both', labelsize=font_size-2)
    ax_inset.grid(True, linestyle='--', alpha=0.75, color='gray')

    # Path out to save figure
    path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11/figs/paper'
    name_out = f'pyplt_maps_study_case_cyclone_track_{dt}.png'
    plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    main()
