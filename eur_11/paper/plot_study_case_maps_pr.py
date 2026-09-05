# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 28, 2026"
__description__ = "This script plot the study case maps"

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


def import_regcm5(path, param, exp):

    arq = f"{path}/RegCM5/{param}_{exp}_day_2006sep26.nc"
    files = glob.glob(arq)

    if len(files) == 0:
        raise FileNotFoundError(f"No file found: {arq}")

    data = netCDF4.Dataset(files[0])
    var = np.squeeze(data.variables[param][:])

    if np.ma.isMaskedArray(var):
        var = var.filled(np.nan)

    lat = np.squeeze(data.variables["xlat"][:])
    lon = np.squeeze(data.variables["xlon"][:])

    data.close()

    return var, lat, lon


def import_obs(path, param, obs_name):

    arq = f"{path}/ERA5/{param}_{obs_name}_day_2006sep26.nc"
    files = glob.glob(arq)

    if len(files) == 0:
        raise FileNotFoundError(f"No file found: {arq}")

    data = netCDF4.Dataset(files[0])
    var = np.squeeze(data.variables[param][:])

    if np.ma.isMaskedArray(var):
        var = var.filled(np.nan)

    lat = np.squeeze(data.variables["latitude"][:])
    lon = np.squeeze(data.variables["longitude"][:])

    data.close()

    return var, lat, lon


def import_regcm5_pressure(path, exp, lon_min, lon_max, lat_min, lat_max):

    arq = f"{path}/RegCM5/psl_{exp}_1hr_2006sep26.nc"
    files = glob.glob(arq)

    if len(files) == 0:
        raise FileNotFoundError(f"No file found: {arq}")

    data = netCDF4.Dataset(files[0])
    psl = np.squeeze(data.variables["psl"][:])

    if np.ma.isMaskedArray(psl):
        psl = psl.filled(np.nan)

    lat = np.squeeze(data.variables["xlat"][:])
    lon = np.squeeze(data.variables["xlon"][:])

    data.close()

    # Convert Pa to hPa if necessary
    if np.nanmean(psl) > 10000:
        psl = psl / 100.0

    # Mask outside spatial extent to find local domain minimum pressure and time
    mask = (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)
    psl_domain = np.where(mask[np.newaxis, :, :], psl, np.nan)

    # Calculate spatial minimum at each hourly time step
    min_per_time = np.nanmin(psl_domain, axis=(1, 2))
    min_time_idx = int(np.argmin(min_per_time))

    slp_min_field = psl[min_time_idx, :, :]

    return slp_min_field, min_time_idx, lat, lon


def import_obs_pressure(path, obs_name, lon_min, lon_max, lat_min, lat_max):

    arq = f"{path}/ERA5/msl_{obs_name}_1hr_2006sep26.nc"
    files = glob.glob(arq)

    if len(files) == 0:
        raise FileNotFoundError(f"No file found: {arq}")

    data = netCDF4.Dataset(files[0])
    msl = np.squeeze(data.variables["msl"][:])

    if np.ma.isMaskedArray(msl):
        msl = msl.filled(np.nan)

    lat = np.squeeze(data.variables["latitude"][:])
    lon = np.squeeze(data.variables["longitude"][:])

    data.close()

    # Convert Pa to hPa if necessary
    if np.nanmean(msl) > 10000:
        msl = msl / 100.0

    lon_2d, lat_2d = np.meshgrid(lon, lat)

    # Mask outside spatial extent to find local domain minimum pressure and time
    mask = (lon_2d >= lon_min) & (lon_2d <= lon_max) & (lat_2d >= lat_min) & (lat_2d <= lat_max)
    msl_domain = np.where(mask[np.newaxis, :, :], msl, np.nan)

    # Calculate spatial minimum at each hourly time step
    min_per_time = np.nanmin(msl_domain, axis=(1, 2))
    min_time_idx = int(np.argmin(min_per_time))

    slp_min_field = msl[min_time_idx, :, :]

    return slp_min_field, min_time_idx, lat, lon


def main():

    path = pathlib.Path("/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11/postproc/paper/cyc")

    obs_name = "ERA5"
    dt = "2006sep26"

    # Get spatial limits
    lon_min, lon_max = -10, 30
    lat_min, lat_max = 30, 60

    # Import precipitation data
    pr_obs, lat_obs, lon_obs = import_obs(path, "tp", obs_name)
    pr_exp_1, lat, lon = import_regcm5(path, "pr", "NoTo-EUR")
    pr_exp_2, lat, lon = import_regcm5(path, "pr", "WSM5-EUR")
    pr_exp_3, lat, lon = import_regcm5(path, "pr", "WSM7-EUR")
    pr_exp_4, lat, lon = import_regcm5(path, "pr", "WDM7-EUR")

    # Import pressure fields and find time of minimum pressure
    slp_obs, time_obs, lat_slp_obs, lon_slp_obs = import_obs_pressure(path, obs_name, lon_min, lon_max, lat_min, lat_max)
    slp_exp_1, time_exp_1, _, _ = import_regcm5_pressure(path, "NoTo-EUR", lon_min, lon_max, lat_min, lat_max)
    slp_exp_2, time_exp_2, _, _ = import_regcm5_pressure(path, "WSM5-EUR", lon_min, lon_max, lat_min, lat_max)
    slp_exp_3, time_exp_3, _, _ = import_regcm5_pressure(path, "WSM7-EUR", lon_min, lon_max, lat_min, lat_max)
    slp_exp_4, time_exp_4, _, _ = import_regcm5_pressure(path, "WDM7-EUR", lon_min, lon_max, lat_min, lat_max)

    # Apply units conversion
    pr_obs = pr_obs * 1000.0
    pr_exp_1 = pr_exp_1 * 3600.0
    pr_exp_2 = pr_exp_2 * 3600.0
    pr_exp_3 = pr_exp_3 * 3600.0
    pr_exp_4 = pr_exp_4 * 3600.0

    # Define contour levels and colormap
    levels = np.linspace(0, 50, 17)
    color = ['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff',\
             '#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff',\
             '#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff',\
             '#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']
    cmap = matplotlib.colors.ListedColormap(color)
    norm = matplotlib.colors.BoundaryNorm(levels, cmap.N)

    # Define pressure contour levels
    p_levels = np.arange(970, 1040, 4)

    # Plot figure setup
    fig = plt.figure(figsize=(12, 6))
    font_size = 10

    # Panel 1: ERA5
    ax1 = fig.add_subplot(2, 3, 1, projection=ccrs.PlateCarree())
    ax1.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax1.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax1.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle=':')
    cf1 = ax1.contourf(lon_obs, lat_obs, pr_obs, levels=levels, cmap=cmap, norm=norm, extend='max', transform=ccrs.PlateCarree())
    cs1 = ax1.contour(lon_slp_obs, lat_slp_obs, slp_obs, levels=p_levels, colors='black', linewidths=0.7, transform=ccrs.PlateCarree())
    ax1.clabel(cs1, cs1.levels, inline=True, fmt='%d', fontsize=7)
    ax1.set_title(f"(a) {obs_name} ({time_obs:02d} UTC)", loc='left', fontsize=font_size, fontweight='bold')
    gl1 = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl1.top_labels = False
    gl1.right_labels = False

    # Panel 2: NoTo-EUR
    ax2 = fig.add_subplot(2, 3, 2, projection=ccrs.PlateCarree())
    ax2.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax2.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax2.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle=':')
    cf2 = ax2.contourf(lon, lat, pr_exp_1, levels=levels, cmap=cmap, norm=norm, extend='max', transform=ccrs.PlateCarree())
    cs2 = ax2.contour(lon, lat, slp_exp_1, levels=p_levels, colors='black', linewidths=0.7, transform=ccrs.PlateCarree())
    ax2.clabel(cs2, cs2.levels, inline=True, fmt='%d', fontsize=7)
    ax2.set_title(f"(b) NoTo ({time_exp_1:02d} UTC)", loc='left', fontsize=font_size, fontweight='bold')
    gl2 = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl2.top_labels = False
    gl2.right_labels = False

    # Panel 3: WSM5-EUR
    ax3 = fig.add_subplot(2, 3, 3, projection=ccrs.PlateCarree())
    ax3.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax3.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax3.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle=':')
    cf3 = ax3.contourf(lon, lat, pr_exp_2, levels=levels, cmap=cmap, norm=norm, extend='max', transform=ccrs.PlateCarree())
    cs3 = ax3.contour(lon, lat, slp_exp_2, levels=p_levels, colors='black', linewidths=0.7, transform=ccrs.PlateCarree())
    ax3.clabel(cs3, cs3.levels, inline=True, fmt='%d', fontsize=7)
    ax3.set_title(f"(c) WSM5 ({time_exp_2:02d} UTC)", loc='left', fontsize=font_size, fontweight='bold')
    gl3 = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl3.top_labels = False
    gl3.right_labels = False

    # Panel 4: WSM7-EUR
    ax4 = fig.add_subplot(2, 3, 4, projection=ccrs.PlateCarree())
    ax4.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax4.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax4.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle=':')
    cf4 = ax4.contourf(lon, lat, pr_exp_3, levels=levels, cmap=cmap, norm=norm, extend='max', transform=ccrs.PlateCarree())
    cs4 = ax4.contour(lon, lat, slp_exp_3, levels=p_levels, colors='black', linewidths=0.7, transform=ccrs.PlateCarree())
    ax4.clabel(cs4, cs4.levels, inline=True, fmt='%d', fontsize=7)
    ax4.set_title(f"(d) WSM7 ({time_exp_3:02d} UTC)", loc='left', fontsize=font_size, fontweight='bold')
    gl4 = ax4.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl4.top_labels = False
    gl4.right_labels = False

    # Panel 5: WDM7-EUR
    ax5 = fig.add_subplot(2, 3, 5, projection=ccrs.PlateCarree())
    ax5.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax5.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax5.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle=':')
    cf5 = ax5.contourf(lon, lat, pr_exp_4, levels=levels, cmap=cmap, norm=norm, extend='max', transform=ccrs.PlateCarree())
    cs5 = ax5.contour(lon, lat, slp_exp_4, levels=p_levels, colors='black', linewidths=0.7, transform=ccrs.PlateCarree())
    ax5.clabel(cs5, cs5.levels, inline=True, fmt='%d', fontsize=7)
    ax5.set_title(f"(e) WDM7 ({time_exp_4:02d} UTC)", loc='left', fontsize=font_size, fontweight='bold')
    gl5 = ax5.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl5.top_labels = False
    gl5.right_labels = False

    # Panel 6: Empty
    ax6 = fig.add_subplot(2, 3, 6)
    ax6.axis('off')

    # Colorbar
    cbar_ax = fig.add_axes([0.15, 0.03, 0.7, 0.02])
    cbar = fig.colorbar(cf1, cax=cbar_ax, orientation='horizontal', label='Precipitation acc (mm/d)')

    # Path out to save figure
    path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11/figs/paper'
    name_out = f'pyplt_maps_study_case_RegCM5_EUR-11_pr_{dt}.png'
    plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    main()
