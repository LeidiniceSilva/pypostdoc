# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 28, 2026"
__description__ = "This script plot the study case wind speed maps with vectors"

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


def import_regcm5_wind(path, exp, lon_min, lon_max, lat_min, lat_max):

    # Import psl to find time of minimum pressure
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
    min_per_time = np.nanmin(psl_domain, axis=(1, 2))
    min_time_idx = int(np.argmin(min_per_time))

    # Import uas and vas at the selected time step
    arq_uas = f"{path}/RegCM5/uas_{exp}_1hr_2006sep26.nc"
    files_uas = glob.glob(arq_uas)
    if len(files_uas) == 0:
        raise FileNotFoundError(f"No file found: {arq_uas}")

    data_uas = netCDF4.Dataset(files_uas[0])
    uas = np.squeeze(data_uas.variables["uas"][min_time_idx, :, :])
    if np.ma.isMaskedArray(uas):
        uas = uas.filled(np.nan)
    data_uas.close()

    arq_vas = f"{path}/RegCM5/vas_{exp}_1hr_2006sep26.nc"
    files_vas = glob.glob(arq_vas)
    if len(files_vas) == 0:
        raise FileNotFoundError(f"No file found: {arq_vas}")

    data_vas = netCDF4.Dataset(files_vas[0])
    vas = np.squeeze(data_vas.variables["vas"][min_time_idx, :, :])
    if np.ma.isMaskedArray(vas):
        vas = vas.filled(np.nan)
    data_vas.close()

    # Calculate wind speed magnitude (m/s)
    ws = np.sqrt(uas**2 + vas**2)

    return ws, uas, vas, min_time_idx, lat, lon


def import_obs_wind(path, obs_name, lon_min, lon_max, lat_min, lat_max):

    # Import msl to find time of minimum pressure
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
    min_per_time = np.nanmin(msl_domain, axis=(1, 2))
    min_time_idx = int(np.argmin(min_per_time))

    # Import u10 and v10 at the selected time step
    arq_u10 = f"{path}/ERA5/u10_{obs_name}_1hr_2006sep26.nc"
    files_u10 = glob.glob(arq_u10)
    if len(files_u10) == 0:
        raise FileNotFoundError(f"No file found: {arq_u10}")

    data_u10 = netCDF4.Dataset(files_u10[0])
    u10 = np.squeeze(data_u10.variables["u10"][min_time_idx, :, :])
    if np.ma.isMaskedArray(u10):
        u10 = u10.filled(np.nan)
    u10 = u10[:, lon_idx]
    data_u10.close()

    arq_v10 = f"{path}/ERA5/v10_{obs_name}_1hr_2006sep26.nc"
    files_v10 = glob.glob(arq_v10)
    if len(files_v10) == 0:
        raise FileNotFoundError(f"No file found: {arq_v10}")

    data_v10 = netCDF4.Dataset(files_v10[0])
    v10 = np.squeeze(data_v10.variables["v10"][min_time_idx, :, :])
    if np.ma.isMaskedArray(v10):
        v10 = v10.filled(np.nan)
    v10 = v10[:, lon_idx]
    data_v10.close()

    # Calculate wind speed magnitude (m/s)
    ws = np.sqrt(u10**2 + v10**2)

    return ws, u10, v10, min_time_idx, lat, lon


def main():

    path = pathlib.Path("/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11/postproc/paper/cyc")

    obs_name = "ERA5"
    dt = "2006sep26"

    # Get spatial limits
    lon_min, lon_max = -10, 30
    lat_min, lat_max = 30, 60

    # Import wind speed and vector components at minimum MSLP time
    ws_obs, u_obs, v_obs, time_obs, lat_obs, lon_obs = import_obs_wind(path, obs_name, lon_min, lon_max, lat_min, lat_max)
    ws_exp_1, u_exp_1, v_exp_1, time_exp_1, lat, lon = import_regcm5_wind(path, "NoTo-EUR", lon_min, lon_max, lat_min, lat_max)
    ws_exp_2, u_exp_2, v_exp_2, time_exp_2, _, _ = import_regcm5_wind(path, "WSM5-EUR", lon_min, lon_max, lat_min, lat_max)
    ws_exp_3, u_exp_3, v_exp_3, time_exp_3, _, _ = import_regcm5_wind(path, "WSM7-EUR", lon_min, lon_max, lat_min, lat_max)
    ws_exp_4, u_exp_4, v_exp_4, time_exp_4, _, _ = import_regcm5_wind(path, "WDM7-EUR", lon_min, lon_max, lat_min, lat_max)

    # 2D coordinates for vector quiver plotting
    lon_obs_2d, lat_obs_2d = np.meshgrid(lon_obs, lat_obs)

    # Define contour levels and colormap for wind speed (m/s)
    levels = np.arange(0, 31, 1)
    cmap = plt.cm.get_cmap("gist_ncar_r", len(levels) - 1)
    norm = matplotlib.colors.BoundaryNorm(levels, cmap.N)

    # Plot figure setup
    fig = plt.figure(figsize=(12, 6))
    font_size = 10
    skip = 12  # Thinning factor for wind arrows

    # Panel 1: ERA5
    ax1 = fig.add_subplot(2, 3, 1, projection=ccrs.PlateCarree())
    ax1.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax1.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax1.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle=':')
    cf1 = ax1.contourf(lon_obs_2d, lat_obs_2d, ws_obs, levels=levels, cmap=cmap, norm=norm, extend='max', transform=ccrs.PlateCarree())
    q1 = ax1.quiver(lon_obs_2d[::6, ::6], lat_obs_2d[::6, ::6], u_obs[::6, ::6], v_obs[::6, ::6], scale=300, color='black', transform=ccrs.PlateCarree())
    ax1.quiverkey(q1, 0.85, 1.03, 15, '15 m/s', labelpos='E', coordinates='axes', fontproperties={'size': 7})
    ax1.set_title(f"(a) {obs_name} ({time_obs:02d} UTC)", loc='left', fontsize=font_size, fontweight='bold')
    gl1 = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl1.top_labels = False
    gl1.right_labels = False

    # Panel 2: NoTo-EUR
    ax2 = fig.add_subplot(2, 3, 2, projection=ccrs.PlateCarree())
    ax2.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax2.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax2.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle=':')
    cf2 = ax2.contourf(lon, lat, ws_exp_1, levels=levels, cmap=cmap, norm=norm, extend='max', transform=ccrs.PlateCarree())
    q2 = ax2.quiver(lon[::skip, ::skip], lat[::skip, ::skip], u_exp_1[::skip, ::skip], v_exp_1[::skip, ::skip], scale=300, color='black', transform=ccrs.PlateCarree())
    ax2.quiverkey(q2, 0.85, 1.03, 15, '15 m/s', labelpos='E', coordinates='axes', fontproperties={'size': 7})
    ax2.set_title(f"(b) NoTo ({time_exp_1:02d} UTC)", loc='left', fontsize=font_size, fontweight='bold')
    gl2 = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl2.top_labels = False
    gl2.right_labels = False

    # Panel 3: WSM5-EUR
    ax3 = fig.add_subplot(2, 3, 3, projection=ccrs.PlateCarree())
    ax3.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax3.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax3.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle=':')
    cf3 = ax3.contourf(lon, lat, ws_exp_2, levels=levels, cmap=cmap, norm=norm, extend='max', transform=ccrs.PlateCarree())
    q3 = ax3.quiver(lon[::skip, ::skip], lat[::skip, ::skip], u_exp_2[::skip, ::skip], v_exp_2[::skip, ::skip], scale=300, color='black', transform=ccrs.PlateCarree())
    ax3.quiverkey(q3, 0.85, 1.03, 15, '15 m/s', labelpos='E', coordinates='axes', fontproperties={'size': 7})
    ax3.set_title(f"(c) WSM5 ({time_exp_2:02d} UTC)", loc='left', fontsize=font_size, fontweight='bold')
    gl3 = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl3.top_labels = False
    gl3.right_labels = False

    # Panel 4: WSM7-EUR
    ax4 = fig.add_subplot(2, 3, 4, projection=ccrs.PlateCarree())
    ax4.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax4.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax4.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle=':')
    cf4 = ax4.contourf(lon, lat, ws_exp_3, levels=levels, cmap=cmap, norm=norm, extend='max', transform=ccrs.PlateCarree())
    q4 = ax4.quiver(lon[::skip, ::skip], lat[::skip, ::skip], u_exp_3[::skip, ::skip], v_exp_3[::skip, ::skip], scale=300, color='black', transform=ccrs.PlateCarree())
    ax4.quiverkey(q4, 0.85, 1.03, 15, '15 m/s', labelpos='E', coordinates='axes', fontproperties={'size': 7})
    ax4.set_title(f"(d) WSM7 ({time_exp_3:02d} UTC)", loc='left', fontsize=font_size, fontweight='bold')
    gl4 = ax4.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl4.top_labels = False
    gl4.right_labels = False

    # Panel 5: WDM7-EUR
    ax5 = fig.add_subplot(2, 3, 5, projection=ccrs.PlateCarree())
    ax5.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax5.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax5.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle=':')
    cf5 = ax5.contourf(lon, lat, ws_exp_4, levels=levels, cmap=cmap, norm=norm, extend='max', transform=ccrs.PlateCarree())
    q5 = ax5.quiver(lon[::skip, ::skip], lat[::skip, ::skip], u_exp_4[::skip, ::skip], v_exp_4[::skip, ::skip], scale=300, color='black', transform=ccrs.PlateCarree())
    ax5.quiverkey(q5, 0.85, 1.03, 15, '15 m/s', labelpos='E', coordinates='axes', fontproperties={'size': 7})
    ax5.set_title(f"(e) WDM7 ({time_exp_4:02d} UTC)", loc='left', fontsize=font_size, fontweight='bold')
    gl5 = ax5.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl5.top_labels = False
    gl5.right_labels = False

    # Panel 6: Empty
    ax6 = fig.add_subplot(2, 3, 6)
    ax6.axis('off')

    # Colorbar 
    cbar_ax = fig.add_axes([0.15, 0.03, 0.7, 0.02])
    cbar = fig.colorbar(cf1, cax=cbar_ax, orientation='horizontal', label='10-m wind speed (m/s)')

    # Path out to save figure
    path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11/figs/paper'
    name_out = f'pyplt_maps_study_case_RegCM5_EUR-11_ws_{dt}.png'
    plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    main()
