# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 28, 2026"
__description__ = "This script plot the climatology maps"

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


def import_regcm5_srf(path, param, model, exp, season, dt, res):

    arq = f"{path}/{param}_{model}_{exp}_{season}_{dt}_{res}.nc"
    files = glob.glob(arq)

    if len(files) == 0:
        raise FileNotFoundError(f"No file found: {arq}")

    data = netCDF4.Dataset(files[0])
    var = np.squeeze(data.variables[param][:])

    if np.ma.isMaskedArray(var):
        var = var.filled(np.nan)

    lat = np.squeeze(data.variables["lat"][:])
    lon = np.squeeze(data.variables["lon"][:])

    data.close()

    return var, lat, lon


def import_obs_srf(path, param, obs_name, season, dt, res):

    arq = f"{path}/{param}_{obs_name}_{season}_{dt}_{res}.nc"
    files = glob.glob(arq)

    if len(files) == 0:
        raise FileNotFoundError(f"No file found: {arq}")

    data = netCDF4.Dataset(files[0])
    var = np.squeeze(data.variables[param][:])

    if np.ma.isMaskedArray(var):
        var = var.filled(np.nan)

    lat = np.squeeze(data.variables["lat"][:])
    lon = np.squeeze(data.variables["lon"][:])

    data.close()

    return var, lat, lon


def main():

    path = pathlib.Path("/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11/postproc/paper")

    seasons = ["DJF", "MAM", "JJA", "SON", "ANN"]
    variables = ["pr", "tas", "clt"]

    dt = "2000-2009"
    res = "0.25"

    model = "RegCM5"
    exp1 = "NoTo-EUR"
    exp2 = "WSM5-EUR"
    exp3 = "WSM7-EUR"
    exp4 = "WDM7-EUR"

    # Generate letters list: ['(a)', '(b)', '(c)', ...]
    labels = [f"({letter})" for letter in string.ascii_lowercase]

    for param in variables:
        print(f"\nProcessing variable: {param}")

        # Define obs dataset and vars
        if param == "pr":
            obs_param = "precip"
            obs_name = "CPC"
            unit = "mm/day"
        elif param == "tas":
            obs_param = "t2m"
            obs_name = "ERA5"
            unit = "°C"
        elif param == "clt":
            obs_param = "tcc"
            obs_name = "ERA5"
            unit = "0-1"
        else:
            raise ValueError(f"Unknown variable: {param}")

        # Figure
        fig, axes = plt.subplots(nrows=5, ncols=5, figsize=(17, 15), subplot_kw={"projection": ccrs.PlateCarree()})

        # Store climatology and bias data
        obs_data = []
        exp_data = {
            exp1: [],
            exp2: [],
            exp3: [],
            exp4: []
        }

        # Read data
        for season in seasons:
            obs, lat, lon = import_obs_srf(path, obs_param, obs_name, season, dt, res)
            exp_1, lat, lon = import_regcm5_srf(path, param, model, exp1, season, dt, res)
            exp_2, lat, lon = import_regcm5_srf(path, param, model, exp2, season, dt, res)
            exp_3, lat, lon = import_regcm5_srf(path, param, model, exp3, season, dt, res)
            exp_4, lat, lon = import_regcm5_srf(path, param, model, exp4, season, dt, res)

            obs_data.append(obs)
            exp_data[exp1].append(exp_1 - obs)
            exp_data[exp2].append(exp_2 - obs)
            exp_data[exp3].append(exp_3 - obs)
            exp_data[exp4].append(exp_4 - obs)

        # Get spatial limits
        lon_min, lon_max = lon.min(), lon.max()
        lat_min, lat_max = lat.min(), lat.max()

        # Define contour levels
        if param == "pr":
            levels_obs = np.arange(0, 8.5, 0.5)
            levels_bias = np.arange(-4, 4.5, 0.5)
            color=['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff',\
                   '#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff',\
                   '#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff',\
                   '#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']
            cmap_obs = matplotlib.colors.ListedColormap(color)
            cmap_bias = "BrBG"

        elif param == "tas":
            levels_obs = np.arange(-4, 30, 2)
            levels_bias = np.arange(-5, 5.5, 0.5)
            cmap_obs = "jet"
            cmap_bias = "bwr"

        elif param == "clt":
            levels_obs = np.arange(0, 1.05, 0.05)
            levels_bias = np.arange(-0.5, 0.55, 0.05)
            cmap_obs = "Greys"
            cmap_bias = "RdGy"

        # Track subplot index for sequential labeling (a)-(y)
        panel_idx = 0

        # OBS climatology (Row 0)
        for j, season in enumerate(seasons):

            ax = axes[0, j]
            data = obs_data[j]
            im_obs = ax.contourf(lon, lat, data, levels=levels_obs, cmap=cmap_obs, extend="max", transform=ccrs.PlateCarree())
            ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
            ax.coastlines(linewidth=0.7)
            ax.add_feature(cfeature.BORDERS, linewidth=0.5)
            
            # Title only with item letter
            ax.set_title(labels[panel_idx], loc="left", fontsize=12, fontweight="bold")
            panel_idx += 1

            if j == 0:
                ax.text(-0.12, 0.5, f"{obs_name}", transform=ax.transAxes, rotation=90, va="center", ha="center", fontsize=12, fontweight="bold")

        # BIAS (Rows 1 to 4)
        experiments = [exp1, exp2, exp3, exp4]

        for i, exp in enumerate(experiments, start=1):
            for j, season in enumerate(seasons):

                ax = axes[i, j]
                data = exp_data[exp][j]
                im_bias = ax.contourf(lon, lat, data, levels=levels_bias, cmap=cmap_bias, extend="both", transform=ccrs.PlateCarree())
                ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
                ax.coastlines(linewidth=0.7)
                ax.add_feature(cfeature.BORDERS, linewidth=0.5)

                # Title with item letter
                ax.set_title(labels[panel_idx], loc="left", fontsize=12, fontweight="bold")
                panel_idx += 1

                if j == 0:
                    ax.text(-0.12, 0.5, exp, transform=ax.transAxes, rotation=90, va="center", ha="center", fontsize=12, fontweight="bold")
                if i == 4:
                    ax.text(0.50, -0.22, season, transform=ax.transAxes, ha="center", va="bottom", fontsize=12, fontweight="bold")

        # Layout adjustment for horizontal spacing
        plt.subplots_adjust(left=0.08, right=0.98, bottom=0.10, top=0.95, wspace=0.03)

        # Custom vertical positions for each row:
        # Distance between Row 0 (0.83) and Row 1 (0.65) is wider to leave room for the top colorbar.
        # Distance between Rows 1, 2, 3, and 4 is tight (step of ~0.115).
        row_y_positions = [0.83, 0.65, 0.535, 0.42, 0.305]

        for i in range(5):
            for j in range(5):
                pos = axes[i, j].get_position()
                axes[i, j].set_position([pos.x0, row_y_positions[i], pos.width, pos.height])

        # Colorbar for OBS (placed between Row 0 and Row 1)
        cbar_ax_obs = fig.add_axes([0.20, 0.8, 0.60, 0.012])
        cbar_obs = fig.colorbar(im_obs, cax=cbar_ax_obs, orientation="horizontal")
        cbar_obs.set_label(f"Mean ({unit})", fontsize=12)

        # Colorbar for BIAS (placed below Row 4)
        cbar_ax_bias = fig.add_axes([0.20, 0.265, 0.60, 0.012])
        cbar_bias = fig.colorbar(im_bias, cax=cbar_ax_bias, orientation="horizontal")
        cbar_bias.set_label(f"Bias ({unit})", fontsize=12)

        # Path out to save figure
        path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11/figs/paper'
        name_out = 'pyplt_maps_clim_{0}_RegCM5_EUR-11_{1}.png'.format(param, dt)
        plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
        plt.close(fig)

if __name__ == "__main__":
    main()
