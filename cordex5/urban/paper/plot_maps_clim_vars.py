# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 28, 2026"
__description__ = "This script plot the climatology maps"


import pathlib
import glob
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def import_regcm5_srf(path, param, domain, model, exp, season, dt, res):

    arq = f"{path}/{param}_{domain}_{model}_{exp}_{season}_{dt}_{res}.nc"
    files = glob.glob(arq)

    if len(files) == 0:
        raise FileNotFoundError(f"No file found: {arq}")

    data = netCDF4.Dataset(files[0])
    var = data.variables[param][:]

    if np.ma.isMaskedArray(var):
        var = var.filled(np.nan)

    lat = data.variables["lat"][:]
    lon = data.variables["lon"][:]

    data.close()

    return var, lat, lon


def import_obs_srf(path, param, domain, obs_name, season, dt, res):

    arq = f"{path}/{param}_{domain}_{obs_name}_{season}_{dt}_{res}.nc"
    files = glob.glob(arq)

    if len(files) == 0:
        raise FileNotFoundError(f"No file found: {arq}")

    data = netCDF4.Dataset(files[0])
    var = data.variables[param][:]

    if np.ma.isMaskedArray(var):
        var = var.filled(np.nan)

    lat = data.variables["lat"][:]
    lon = data.variables["lon"][:]

    data.close()

    return var, lat, lon


def main():

    path = pathlib.Path("/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/postproc/urban/paper")
    output_dir = pathlib.Path("/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/figs/urban/paper")
    output_dir.mkdir(parents=True, exist_ok=True)

    seasons = ["DJF", "MAM", "JJA", "SON", "ANN"]
    variables = ["hfls", "hfss"]

    dt = "2000-2009"
    domain = "CSAM-3"
    res = "0.25"

    model = "RegCM5-ERA5"

    exp1 = "CTRL"
    exp2 = "URBAN"

    obs_name = "ERA5"

    for param in variables:

        print(f"\nProcessing variable: {param}")

        if param == "hfls":
            obs_param = "mslhf"
            variable_name = "Latent heat flux"
        elif param == "hfss":
            obs_param = "msshf"
            variable_name = "Sensible heat flux"
        else:
            raise ValueError(f"Variable not defined: {param}")

        unit = "W m$^{-2}$"

        # Figure: 5 rows x 3 columns
        fig, axes = plt.subplots(
            nrows=5,
            ncols=3,
            figsize=(15, 15),
            subplot_kw={
                "projection": ccrs.PlateCarree()
            }
        )

        letters = list("abcdefghijklmno")

        # Column names
        column_names = [f"{exp1} - ERA5", f"{exp2} - ERA5", f"{exp2} - {exp1}"]

        # Loop over seasons
        for i, season in enumerate(seasons):

            print(f"  Processing season: {season}")

            # Model experiment 1
            exp1_var, exp1_lat, exp1_lon = import_regcm5_srf(
                path,
                param,
                domain,
                model,
                exp1,
                season,
                dt,
                res
            )

            # Model experiment 2
            exp2_var, exp2_lat, exp2_lon = import_regcm5_srf(
                path,
                param,
                domain,
                model,
                exp2,
                season,
                dt,
                res
            )

            # Observation
            obs_var, obs_lat, obs_lon = import_obs_srf(
                path,
                obs_param,
                domain,
                obs_name,
                season,
                dt,
                res
            )

            # Temporal mean
            if exp1_var.ndim == 3:
                exp1_var = np.nanmean(exp1_var, axis=0)
            if exp2_var.ndim == 3:
                exp2_var = np.nanmean(exp2_var, axis=0)
            if obs_var.ndim == 3:
                obs_var = np.nanmean(obs_var, axis=0)

            # Multiply ERA5 by -1
            obs_var = -obs_var

            # Bias
            bias_exp1 = (exp1_var - obs_var)
            bias_exp2 = (exp2_var - obs_var)

            # Difference between experiments
            diff = (exp2_var - exp1_var)

            maps = [bias_exp1, bias_exp2, diff]

            # Plot each column
            for j in range(3):

                ax = axes[i, j]

                im = ax.pcolormesh(
                    exp1_lon,
                    exp1_lat,
                    maps[j],
                    transform=ccrs.PlateCarree(),
                    cmap="RdBu_r",
                    vmin=-100,
                    vmax=100,
                    shading="auto"
                )

                # Coastlines
                ax.coastlines(
                    linewidth=0.8
                )

                ax.add_feature(
                    cfeature.BORDERS,
                    linewidth=0.5
                )

                # Latitude / Longitude
                gl = ax.gridlines(
                    draw_labels=True,
                    linewidth=0.5,
                    linestyle="--",
                    alpha=0.5
                )

                gl.top_labels = False
                gl.right_labels = False

                gl.bottom_labels = (
                    i == 4
                )

                gl.left_labels = (
                    j == 0
                )

                # Figure letter
                ax.set_title(
                    f"({letters[i * 3 + j]})",
                    fontsize=11,
                    fontweight="bold",
                    loc="left"
                )

            # Season on the left side
            axes[i, 0].text(
                -0.20,
                0.50,
                season,
                transform=axes[i, 0].transAxes,
                rotation=90,
                va="center",
                ha="center",
                fontsize=12,
                fontweight="bold"
            )

        # Column names below the maps
        for j, name in enumerate(column_names):

            axes[4, j].text(
                0.50,
                -0.20,
                name,
                transform=axes[4, j].transAxes,
                ha="center",
                va="top",
                fontsize=12,
                fontweight="bold"
            )

        # Layout
        plt.tight_layout(rect=[0.03, 0.04, 1, 0.97])

        # Colorbars
        colorbar_positions = [
            [0.11, 0.005, 0.25, 0.018],
            [0.41, 0.005, 0.25, 0.018],
            [0.72, 0.005, 0.25, 0.018]
        ]

        colorbar_labels = [
            f"{variable_name} ({unit})",
            f"{variable_name} ({unit})",
            f"{variable_name} ({unit})"
        ]

        for j in range(3):
            cax = fig.add_axes(colorbar_positions[j])
            cbar = fig.colorbar(axes[0, j].collections[0], cax=cax, orientation="horizontal")
            cbar.set_label(colorbar_labels[j], fontsize=10)
            cbar.set_ticks([-100, -50, 0, 50, 100])

        # Save figure
        output_file = (output_dir / f"pyplt_maps_{param}_RegCM5_{domain}_{dt}.png")
        fig.savefig(output_file, dpi=400, bbox_inches="tight")
        plt.close(fig)

        print(f"Figure saved: {output_file}")

if __name__ == "__main__":
    main()

