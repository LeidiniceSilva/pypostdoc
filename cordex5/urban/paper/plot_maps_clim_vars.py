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

        # FIGURE: 6 ROWS x 5 COLUMNS
        fig, axes = plt.subplots(nrows=6,ncols=5,figsize=(17, 17),subplot_kw={"projection": ccrs.PlateCarree()})

        # MAP EXTENTS
        buenos_aires_extent = [-59.5, -57, -35.5, -33.5]
        sao_paulo_extent = [-48.0, -45.0, -24.5, -22.0]

        # ROW NAMES
        row_names = [f"{exp1} - ERA5",
            f"{exp2} - ERA5",
            f"{exp2} - {exp1}",
            f"{exp1} - ERA5",
            f"{exp2} - ERA5",
            f"{exp2} - {exp1}"]

        # REGION NAMES
        region_names = ["Buenos Aires", "São Paulo"]

        # FIGURE LETTERS
        letters = ["(a)", "(b)", "(c)", "(d)", "(e)",
            "(f)", "(g)", "(h)", "(i)", "(j)",
            "(k)", "(l)", "(m)", "(n)", "(o)",
            "(p)", "(q)", "(r)", "(s)", "(t)",
            "(u)", "(v)", "(w)", "(x)", "(y)",
            "(z)", "(aa)", "(ab)", "(ac)", "(ad)"]

        # LOOP OVER SEASONS
        for j, season in enumerate(seasons):
            print(f"  Processing season: {season}")

            # Model experiment 1
            exp1_var, exp1_lat, exp1_lon = import_regcm5_srf(path, param, domain, model, exp1, season, dt, res)

            # Model experiment 2
            exp2_var, exp2_lat, exp2_lon = import_regcm5_srf(path, param, domain, model, exp2, season, dt, res)

            # Observation
            obs_var, obs_lat, obs_lon = import_obs_srf(path, obs_param, domain, obs_name, season, dt, res)

            # Temporal mean
            if exp1_var.ndim == 3:
                exp1_var = np.nanmean(exp1_var, axis=0)

            if exp2_var.ndim == 3:
                exp2_var = np.nanmean(exp2_var, axis=0)

            if obs_var.ndim == 3:
                obs_var = np.nanmean(obs_var, axis=0)

            # Multiply ERA5 by -1
            obs_var = -obs_var

            # Bias and difference
            bias_exp1 = exp1_var - obs_var
            bias_exp2 = exp2_var - obs_var

            diff = exp2_var - exp1_var

            maps = [bias_exp1, bias_exp2, diff]

            # PLOT BUENOS AIRES
            for i in range(3):

                row = i
                ax = axes[row, j]
                ax.set_extent(buenos_aires_extent, crs=ccrs.PlateCarree())
                im = ax.pcolormesh(exp1_lon, exp1_lat, maps[i], transform=ccrs.PlateCarree(), cmap="bwr", vmin=-100, vmax=100, shading="auto")
                ax.coastlines(linewidth=0.8)
                ax.add_feature(cfeature.BORDERS, linewidth=0.5)

                gl = ax.gridlines(draw_labels=True, linewidth=0.5, linestyle="--", alpha=0.5)
                gl.top_labels = False
                gl.right_labels = False
                gl.bottom_labels = (row == 2)
                gl.left_labels = (j == 0)

                # Figure letter
                ax.text(0.02, 0.97, letters[row * 5 + j], transform=ax.transAxes, fontsize=11, fontweight="bold", va="top", ha="left")

            # PLOT SÃO PAULO
            for i in range(3):

                row = i + 3
                ax = axes[row, j]
                ax.set_extent(sao_paulo_extent, crs=ccrs.PlateCarree())
                im = ax.pcolormesh(exp1_lon, exp1_lat, maps[i], transform=ccrs.PlateCarree(), cmap="bwr", vmin=-100, vmax=100, shading="auto")
                ax.coastlines(linewidth=0.8)
                ax.add_feature(cfeature.BORDERS, linewidth=0.5)

                gl = ax.gridlines(draw_labels=True, linewidth=0.5, linestyle="--", alpha=0.5)
                gl.top_labels = False
                gl.right_labels = False
                gl.bottom_labels = (row == 5)
                gl.left_labels = (j == 0)

                # Figure letter
                ax.text(0.02, 0.97, letters[row * 5 + j], transform=ax.transAxes, fontsize=11, fontweight="bold", va="top", ha="left")

        # SEASON NAMES
        for j, season in enumerate(seasons):
            axes[0, j].text(0.50, 1.08, season, transform=axes[0, j].transAxes, ha="center", va="bottom", fontsize=13, fontweight="bold")

        # ROW NAMES
        for i, name in enumerate(row_names):
            axes[i, 0].text(-0.3, 0.50, name, transform=axes[i, 0].transAxes, rotation=90, va="center", ha="center", fontsize=10, fontweight="bold")

        # REGION LABELS
        axes[0, 0].text(-0.45, -0.7, region_names[0], transform=axes[0, 0].transAxes, rotation=90, va="center", ha="center", fontsize=13, fontweight="bold")
        axes[3, 0].text(-0.45, -0.7, region_names[1], transform=axes[3, 0].transAxes, rotation=90, va="center", ha="center", fontsize=13, fontweight="bold")

        # LAYOUT
        plt.subplots_adjust(left=0.12, right=0.90, bottom=0.08, top=0.92, wspace=0.08, hspace=0.20)

        # VERTICAL COLORBAR
        cax = fig.add_axes([0.91, 0.25, 0.018, 0.50])
        cbar = fig.colorbar(axes[0, 0].collections[0], cax=cax, orientation="vertical")
        cbar.set_label(f"{variable_name} ({unit})", fontsize=12)

        # SAVE FIGURE
        output_file = (output_dir / f"pyplt_maps_{param}_RegCM5_{domain}_{dt}.png")
        fig.savefig(output_file, dpi=400, bbox_inches="tight")
        plt.close(fig)

if __name__ == "__main__":
    main()
