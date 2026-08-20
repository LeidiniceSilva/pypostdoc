# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 28, 2026"
__description__ = "This script plots the annual cycle"


import pathlib
import glob
import numpy as np
import netCDF4
import matplotlib.pyplot as plt


def import_model(path, param, domain, model, exp, dt, res):

    arq = f"{path}/{param}_{domain}_{model}_{exp}_AC_{dt}_{res}_box.nc"
    files = glob.glob(arq)

    if len(files) == 0:
        raise FileNotFoundError(f"No file found: {arq}")

    data = netCDF4.Dataset(files[0])
    var = data.variables[param][:]

    if np.ma.isMaskedArray(var):
        var = var.filled(np.nan)

    data.close()

    return var


def import_obs(path, param, domain, obs, dt, res):

    arq = f"{path}/{param}_{domain}_{obs}_AC_{dt}_{res}_box.nc"
    files = glob.glob(arq)

    if len(files) == 0:
        raise FileNotFoundError(f"No file found: {arq}")

    data = netCDF4.Dataset(files[0])
    var = data.variables[param][:]

    if np.ma.isMaskedArray(var):
        var = var.filled(np.nan)

    data.close()

    return var


def main():

    path = pathlib.Path("/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/postproc/urban/paper")
    output_dir = pathlib.Path("/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/figs/urban/paper")
    output_dir.mkdir(parents=True, exist_ok=True)

    domain = "CSAM-3"
    dt = "2000-2009"
    res = "0.25"
    model = "RegCM5-ERA5"

    exp1 = "CTRL"
    exp2 = "URBAN"

    months = np.arange(1, 13)

    month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    # MODEL VARIABLE
    variables = {"pr": ("precip", "CPC"),
        "tasmax":  ("tmax",   "CPC"),
        "tasmin":  ("tmin",   "CPC"),
        "huss":    ("huss",   "ERA5"),
        "sfcWind": ("si10",   "ERA5"),
        "hfss":    ("msshf",  "ERA5"),
        "hfls":    ("mslhf",  "ERA5")}

    units = {"pr": r"Precipitation (mm day$^{-1}$)",
        "tasmax":  r"Maximum air temperature ($^\circ$C)",
        "tasmin":  r"Minimum air temperature ($^\circ$C)",
        "huss":    r"Specific humidity (kg kg$^{-1}$)",
        "sfcWind": r"Surface wind speed (m s$^{-1}$)",
        "hfss":    r"Sensible heat flux (W m$^{-2}$)",
        "hfls":    r"Latent heat flux (W m$^{-2}$)"}

    # FIGURE
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(16, 12))
    axes = axes.flatten()

    letters = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"]

    for n, param in enumerate(variables):
        print(f"Processing variable: {param}")

        ax = axes[n]
        obs_param, obs_name = variables[param]
        ctrl = import_model(path, param, domain, model, exp1, dt, res)
        urban = import_model(path, param, domain, model, exp2, dt, res)
        obs = import_obs(path, obs_param, domain, obs_name, dt, res)

        # SPATIAL AVERAGE
        obs = np.nanmean(obs, axis=(1, 2))
        ctrl = np.nanmean(ctrl, axis=(1, 2))
        urban = np.nanmean(urban, axis=(1, 2))

        # MULTIPLY ERA5 FLUXES BY -1
        if param in ["hfss", "hfls"]:
            obs = -obs

        # PLOT
        ax.plot(months, obs, color='black', linestyle='-', linewidth=2, label=obs_name)
        ax.plot(months, ctrl, color='blue', linestyle='-', linewidth=2, label="CTRL")
        ax.plot(months, urban, color='red', linestyle='-', linewidth=2, label="URBAN")
        ax.set_title(letters[n], loc='left', fontweight="bold")
        ax.set_xticks(months)
        ax.set_xticklabels(month_names)
        ax.set_xlim(1, 12)
        ax.set_ylabel(units[param], fontsize=10)
        ax.set_xlabel('Months', fontsize=10)
        ax.grid(linestyle="--", linewidth=0.5, alpha=0.5)
        ax.legend(fontsize=9, frameon=False)

    axes[7].axis("off")
    axes[8].axis("off")

    plt.subplots_adjust(left=0.08, right=0.97, bottom=0.08, top=0.93, wspace=0.25, hspace=0.30)

    # SAVE
    output_file = (output_dir / f"pyplt_annual_cycle_RegCM5_{domain}_{dt}.png")
    fig.savefig(output_file, dpi=400, bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    main()
