"""
Plot Tb (gray_r) with precipitation contourf (jet) overlaid — 2×2 subplots.
Panels: GPM | ERA5 | RegCM5-CSAM-3 | RegCM5-SAM-12
Date: 2007-12-05, timestep index 10 of that day.
No title. Dataset name + panel letter inside each subplot.
Subplots glued together, no ticks, no axis labels.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import xarray as xr

OUTFIG = (
    "/leonardo/home/userexternal/mdasilva/leonardo_work/"
    "TRACK-CYCLONE/CPRCM-TF/figs/"
    "pyplt_maps_event_CSAM-3_20071205.png"
)

DATADIR    = "/leonardo/home/userexternal/mdasilva/leonardo_work/TRACK-CYCLONE/CPRCM-TF/data/Dataset/"
TARGET_DAY = "2007-12-05"
TIME_STEP  = 8          # 0-based index within TARGET_DAY

DOMAIN = [-66, -38, -36, -14]   # [lon_min, lon_max, lat_min, lat_max]

TB_LEVELS   = np.arange(200, 320, 5)                         # K
PREC_LEVELS = [0.1, 0.5, 1, 2, 4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50]      # mm h⁻¹
TB_CMAP     = "gray_r"
PREC_CMAP   = "jet"

DATASETS = {
    "GPM": {
        "tb_file" : DATADIR + "Tb_CSAM-3_Kelvin_GPM_20071205.nc",
        "pr_file" : DATADIR + "PR_CSAM-3_mmhr_GPM_20071205.nc",
        "tb_var"  : "Tb",
        "pr_var"  : "PR",
    },
    "ERA5": {
        "tb_file" : DATADIR + "Tb_CSAM-3_Kelvin_ERA5_20071205.nc",
        "pr_file" : DATADIR + "tp_CSAM-3_mmhr_ERA5_20071205.nc",
        "tb_var"  : "Tb",
        "pr_var"  : "tp",
    },
    "RegCM5 CSAM-3": {
        "tb_file" : DATADIR + "Tb_CSAM-3_Kelvin_RegCM5_20071205.nc",
        "pr_file" : DATADIR + "pr_CSAM-3_mmhr_RegCM5_20071205.nc",
        "tb_var"  : "Tb",
        "pr_var"  : "pr",
    },
    "RegCM5 SAM-12": {
        "tb_file" : DATADIR + "Tb_SAM-12_Kelvin_RegCM5_20071205.nc",
        "pr_file" : DATADIR + "pr_SAM-12_mmhr_RegCM5_20071205.nc",
        "tb_var"  : "Tb",
        "pr_var"  : "pr",
    },
}


def find_coord(ds, candidates):
    for c in candidates:
        if c in ds.coords or c in ds.data_vars:
            return c
    raise KeyError(
        f"None of {candidates} found. Available: {list(ds.coords) + list(ds.data_vars)}"
    )


def get_latlon(ds):
    lat_name = find_coord(ds, ["lat", "latitude", "xlat", "rlat", "LAT", "Lat"])
    lon_name = find_coord(ds, ["lon", "longitude", "xlon", "rlon", "LON", "Lon"])
    lat = ds[lat_name].values
    lon = ds[lon_name].values
    if lat.ndim == 1 and lon.ndim == 1:
        lon, lat = np.meshgrid(lon, lat)
    return lat, lon


def select_step(ds, var_name, target_day, step_index):
    tname = find_coord(ds, ["time", "valid_time", "TIME", "Time"])
    times = pd.DatetimeIndex(ds[tname].values)
    mask  = times.date == pd.Timestamp(target_day).date()
    idxs  = np.where(mask)[0]
    if len(idxs) == 0:
        raise ValueError(f"No steps found for {target_day}.")
    if step_index >= len(idxs):
        raise IndexError(f"step_index={step_index} but only {len(idxs)} steps on {target_day}.")
    da = ds[var_name].isel({tname: idxs[step_index]})
    for dim in list(da.dims):
        if dim not in ("lat","lon","latitude","longitude",
                       "xlat","xlon","rlat","rlon","x","y",
                       "south_north","west_east","ncells"):
            da = da.isel({dim: 0})
    return da


def load(file_path, var_name):
    ds = xr.open_dataset(file_path, decode_times=True)
    if var_name not in ds.data_vars:
        candidates = list(ds.data_vars)
        print(f"  [warn] '{var_name}' not in {file_path}; trying '{candidates[0]}'")
        var_name = candidates[0]
    da  = select_step(ds, var_name, TARGET_DAY, TIME_STEP)
    lat, lon = get_latlon(ds)
    return lat, lon, da.values


def configure_subplot(ax):
    ax.set_extent(DOMAIN, crs=ccrs.PlateCarree())
    ax.add_feature(cfeat.BORDERS, linewidth=0.75, color="black")
    ax.coastlines(linewidth=0.75, color="black")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)


def main():

    tb_norm   = mcolors.BoundaryNorm(TB_LEVELS,   plt.colormaps[TB_CMAP].N)

    labels = list(DATASETS.keys())
    nrows, ncols = 2, 2

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(12, 10),
        subplot_kw={"projection": ccrs.PlateCarree()},
        gridspec_kw={"hspace": 0.01, "wspace": 0.01},
    )
    axes_flat = axes.flatten()

    tb_pcm   = None
    prec_cf  = None          # contourf handle for colorbar

    panel_letters = ["Satellite estimation (10 km)", "Reanalysis (25 km)", "Convection-permitting regional climate model (3 km)", "Regional climate model (12 km)"]

    for i, (label, cfg) in enumerate(DATASETS.items()):
        ax = axes_flat[i]
        configure_subplot(ax)

        try:
            lat_tb, lon_tb, tb = load(cfg["tb_file"], cfg["tb_var"])
            tb_pcm = ax.pcolormesh(
                lon_tb, lat_tb, tb,
                cmap=TB_CMAP, norm=tb_norm,
                transform=ccrs.PlateCarree(),
                shading="auto", zorder=1,
            )
        except Exception as e:
            print(f"[{label}] Tb error: {e}")

        try:
            lat_pr, lon_pr, pr = load(cfg["pr_file"], cfg["pr_var"])

            # contourf needs regular or 2-D grids; mask values below first level
            pr_plot = np.where(pr >= PREC_LEVELS[0], pr, np.nan)

            prec_cf = ax.contourf(
                lon_pr, lat_pr, pr_plot,
                levels=PREC_LEVELS,
                cmap=PREC_CMAP,
                extend="max",
                transform=ccrs.PlateCarree(),
                zorder=2,
                alpha=0.85,
            )
        except Exception as e:
            print(f"[{label}] Precip error: {e}")

        ax.text(
            0.97, 0.97,
            f"{panel_letters[i]}",
            transform=ax.transAxes,
            fontsize=10,
            va="top", ha="right",
            color="white",
            bbox=dict(facecolor="black", alpha=0.75, edgecolor="none", pad=3),
            zorder=6,
        )

    if tb_pcm is not None:
        cax_tb = fig.add_axes([0.91, 0.53, 0.015, 0.33])
        cb_tb  = fig.colorbar(tb_pcm, cax=cax_tb,
                               extend="both")
        cb_tb.set_label("Brightness Temperature (K)", fontsize=8)
        cb_tb.ax.tick_params(labelsize=8)

    if prec_cf is not None:
        cax_pr = fig.add_axes([0.91, 0.14, 0.015, 0.33])
        cb_pr  = fig.colorbar(prec_cf, cax=cax_pr)
        cb_pr.set_label("Precipitation (mm h⁻¹)", fontsize=8)
        cb_pr.set_ticks(PREC_LEVELS)
        cb_pr.ax.tick_params(labelsize=8)

    import os
    os.makedirs(os.path.dirname(OUTFIG), exist_ok=True)
    plt.savefig(OUTFIG, dpi=400, bbox_inches="tight")
    print(f"Saved → {OUTFIG}")
    plt.show()


if __name__ == "__main__":
    main()
