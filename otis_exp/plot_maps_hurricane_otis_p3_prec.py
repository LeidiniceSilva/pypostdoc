# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 17, 2026"
__description__ = "This script plot Otis hurricane precipitation"

import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import warnings
warnings.filterwarnings('ignore')

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--domain", default="large")
parser.add_argument("--experiment", default="exps_v1")
parser.add_argument("--date", default="2023102000")
args = parser.parse_args()

domain = args.domain
experiment = args.experiment
date = args.date

MAP_EXTENT = [-110, -90, 5, 25]

EXPERIMENTS = ["ctrl", "holt_r2", "holt_r3", "uw_r2", "uw_r3"]

START_DATE = "2023-10-24 00:00:00"
END_DATE   = "2023-10-25 23:00:00"

TS_START_DATE = "2023-10-24 15:00:00"
TS_END_DATE   = "2023-10-25 13:00:00"

ACAPULCO_LON = -99.91
ACAPULCO_LAT = 16.85

# Paths
DATA_DIR_MODELS = (f"/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/exps/{experiment}/domain_{domain}_regridded/")
DATA_DIR_CMORPH = ("/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/cmorph")
OUTPUT_PATH = (f"/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/figs/{experiment}/")
os.makedirs(OUTPUT_PATH, exist_ok=True)


# Functions
def fix_longitude(ds):

    if np.any(ds.lon > 180):
        ds['lon'] = (ds.lon + 180) % 360 - 180
        ds = ds.sortby('lon')

    return ds


def load_cmorph():

    nc_file = os.path.join(
        DATA_DIR_CMORPH,
        "CMORPH_cut_domain_20251022_20251025.nc"
    )

    print(f'Opening {nc_file}')

    ds = xr.open_dataset(nc_file)
    ds = fix_longitude(ds)
    ds = ds.sel(time=slice(START_DATE, END_DATE))

    pr_acc = ds['cmorph'].sum(dim='time')

    pr_ts = ds['cmorph'].sel(
        lon=ACAPULCO_LON,
        lat=ACAPULCO_LAT,
        method='nearest'
    )

    return pr_acc, pr_ts


def load_model(exp):

    file_name = os.path.join(
        DATA_DIR_MODELS,
        exp,
        f'pr_{exp}_{date}.nc'
    )

    print(f'Opening {file_name}')

    ds = xr.open_dataset(file_name)
    ds = fix_longitude(ds)
    ds = ds.sel(time=slice(START_DATE, END_DATE))

    pr_acc = ds['pr'].sum(dim='time') * 3600.

    pr_ts = ds['pr'].sel(
        lon=ACAPULCO_LON,
        lat=ACAPULCO_LAT,
        method='nearest'
    ) * 3600.

    return pr_acc, pr_ts


# Load data
all_precip = {}
all_ts = {}

print("Loading CMORPH...")
all_precip["CMORPH"], cmorph_ts = load_cmorph()

# Station
station = pd.read_csv("/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/ws/station_acapulco_timeseries.csv")
station['time'] = pd.to_datetime(station['time'])
station = station.set_index('time')
station = station.loc[TS_START_DATE:TS_END_DATE]
station = station['pr']
station_value = station.sum()
print(f"Loaded station: {len(station)} hours")

for exp in EXPERIMENTS:
    print(f"Loading {exp}...")
    pr_acc, pr_ts = load_model(exp)
    all_precip[exp.upper()] = pr_acc
    all_ts[exp.upper()] = pr_ts

# Plot
print("Plotting")
fig, axes = plt.subplots(2, 3, figsize=(14, 9), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()

levels = np.arange(0, 305, 5)

base = plt.cm.terrain_r(np.linspace(0, 1, len(levels)-1))
base[0] = [1, 1, 1, 1]
cmap = mcolors.ListedColormap(base)
norm = mcolors.BoundaryNorm(levels, ncolors=cmap.N + 1, extend='max')

titles = ['CMORPH', 'CTRL', 'HOLT_R2', 'HOLT_R3', 'UW_R2', 'UW_R3']
labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']

for i, name in enumerate(titles):

    ax = axes[i]
    pr = all_precip[name]
    cf = ax.contourf(pr.lon, pr.lat, pr, levels=levels, cmap=cmap, extend='max', transform=ccrs.PlateCarree())

    if name == 'CMORPH':
        ax.scatter(ACAPULCO_LON, ACAPULCO_LAT,s=60, c=[station_value], cmap=cmap, norm=norm, edgecolor='black', transform=ccrs.PlateCarree(), zorder=10)

    ax.set_extent(MAP_EXTENT)
    ax.coastlines(linewidth=0.75)
    ax.add_feature(cfeature.BORDERS, linewidth=0.75)

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, linestyle='--', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False

    # Inset time series
    axins = inset_axes(ax, width="40%", height="30%", loc='upper left', bbox_transform=ax.transAxes)
    axins.patch.set_facecolor('lightgray')
    axins.patch.set_alpha(1.0)

    if name == 'CMORPH':
        axins.plot(station.index, station.values, color='black', linewidth=1.5, label='Station')
        axins.plot(cmorph_ts.time, cmorph_ts, color='blue', linewidth=1.5, label='CMORPH')
        axins.legend(fontsize=6, frameon=False)
    else:
        axins.plot(all_ts[name].time, all_ts[name], color='red', linewidth=1.5)

    axins.set_ylim(0, 35)
    axins.yaxis.set_label_position("right")
    axins.tick_params(axis='both', labelsize=6)
    axins.tick_params(axis='x', rotation=90)
    axins.tick_params(axis='y', left=False, labelleft=False, right=True, labelright=True)
    axins.set_xticks(axins.get_xticks())
    axins.set_ylabel('mm h$^{-1}$', fontsize=6)

    ax.set_title(f"{labels[i]} {name}", loc='left', fontsize=10, fontweight='bold')

# Colorbar
cax = fig.add_axes([0.999, 0.15, 0.02, 0.7])
cbar = fig.colorbar(cf, cax=cax)
cbar.set_label("Accumulated precipitation 24-25Oct2023 (mm)", fontsize=10, fontweight='bold')

# Save figure
print("Save figure")
plt.tight_layout()
outfile = os.path.join(OUTPUT_PATH, f"pyplt_Hurricane_Otis_precip_accum_{domain}.png")
plt.savefig(outfile, dpi=400, bbox_inches='tight')
plt.show()

print(f"Figure saved: {outfile}")
