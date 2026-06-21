# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 17, 2026"
__description__ = "This script plot Otis hurricane precipitation"

import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import warnings
warnings.filterwarnings('ignore')

# Domain type
DOMAIN = "small"
EXPS_ = "exps_v2"
DATE = "2023101900"
MAP_EXTENT = [-110, -90, 5, 25]

# Experiment names
EXPERIMENTS = ["ctrl", "holt_r2", "holt_r3", "uw_r2", "uw_r3"]

START_DATE = "2023-10-24 15:00:00"
END_DATE   = "2023-10-25 13:00:00"

# Paths
DATA_DIR_MODELS = (f"/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/exps/{EXPS_}/domain_{DOMAIN}_regridded/")
DATA_DIR_CMORPH = ("/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/cmorph")
OUTPUT_PATH = (f"/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/figs/{EXPS_}/")
os.makedirs(OUTPUT_PATH, exist_ok=True)


# Function load dataset
def load_cmorph_precip():

    nc_files = os.path.join(DATA_DIR_CMORPH, "CMORPH_cut_domain_20251022_20251025.nc")
    ds = xr.open_dataset(nc_files)
    print(f'Opening {nc_files}')

    # Fix longitude
    if np.any(ds.lon > 180):
        ds['lon'] = (ds.lon + 180) % 360 - 180
        ds = ds.sortby('lon')

    ds = ds.sel(time=slice(START_DATE, END_DATE))

    # CMORPH total precipitation mm/h
    pr = ds['cmorph'].sum(dim='time') 

    return pr


# Function load dataset
def load_model_precip(exp):
    """Load accumulated precipitation for one experiment"""

    file_name = os.path.join(DATA_DIR_MODELS, exp, f'pr_{exp}_{DATE}.nc')
    ds = xr.open_dataset(file_name)
    print(f'Opening {file_name}')

    # Fix longitude
    if np.any(ds.lon > 180):
        ds['lon'] = (ds.lon + 180) % 360 - 180
        ds = ds.sortby('lon')

    ds = ds.sel(time=slice(START_DATE, END_DATE))

    # Sum hourly precipitation over the 48 hours
    pr_acc = ds['pr'].sum(dim='time') * 3600.

    return pr_acc


# Load data
all_precip = {}

print("Loading CMOPRH...")
all_precip["CMORPH"] = load_cmorph_precip()

# Load station data: acapulco_lon, acapulco_lat = -99.91, 16.85
station = pd.read_csv("/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/extracted_timeseries/station_acapulco_timeseries.csv")
station['time'] = pd.to_datetime(station['time'])
station = station.set_index('time')
station = station['pr']
station_value = station.sum()
print(f"  Loaded station: {len(station)} hours")

for exp in EXPERIMENTS:

    print(f"Loading {exp}...")
    all_precip[exp.upper()] = load_model_precip(exp)

# Plot
print("Ploting")
fig, axes = plt.subplots(2, 3, figsize=(14, 10), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()

levels = np.arange(0, 310, 10)

# Create colormap starting with white
base = plt.cm.terrain_r(np.linspace(0, 1, len(levels)-1))
base[0] = [1, 1, 1, 1]   
cmap = mcolors.ListedColormap(base)

norm = mcolors.BoundaryNorm(levels, ncolors=cmap.N, extend='max')

titles = ['CMORPH', 'CTRL', 'HOLT_R2', 'HOLT_R3', 'UW_R2', 'UW_R3']
labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']

for i, name in enumerate(titles):

    ax = axes[i]

    pr = all_precip[name]
    cf = ax.contourf(pr.lon, pr.lat, pr, levels=levels, cmap=cmap, extend='max', transform=ccrs.PlateCarree())

    if name == 'CMORPH':
        ax.scatter(-99.91, 16.85, s=60, c=[station_value], cmap=cmap, norm=norm, edgecolor='black', transform=ccrs.PlateCarree(), zorder=10)

    ax.set_extent(MAP_EXTENT, crs=ccrs.PlateCarree())
    ax.coastlines(linewidth=0.75)
    ax.add_feature(cfeature.BORDERS, linewidth=0.75)

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, linestyle='--', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False

    ax.set_title(f"{labels[i]} {name}", fontsize=10, fontweight='bold')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

# Colorbar
cax = fig.add_axes([0.999, 0.15, 0.02, 0.7])
cbar = fig.colorbar(cf, cax=cax)
cbar.set_label("Accumulated precipitation 24-25Oct2023 (mm)", fontsize=10, fontweight='bold')

# Save figure
print("Save figure")
plt.tight_layout()
outfile = os.path.join(OUTPUT_PATH, f"pyplt_Hurricane_Otis_precip_accum_{DOMAIN}.png")
plt.savefig(outfile, dpi=400, bbox_inches='tight')
plt.show()

print(f"Figure saved: {outfile}")
