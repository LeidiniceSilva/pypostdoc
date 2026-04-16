# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "May 14, 2025"
__description__ = "This script plot study area"

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Paths
path = "/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/era5/"
out_path = "/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/figs/"
os.makedirs(out_path, exist_ok=True)


def load_var(file, varname, time_sel=None):
    ds = xr.open_dataset(path + file)

    # corrigir longitude
    if ds.longitude.max() > 180:
        ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180)).sortby("longitude")

    # selecionar tempo (se fornecido)
    if time_sel is not None:
        ds = ds.sel(valid_time=time_sel)

    return ds[varname], ds["longitude"], ds["latitude"]


# Import vars
mslp, lon, lat = load_var("msl_ERA5_Otis_6hr_Oct2023.nc", "msl", "2023-10-24T06:00")
u10, _, _ = load_var("u10_ERA5_Otis_6hr_Oct2023.nc", "u10", "2023-10-24T06:00")
v10, _, _ = load_var("v10_ERA5_Otis_6hr_Oct2023.nc", "v10", "2023-10-24T06:00")
tp, _, _ = load_var("tp_ERA5_Otis_1hr_Oct2023.nc", "tp")

mslp = mslp / 100.0

wspd = np.sqrt(u10**2 + v10**2)

tp = tp.sel(valid_time=slice("2023-10-20", "2023-10-26"))
tp = tp.sum(dim="valid_time") * 1000.0  

# Plot figure
plt.figure(figsize=(15, 5))
font_size = 8

# Domain
lon_min, lon_max = -120, -80
lat_min, lat_max = 4, 28

def setup(ax):
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    
    # ticks
    ax.set_xticks(np.arange(lon_min, lon_max, 4), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(lat_min, lat_max, 4), crs=ccrs.PlateCarree())
    
    # format
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(font_size)
    
    # grid + features
    ax.grid(c='gray', ls='--')
    ax.coastlines(linewidth=0.75)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)

mslp_levels = np.arange(990, 1020, 2)   # hPa
wind_levels = np.arange(0, 20, 2)       # m/s
prec_levels = np.arange(0, 300, 30)     # mm/h

# MSLP
plt.subplot(1, 3, 1, projection=ccrs.PlateCarree())
ax = plt.gca()
setup(ax)
cf = plt.contourf(lon, lat, mslp, levels=mslp_levels, cmap="bwr", transform=ccrs.PlateCarree())
plt.title("(a)", loc='left', fontsize=font_size, fontweight='bold')
plt.colorbar(cf)

# Wind speed
plt.subplot(1, 3, 2, projection=ccrs.PlateCarree())
ax = plt.gca()
setup(ax)
cf = plt.contourf(lon, lat, wspd, levels=wind_levels, cmap="jet", transform=ccrs.PlateCarree())
plt.title("(b)", loc='left', fontsize=font_size, fontweight='bold')
plt.colorbar(cf)

# Total precipitation 
plt.subplot(1, 3, 3, projection=ccrs.PlateCarree())
ax = plt.gca()
setup(ax)
cf = plt.contourf(lon, lat, tp, levels=prec_levels, cmap="Blues", transform=ccrs.PlateCarree())
plt.title("(c)", loc='left', fontsize=font_size, fontweight='bold')
plt.colorbar(cf, ax=ax)

# Save figure
plt.tight_layout()
plt.savefig(out_path + "pyplt_Hurricane_Otis_synotic.png", dpi=400, bbox_inches="tight")
plt.show()
exit()


