# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 16, 2023"
__description__ = "This script plot urban mask"

import os
import numpy as np
import netCDF4 as nc
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker

from matplotlib.colors import ListedColormap, BoundaryNorm

var='sftlaf'


def load_data(fname):

    with nc.Dataset(fname, "r") as f:
        var_data = f.variables[var][:]
        lat = f.variables["lat"][:]
        lon = f.variables["lon"][:]

    return var_data, lat, lon


def load_data_urb(fname):

    with nc.Dataset(fname, "r") as f:
        var_data = f.variables[var][:]
        lat_urb = f.variables["lat"][:]
        lon_urb = f.variables["lon"][:]

    return var_data, lat, lon


outp = "/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/figs/urb"

var_data, lat, lon = load_data("/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/ERA5/ERA5-CSAM/CORDEX-CMIP6/DD/CSAM-3/ICTP/ERA5/evaluation/r1i1p1f1/RegCM5-0/v1-r1/fx/{0}/{0}_CSAM-3_ERA5_evaluation_r1i1p1f1_ICTP_RegCM5-0_v1-r1_fx.nc".format(var))

var_data_urb, lat_urb, lon_urb = load_data_urb("/leonardo_scratch/large/userexternal/ggiulian/urban/output/CORDEX-CMIP6/DD/CSAM-3/ICTP/ERA5/evaluation/r1i1p1f1/RegCM5-0/v1-r1/fx/{0}/{0}_CSAM-3_ERA5_evaluation_r1i1p1f1_ICTP_RegCM5-0_v1-r1_fx.nc".format(var))

# Plot both subplots in one figure
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
font_size = 10

# Urban Fraction Plot
ax = axes[0]
ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()])
ax.add_feature(cfeature.COASTLINE, edgecolor='black')
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.set_title("(a) RegCM5 {0}".format(var), loc='left', fontsize=font_size, fontweight='bold')
ax.text(-38, -38, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
im = ax.contourf(lon, lat, var_data, cmap="gist_ncar_r", transform=ccrs.PlateCarree())
cbar = plt.colorbar(im, ax=ax, orientation='horizontal')

# Urban Fraction Plot
ax = axes[1]
ax.set_extent([lon_urb.min(), lon_urb.max(), lat_urb.min(), lat_urb.max()])
ax.add_feature(cfeature.COASTLINE, edgecolor='black')
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.set_title("(b) RegCM5 URB {0}".format(var), loc='left', fontsize=font_size, fontweight='bold')
ax.text(-38, -38, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
im = ax.contourf(lon, lat, var_data_urb, cmap="gist_ncar_r", transform=ccrs.PlateCarree())
cbar = plt.colorbar(im, ax=ax, orientation='horizontal')

# Save and close figure
plt.tight_layout()
os.makedirs(outp, exist_ok=True)
plt.savefig(os.path.join(outp, "pyplt_maps_{0}_RegCM5_CSAM-3_2000-2009.png".format(var)), dpi=400, bbox_inches='tight')
plt.show()
exit()

