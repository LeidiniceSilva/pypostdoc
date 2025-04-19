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

# Load environment variables
fname = "/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/ERA5/icbc/CSAM-3_CLM45_surface.nc"
snum = os.getenv("snum")
outp = "/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/figs" 

# Read dataset
f = nc.Dataset(fname)
urb_2d = f.variables["urb_2d"][:]
urb_frac = np.sum(urb_2d, axis=0)  # Summing over time (assuming time is first dimension)

lat = f.variables["xlat"][:]
lon = f.variables["xlon"][:]

# Convert longitude if necessary
if snum in ["Australasia", "EastAsia"]:
    lon = np.where(lon < 0, lon + 360, lon)

# Check resolution
rr = abs(lat[20, 20] - lat[20, 21])
rr_03_diff = abs(rr - 0.03)
rr_12_diff = abs(rr - 0.11)
rr_25_diff = abs(rr - 0.22)

# Determine urban mask threshold
if rr_03_diff < rr_12_diff:
    print("3 km resolution")
    this_frac = 70
elif rr_12_diff < rr_25_diff:
    print("12 km resolution")
    this_frac = 40
else:
    print("25 km resolution")
    this_frac = 10

# Create urban mask
urb_mask = np.where(urb_frac > this_frac, 1, 0)

# Plot both subplots in one figure
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree()})
font_size = 10

# Urban Fraction Plot
ax = axes[0]
ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()])
ax.add_feature(cfeature.COASTLINE, edgecolor='black')
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.set_title("(a) Urban Fraction", loc='left', fontsize=font_size, fontweight='bold')
ax.text(-38, -38, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
im = ax.contourf(lon, lat, urb_frac, levels=np.linspace(0, 100, 11), cmap="Blues", transform=ccrs.PlateCarree())
cbar = plt.colorbar(im, ax=ax, orientation='horizontal')
cbar.set_label('%')

gl = ax.gridlines(draw_labels=True, linestyle="-")
gl.top_labels = False
gl.right_labels = False
gl.xformatter = cticker.LongitudeFormatter()
gl.yformatter = cticker.LatitudeFormatter()

# Urban Mask Plot
ax = axes[1]
ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()])
ax.add_feature(cfeature.COASTLINE, edgecolor='black')
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.set_title(f"(b) Urban Mask (frac > {this_frac})", loc='left', fontsize=font_size, fontweight='bold')
ax.text(-38, -38, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')

colors = ["white", "red"]
cmap = ListedColormap(colors)
norm = BoundaryNorm([0, 0.7, 1], cmap.N)

im = ax.contourf(lon, lat, urb_mask, levels=[0, 0.7, 1], cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cbar = plt.colorbar(im, ax=ax, orientation='horizontal', ticks=[0, 0.7, 1])

gl = ax.gridlines(draw_labels=True, linestyle="-")
gl.top_labels = False
gl.right_labels = False
gl.xformatter = cticker.LongitudeFormatter()
gl.yformatter = cticker.LatitudeFormatter()

# Save and close figure
plt.tight_layout()
os.makedirs(outp, exist_ok=True)
plt.savefig(os.path.join(outp, "pyplt_maps_urban_mask_RegCM5_CSAM-3_2000-2009.png"), dpi=400, bbox_inches='tight')
plt.show()
exit()

