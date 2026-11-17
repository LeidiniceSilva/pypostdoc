# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "April 14, 2026"
__description__ = "This script plots ALPX-3 topography"

import os
import glob
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

from netCDF4 import Dataset as nc

# Paths 
dirnc    = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/dataset/domains'
number_i = 29

# Region of interest 
LON_MIN, LON_MAX = -15, 25
LAT_MIN, LAT_MAX = 35, 55

# Reader 
def import_rcm(file, number):

    RCMf = nc(file, mode='r')

    if 'xlat' in RCMf.variables:
        lat = RCMf.variables['xlat'][:]
        lon = RCMf.variables['xlon'][:]
    else:
        lat = RCMf.variables['lat'][:]
        lon = RCMf.variables['lon'][:]

    if 'orog' in RCMf.variables:
        topo = RCMf.variables['orog'][:]
    elif 'topo' in RCMf.variables:
        topo = RCMf.variables['topo'][:]
    elif 'z' in RCMf.variables:
        topo = RCMf.variables['z'][:]
    else:
        topo = np.zeros_like(lat)

    RCMf.close()

    lat  = np.squeeze(lat)
    lon  = np.squeeze(lon)
    topo = np.squeeze(topo)

    if lat.ndim == 1 and lon.ndim == 1:
        lon, lat = np.meshgrid(lon, lat)

    topo = np.where(topo < 0, 0, topo)

    ny, nx = topo.shape
    border_mask = np.full((ny, nx), np.nan)
    border_mask[:number, :]  = 1
    border_mask[-number:, :] = 1
    border_mask[:, :number]  = 1
    border_mask[:, -number:] = 1

    return lat, lon, topo, border_mask


# Load ALPX-3 
files = sorted(glob.glob(os.path.join(dirnc, 'orog_*.nc')))
domains_data = {}

for file in files:
    domain = os.path.basename(file).split('_')[1]
    lat, lon, topo, border_mask = import_rcm(file, number_i)
    domains_data[domain] = (lat, lon, topo, border_mask)

print("Available domains:", list(domains_data.keys()))

alpx_key = None
for candidate in ['ALPX-3', 'ALPX3', 'ALPX_3', 'alpx-3']:
    if candidate in domains_data:
        alpx_key = candidate
        break

if alpx_key is None:
    raise KeyError(f"ALPX-3 not found. Available keys: {list(domains_data.keys())}")

lat_a, lon_a, topo_a, border_a = domains_data[alpx_key]

# Figure 
fig = plt.figure(figsize=(10, 8))
proj = ccrs.PlateCarree()

levels    = np.linspace(0, 4000, 21)
font_size = 11

ax   = fig.add_subplot(1, 1, 1, projection=proj)
ax.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.OCEAN, facecolor='#cce6f4', zorder=0)
ax.add_feature(cfeature.LAND,  facecolor='#f2efe9', zorder=1)
cs = ax.contourf(lon_a, lat_a, topo_a, levels=levels, cmap='terrain', extend='max', transform=ccrs.PlateCarree(), zorder=2)

# Gray border band
ax.contourf(lon_a, lat_a, border_a, colors='gray', levels=[0, 1], alpha=0.85, transform=ccrs.PlateCarree(), zorder=3)

# Domain outline
ax.contour(lon_a, lat_a, border_a, levels=[0.5], colors='dimgray', linewidths=1.0, transform=ccrs.PlateCarree(), zorder=4)

# Country borders and coastlines
ax.add_feature(cfeature.COASTLINE, linewidth=0.9, edgecolor='#222222', zorder=5)
ax.add_feature(cfeature.BORDERS,   linestyle='-', linewidth=0.7, edgecolor='#444444', zorder=5)

# Gridlines with labels
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, color='gray', linestyle='--',linewidth=0.4, alpha=0.6, zorder=6)

gl.top_labels   = False
gl.right_labels = False
gl.xlabel_style = {'size': font_size - 1}
gl.ylabel_style = {'size': font_size - 1}

# North arrow
ax.text(0.05, 0.15, u'\u25B2\nN', transform=ax.transAxes, ha='center', va='top', fontsize=font_size, fontweight='bold', color='black')

# Colorbar
cbar = fig.colorbar(cs, ax=ax, orientation='horizontal', pad=0.04, shrink=0.80, aspect=28)
cbar.set_label('Topography (m)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size - 1)

#  Save
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/figs/workshop'
os.makedirs(path_out, exist_ok=True)
name_out = 'pyplt_maps_topography_ALPX3_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
