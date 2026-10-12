# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "April 14, 2026"
__description__ = "This script plots domain topography with manual plot order"

import os
import sys
import glob
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

from netCDF4 import Dataset as nc

dirnc = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/dataset/domains'
number_i = 29


def import_rcm(dirnc, file, number):

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

    lat = np.squeeze(lat)
    lon = np.squeeze(lon)
    topo = np.squeeze(topo)

    if lat.ndim == 1 and lon.ndim == 1:
        lon, lat = np.meshgrid(lon, lat)

    topo = np.where(topo < 0, 0, topo)

    ny, nx = topo.shape
    border_mask = np.full((ny, nx), np.nan)
    border_mask[:number, :] = 1
    border_mask[-number:, :] = 1
    border_mask[:, :number] = 1
    border_mask[:, -number:] = 1

    return lat, lon, topo, border_mask


# Import dataset into a dictionary by domain name
files = sorted(glob.glob(os.path.join(dirnc, 'orog_*.nc')))
domains_data = {}

for file in files:
    domain = os.path.basename(file).split('_')[1]
    lat, lon, topo, border_mask = import_rcm(dirnc, file, number_i)
    domains_data[domain] = (lat, lon, topo, border_mask)


# Define your custom plotting order here:
domain_order = ['EURR-3', 'CAR-4', 'CSAM-3', 'NCV-3', 'PFR-3', 'ALPX-3', 'ANT-3']


# Figure Setup
fig, ax = plt.subplots(figsize=(14, 8), subplot_kw={'projection': ccrs.PlateCarree()})
font_size = 9

ax.stock_img()
ax.coastlines(linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)

levels = np.linspace(0, 4000, 21)

# Manual sequential plotting based on chosen order
for dom in domain_order:
    if dom in domains_data:
        lat, lon, topo, border_mask = domains_data[dom]

        # Plot topography
        cs = ax.contourf(
            lon, lat, topo, 
            levels=levels, 
            cmap='terrain', 
            extend='max', 
            transform=ccrs.PlateCarree()
        )

        # Plot white borders
        ax.contourf(
            lon, lat, border_mask, 
            colors='gray', 
            levels=[0, 1], 
            alpha=0.9, 
            transform=ccrs.PlateCarree()
        )

# North Arrow
plt.text(140, -80, u'\u25B2 \nN', color='black', fontsize=font_size + 2, fontweight='bold', transform=ccrs.PlateCarree())

# Colorbar for Topography
cbar = fig.colorbar(cs, ax=ax, orientation='horizontal', pad=0.06, shrink=0.7, aspect=30)
cbar.set_label('Topography (m)', fontsize=font_size + 1, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Gridlines
gridlines = ax.gridlines(draw_labels=True, color='gray', linestyle='--', alpha=0.6)
gridlines.right_labels = False
gridlines.bottom_labels = True
gridlines.left_labels = True
gridlines.top_labels = False
gridlines.xlabel_style = {'size': font_size}
gridlines.ylabel_style = {'size': font_size}

# Save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/figs/paper'
os.makedirs(path_out, exist_ok=True)
name_out = 'pyplt_maps_topography_CP_domains.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
