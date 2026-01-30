# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot study area"

import os
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

domain = 'SAM-22'
dt = '1970-1971'
latlon = [-85, -30, -60, 15]

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/{0}'.format(domain)

# Subdomains (lat_min, lat_max, lon_min, lon_max)
subdomains = {
    "AMZ": (-15, -5,  -68, -48),
    "LPB": (-33, -20, -62, -49),
    "AND": (-45, -35, -73, -69),
    "NEB": (-15, -2,  -45, -35)
}

# Plot figure
fig = plt.figure(figsize=(10, 10))

ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent(latlon, crs=ccrs.PlateCarree())
ax.set_xticks(np.arange(latlon[0], latlon[1], 20), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(latlon[2], latlon[3], 20), crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
ax.grid(c='k', ls='--', alpha=0.4)

ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='lightblue')

# Draw red boxes for subdomains
for name, (lat_min, lat_max, lon_min, lon_max) in subdomains.items():
    lats = [lat_min, lat_max, lat_max, lat_min, lat_min]
    lons = [lon_min, lon_min, lon_max, lon_max, lon_min]
    ax.plot(lons, lats, transform=ccrs.PlateCarree(), color='red', linewidth=2)

    ax.text(lon_min + 0.5, lat_max - 0.5, name, color='red', fontsize=10,
            transform=ccrs.PlateCarree(), bbox=dict(facecolor='white', alpha=0.6, edgecolor='none'))

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_maps_study_area_{0}_RegCM5_{1}.png'.format(domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

