# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 03, 2026"
__description__ = "This script plot MCSs"

import os
import pickle
import netCDF4
import cartopy
import warnings
import numpy as np
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# Import data
data_vars = xr.open_dataset('/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/ERA5/CAR-4/input/CAR-4_ERA5_reanalysis_1hr_2000010100.nc')
time_datetime = pd.to_datetime(np.array(data_vars['time'].values, dtype='datetime64'))

# load MCS characteristics
data_moaap = xr.open_dataset('/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/ERA5/CAR-4/output/20000101_ERA5_evaluation_ObjectMasks__dt-1h_MOAAP-masks.nc')

with open('/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/ERA5/CAR-4/output/MCSs_200001__dt-1h_MOAAP-masks.pkl', 'rb') as f:
    mcs_charac = pickle.load(f)

# Plot figure
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(14,6))

# Set the extent of the map 
ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
ax.coastlines(color='#969696')
ax.gridlines()

# Generate some random data (latitude, longitude, and values)
mcs_mask = np.array(data_moaap['MCS_Tb_Objects'][12,:,:])
mcs_mask[mcs_mask == 0] = np.nan
sc = plt.pcolormesh(data_moaap['lon'],
                    data_moaap['lat'],
                    mcs_mask,
                    cmap = 'nipy_spectral')

# plot MCS tracks
for ii in range(len(mcs_charac.keys())):
  LatLonTrack = mcs_charac[list(mcs_charac.keys())[ii]]['track']
  plt.plot(LatLonTrack[:,1],LatLonTrack[:,0], transform=ccrs.PlateCarree(), lw=1, color='k')

ax.set_extent([-180, 180, -70, 70], crs=ccrs.PlateCarree())
cbar = plt.colorbar(sc, ax=ax, orientation='vertical', shrink=0.7, label='MCS mask')
plt.title('MCS tracks (black lines) and masks at '+str(time_datetime[12])[:16])

# Path out to save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/figs'
name_out = 'pyplt_maps_moaap_mcs_track_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


