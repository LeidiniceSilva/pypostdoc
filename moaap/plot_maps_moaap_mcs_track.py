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


def open_moaap_era5(domain, start='2000-01', end='2009-12'):

    path = f'/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/ERA5/{domain}/output/'

    dates = pd.date_range(start=start, end=end, freq='MS')
    moaap = {}

    for d in dates:

        f = d.strftime('%Y%m%d') + '_ERA5_evaluation_ObjectMasks__dt-1h_MOAAP-masks.nc'
        f = os.path.join(path, f)

        if os.path.exists(f):
            moaap[d.strftime('%Y-%m')] = xr.open_dataset(f)

    return moaap


def open_mcs_era5(domain, start='2000-01', end='2009-12'):

    path = f'/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/ERA5/{domain}/output/'

    dates = pd.date_range(start=start, end=end, freq='MS')
    mcs = {}

    for d in dates:
        f = 'MCSs_' + d.strftime('%Y%m') + '__dt-1h_MOAAP-masks.pkl'
        f = os.path.join(path, f)

        if os.path.exists(f):
            with open(f, 'rb') as file:
                mcs[d.strftime('%Y-%m')] = pickle.load(file)

    return mcs


# load data
data_moaap_car = open_moaap_era5('CAR-4')

# load MCS characteristics
mcs_moaap_car = open_mcs_era5('CAR-4')

# plot
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(14,6))

# =========================
# PLOT ALL MONTHS TRACKS (NEW)
# =========================
for month in mcs_moaap_car.keys():

    mcs_charac = mcs_moaap_car[month]

    for ii in range(len(mcs_charac.keys())):
        LatLonTrack = mcs_charac[list(mcs_charac.keys())[ii]]['track']
        plt.plot(
            LatLonTrack[:,1],
            LatLonTrack[:,0],
            transform=ccrs.PlateCarree(),
            lw=1,
            color='k'
        )

# Path out to save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/figs'
name_out = 'pyplt_maps_moaap_mcs_track_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
