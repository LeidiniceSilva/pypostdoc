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
import matplotlib.gridspec as gridspec

warnings.filterwarnings("ignore")

# load MCS characteristics
data_moaap = xr.open_dataset('/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/ERA5/CAR-4/output/20000101_ERA5_evaluation_ObjectMasks__dt-1h_MOAAP-masks.nc')

with open('/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/ERA5/CAR-4/output/MCSs_200001__dt-1h_MOAAP-masks.pkl', 'rb') as f:
    mcs_charac = pickle.load(f)

# Plot figure
fig = plt.figure(figsize=(15,4))
gs1 = gridspec.GridSpec(1,3)
gs1.update(left=0.01, right=0.98,
           bottom=0.1, top=0.95,
           wspace=0.2, hspace=0.2)

obj = '1'

# MCS anvile size
ax = plt.subplot(gs1[0,0])
plt.title('Development of anvil size of MCS #' + obj)
plt.plot(mcs_charac[obj]['times'], mcs_charac[obj]['size']/1000**2)

ax.set_xlabel('time [hours]')
plt.xticks(rotation=45)
ax.set_ylabel('MCS anvil size [km$^2$]')

# MCS volume precipitation
ax = plt.subplot(gs1[0,1])
plt.title('Precipitation volume of MCS #' + obj)
plt.plot(mcs_charac[obj]['times'], mcs_charac[obj]['tot'])
ax.set_xlabel('time [hours]')
plt.xticks(rotation=45)
ax.set_ylabel('precipitation volume [km$^{3}$ h$^{-1}$]')

# MCS maximum precipitation
ax = plt.subplot(gs1[0,2])
plt.title('Max. precipitation of MCS #' + obj)
plt.plot(mcs_charac[obj]['times'], mcs_charac[obj]['max'])
ax.set_xlabel('time [hours]')
plt.xticks(rotation=45)
ax.set_ylabel('precipitation [mm h$^{-1}$]')

# Path out to save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/figs'
name_out = 'pyplt_maps_moaap_mcs_charac_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


