# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 03, 2026"
__description__ = "This script plot MCSs"

import os
import pickle
import warnings
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

warnings.filterwarnings("ignore")


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


def comp_diurnal_cycle(mcs_charac):

    all_size, all_tot, all_max, all_hours = [], [], [], []

    for obj in mcs_charac.keys():

        for m in mcs_charac[obj].keys():

            times = mcs_charac[obj][m]['times']

            size = mcs_charac[obj][m]['size'] / 1000**2
            tot  = mcs_charac[obj][m]['tot']
            maxv = mcs_charac[obj][m]['max']

            hours = pd.DatetimeIndex(times).hour

            all_size.append(size)
            all_tot.append(tot)
            all_max.append(maxv)
            all_hours.append(hours)

    all_size = np.concatenate(all_size)
    all_tot  = np.concatenate(all_tot)
    all_max  = np.concatenate(all_max)
    all_hours = np.concatenate(all_hours)

    size_clim = np.full(24, np.nan)
    tot_clim  = np.full(24, np.nan)
    max_clim  = np.full(24, np.nan)

    for h in range(24):

        sel = all_hours == h

        if np.any(sel):
            size_clim[h] = np.nanmean(all_size[sel])
            tot_clim[h]  = np.nanmean(all_tot[sel])
            max_clim[h]  = np.nanmean(all_max[sel])

    return size_clim, tot_clim, max_clim


# load data
mcs_charac_car = open_mcs_era5('CAR-4')
mcs_charac_csam = open_mcs_era5('CSAM-3')
mcs_charac_eurr = open_mcs_era5('EURR-3')

size_clim_car, tot_clim_car, max_clim_car = comp_diurnal_cycle(mcs_charac_car)
size_clim_csam, tot_clim_csam, max_clim_csam = comp_diurnal_cycle(mcs_charac_csam)
size_clim_eurr, tot_clim_eurr, max_clim_eurr = comp_diurnal_cycle(mcs_charac_eurr)

# plot
fig = plt.figure(figsize=(17, 14))
time = np.arange(0, 24)
font_size = 10

ax = fig.add_subplot(3, 3, 1)
plt.plot(time, size_clim_car, marker='.', linewidth=1, color='red', label='ERA5')
plt.title('(a)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('CAR-4', fontsize=font_size, fontweight='bold')
plt.ylabel('MCS anvil size (km$^2$)', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.grid(True, linestyle='--', linewidth=1)
plt.legend(fontsize=font_size)

ax = fig.add_subplot(3, 3, 2)
plt.plot(time, size_clim_csam, marker='.', linewidth=1, color='red', label='ERA5')
plt.title('(b)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('CSAM-3', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.grid(True, linestyle='--', linewidth=1)

ax = fig.add_subplot(3, 3, 3)
plt.plot(time, size_clim_eurr, marker='.', linewidth=1, color='red', label='ERA5')
plt.title('(c)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('EURR-3', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.grid(True, linestyle='--', linewidth=1)

ax = fig.add_subplot(3, 3, 4)
plt.plot(time, tot_clim_car, marker='.', linewidth=1, color='red', label='ERA5')
plt.title('(d)', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation volume (km$^{3}$ h$^{-1}$)', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.grid(True, linestyle='--', linewidth=1)

ax = fig.add_subplot(3, 3, 5)
plt.plot(time, tot_clim_csam, marker='.', linewidth=1, color='red', label='ERA5')
plt.title('(e)', loc='left', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.grid(True, linestyle='--', linewidth=1)

ax = fig.add_subplot(3, 3, 6)
plt.plot(time, tot_clim_eurr, marker='.', linewidth=1, color='red', label='ERA5')
plt.title('(f)', loc='left', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.grid(True, linestyle='--', linewidth=1)

ax = fig.add_subplot(3, 3, 7)
plt.plot(time, max_clim_car, marker='.', linewidth=1, color='red', label='ERA5')
plt.title('(g)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('time', fontsize=font_size, fontweight='bold')
plt.ylabel('Max. precipitation (mm h$^{-1}$)', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.grid(True, linestyle='--', linewidth=1)

ax = fig.add_subplot(3, 3, 8)
plt.plot(time, max_clim_csam, marker='.', linewidth=1, color='red', label='ERA5')
plt.title('(h)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('time', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.grid(True, linestyle='--', linewidth=1)

ax = fig.add_subplot(3, 3, 9)
plt.plot(time, max_clim_eurr, marker='.', linewidth=1, color='red', label='ERA5')
plt.title('(i)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('time', fontsize=font_size, fontweight='bold')
plt.xticks(time, ('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
plt.grid(True, linestyle='--', linewidth=1)

# Path out to save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/figs'
name_out = f'pyplt_maps_moaap_mcs_charac_domains_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()



