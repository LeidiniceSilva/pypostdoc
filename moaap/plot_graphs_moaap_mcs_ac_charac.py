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
import cartopy.crs as ccrs
import cartopy.feature as cfeature
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


def open_mcs_cpm(domain, start='2000-01', end='2009-12'):

    path = f'/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/CPMs/{domain}/output/'

    dates = pd.date_range(start=start, end=end, freq='MS')
    mcs = {}

    for d in dates:
        f = 'MCSs_' + d.strftime('%Y%m') + '__dt-1h_MOAAP-masks.pkl'
        f = os.path.join(path, f)

        if os.path.exists(f):
            with open(f, 'rb') as file:
                mcs[d.strftime('%Y-%m')] = pickle.load(file)

    return mcs


def comp_annual_cycle(mcs_charac):

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

    size_clim = np.full(12, np.nan)
    tot_clim  = np.full(12, np.nan)
    max_clim  = np.full(12, np.nan)

    for h in range(12):

        sel = all_hours == h

        if np.any(sel):
            size_clim[h] = np.nanmean(all_size[sel])
            tot_clim[h]  = np.nanmean(all_tot[sel])
            max_clim[h]  = np.nanmean(all_max[sel])

    return size_clim, tot_clim, max_clim



def plot_domain_inset(ax, domain):

    domains = {
        'CAR-4':  [-119.0, -58.25, 9.25, 35.75],
        'CSAM-3': [-78.75, -35.5, -36.5, -12.25],
        'EURR-3': [-25.25, 38.25, 33.5, 64.75]
    }

    lon_min, lon_max, lat_min, lat_max = domains[domain]

    if domain == 'CAR-4':
        inset_ax = ax.inset_axes([0.05, 0.7, 0.4, 0.5], projection=ccrs.PlateCarree())
    elif domain == 'CSAM-3':
        inset_ax = ax.inset_axes([0.05, 0.7, 0.4, 0.5], projection=ccrs.PlateCarree())
    else:
        inset_ax = ax.inset_axes([0.05, 0.7, 0.4, 0.5], projection=ccrs.PlateCarree())

    inset_ax.set_extent([lon_min, lon_max,
                         lat_min, lat_max],
                         crs=ccrs.PlateCarree())

    inset_ax.add_feature(cfeature.OCEAN, facecolor='#a6cbe3')
    inset_ax.add_feature(cfeature.LAND, facecolor='#e6d2b5')

    inset_ax.coastlines(linewidth=0.5)
    inset_ax.add_feature(cfeature.BORDERS, linewidth=0.25)

    inset_ax.plot(
        [lon_min, lon_max, lon_max, lon_min, lon_min],
        [lat_min, lat_min, lat_max, lat_max, lat_min],
        color='black', linewidth=0.5,
        transform=ccrs.PlateCarree()
    )

    inset_ax.set_xticks([])
    inset_ax.set_yticks([])


# load data
mcs_car_era5 = open_mcs_era5('CAR-4')
mcs_csam_era5 = open_mcs_era5('CSAM-3')
mcs_eurr_era5 = open_mcs_era5('EURR-3')

mcs_car_cpm = open_mcs_cpm('CAR-4')
mcs_csam_cpm = open_mcs_cpm('CSAM-3')
mcs_eurr_cpm = open_mcs_cpm('EURR-3')

size_mcs_car_era5, tot_mcs_car_era5, max_mcs_car_era5 = comp_annual_cycle(mcs_car_era5)
size_mcs_csam_era5, tot_mcs_csam_era5, max_mcs_csam_era5 = comp_annual_cycle(mcs_csam_era5)
size_mcs_eurr_era5, tot_mcs_eurr_era5, max_mcs_eurr_era5 = comp_annual_cycle(mcs_eurr_era5)

size_mcs_car_cpm, tot_mcs_car_cpm, max_mcs_car_cpm = comp_annual_cycle(mcs_car_cpm)
size_mcs_csam_cpm, tot_mcs_csam_cpm, max_mcs_csam_cpm = comp_annual_cycle(mcs_csam_cpm)
size_mcs_eurr_cpm, tot_mcs_eurr_cpm, max_mcs_eurr_cpm = comp_annual_cycle(mcs_eurr_cpm)

# plot
fig = plt.figure(figsize=(18, 12))
font_size = 10

width = 0.35
time = np.arange(0, 12)
xtick = ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D')

ax = fig.add_subplot(3, 3, 1)
ax.bar(time - width/2, size_mcs_car_era5, width=width, color='red', alpha=0.75, edgecolor='red', label='ERA5')
ax.bar(time + width/2, size_mcs_car_cpm, width=width, color='blue', alpha=0.75, edgecolor='blue', label='RegCM5')
plot_domain_inset(ax, 'CAR-4')
plt.title('(a)', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('MCS anvil size (km$^2$)', fontsize=font_size, fontweight='bold')
plt.yticks(np.arange(0, 550000, 50000), fontsize=font_size)
plt.xticks(time, xtick, fontsize=font_size)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.legend(loc=1, ncol=2, frameon=False, fontsize=font_size)

ax = fig.add_subplot(3, 3, 2)
ax.bar(time - width/2, size_mcs_csam_era5, width=width, color='red', alpha=0.75, edgecolor='red', label='ERA5')
ax.bar(time + width/2, size_mcs_csam_cpm, width=width, color='blue', alpha=0.75, edgecolor='blue', label='RegCM5')
plot_domain_inset(ax, 'CSAM-3')
plt.title('(b)', loc='left', fontsize=font_size, fontweight='bold')
plt.yticks(np.arange(0, 550000, 50000), fontsize=font_size)
plt.xticks(time, xtick, fontsize=font_size)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax = fig.add_subplot(3, 3, 3)
ax.bar(time - width/2, size_mcs_eurr_era5, width=width, color='red', alpha=0.75, edgecolor='red', label='ERA5')
ax.bar(time + width/2, size_mcs_eurr_cpm, width=width, color='blue', alpha=0.75, edgecolor='blue', label='RegCM5')
plot_domain_inset(ax, 'EURR-3')
plt.title('(c)', loc='left', fontsize=font_size, fontweight='bold')
plt.yticks(np.arange(0, 550000, 50000), fontsize=font_size)
plt.xticks(time, xtick, fontsize=font_size)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax = fig.add_subplot(3, 3, 4)
ax.bar(time - width/2, tot_mcs_car_era5, width=width, color='red', alpha=0.75, edgecolor='red', label='ERA5')
ax.bar(time + width/2, tot_mcs_car_cpm, width=width, color='blue', alpha=0.75, edgecolor='blue', label='RegCM5')
plt.title('(d)', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation volume (km$^{3}$ h$^{-1}$)', fontsize=font_size, fontweight='bold')
plt.yticks(np.arange(0, 800, 100), fontsize=font_size)
plt.xticks(time, xtick, fontsize=font_size)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax = fig.add_subplot(3, 3, 5)
ax.bar(time - width/2, tot_mcs_csam_era5, width=width, color='red', alpha=0.75, edgecolor='red', label='ERA5')
ax.bar(time + width/2, tot_mcs_csam_cpm, width=width, color='blue', alpha=0.75, edgecolor='blue', label='RegCM5')
plt.title('(e)', loc='left', fontsize=font_size, fontweight='bold')
plt.yticks(np.arange(0, 800, 100), fontsize=font_size)
plt.xticks(time, xtick, fontsize=font_size)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax = fig.add_subplot(3, 3, 6)
ax.bar(time - width/2, tot_mcs_eurr_era5, width=width, color='red', alpha=0.75, edgecolor='red', label='ERA5')
ax.bar(time + width/2, tot_mcs_eurr_cpm, width=width, color='blue', alpha=0.75, edgecolor='blue', label='RegCM5')
plt.title('(f)', loc='left', fontsize=font_size, fontweight='bold')
plt.yticks(np.arange(0, 800, 100), fontsize=font_size)
plt.xticks(time, xtick, fontsize=font_size)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax = fig.add_subplot(3, 3, 7)
ax.bar(time - width/2, max_mcs_car_era5, width=width, color='red', alpha=0.75, edgecolor='red', label='ERA5')
ax.bar(time + width/2, max_mcs_car_cpm, width=width, color='blue', alpha=0.75, edgecolor='blue', label='RegCM5')
plt.title('(g)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('time', fontsize=font_size, fontweight='bold')
plt.ylabel('Max. precipitation (mm h$^{-1}$)', fontsize=font_size, fontweight='bold')
plt.yticks(np.arange(0, 30, 5), fontsize=font_size)
plt.xticks(time, xtick, fontsize=font_size)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax = fig.add_subplot(3, 3, 8)
ax.bar(time - width/2, max_mcs_csam_era5, width=width, color='red', alpha=0.75, edgecolor='red', label='ERA5')
ax.bar(time + width/2, max_mcs_csam_cpm, width=width, color='blue', alpha=0.75, edgecolor='blue', label='RegCM5')
plt.title('(h)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('time', fontsize=font_size, fontweight='bold')
plt.ylim(0, 25)
plt.yticks(np.arange(0, 30, 5), fontsize=font_size)
plt.xticks(time, xtick, fontsize=font_size)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax = fig.add_subplot(3, 3, 9)
ax.bar(time - width/2, max_mcs_eurr_era5, width=width, color='red', alpha=0.75, edgecolor='red', label='ERA5')
ax.bar(time + width/2, max_mcs_eurr_cpm, width=width, color='blue', alpha=0.75, edgecolor='blue', label='RegCM5')
plt.title('(i)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('time', fontsize=font_size, fontweight='bold')
plt.yticks(np.arange(0, 30, 5), fontsize=font_size)
plt.xticks(time, xtick, fontsize=font_size)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Path out to save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/figs'
name_out = f'pyplt_graphs_moaap_mcs_charac_ac_domains_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()



