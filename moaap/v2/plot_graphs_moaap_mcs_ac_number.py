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


def comp_annual_cycle_number(mcs_charac):

    all_hours = []
    count_by_month = np.zeros(12)

    for obj in mcs_charac.keys():

        for m in mcs_charac[obj].keys():

            times = mcs_charac[obj][m]['times']
            hours = pd.DatetimeIndex(times).hour
            
            # Count number of MCSs by month
            for h in range(12):
                count_by_month[h] += np.sum(hours == h)

    # Calculate mean number per month (divide by number of years)
    n_years = 10  # 2000-2009
    mean_by_month = count_by_month / n_years

    return mean_by_month


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
print("Loading MCS data...")
mcs_car_era5 = open_mcs_era5('CAR-4')
mcs_csam_era5 = open_mcs_era5('CSAM-3')
mcs_eurr_era5 = open_mcs_era5('EURR-3')

mcs_car_cpm = open_mcs_cpm('CAR-4')
mcs_csam_cpm = open_mcs_cpm('CSAM-3')
mcs_eurr_cpm = open_mcs_cpm('EURR-3')

# Calculate mean number of MCSs by month
print("Calculating MCS counts...")
num_mcs_car_era5 = comp_annual_cycle_number(mcs_car_era5)
num_mcs_csam_era5 = comp_annual_cycle_number(mcs_csam_era5)
num_mcs_eurr_era5 = comp_annual_cycle_number(mcs_eurr_era5)

num_mcs_car_cpm = comp_annual_cycle_number(mcs_car_cpm)
num_mcs_csam_cpm = comp_annual_cycle_number(mcs_csam_cpm)
num_mcs_eurr_cpm = comp_annual_cycle_number(mcs_eurr_cpm)

# Print values
print("\nCAR-4 Domain - Mean number of MCSs per month:")
print("Month    ERA5    RegCM5")
for i in range(12):
    print(f"{i+1:2d}     {num_mcs_car_era5[i]:6.1f}  {num_mcs_car_cpm[i]:6.1f}")

print("\nCSAM-3 Domain - Mean number of MCSs per month:")
print("Month    ERA5    RegCM5")
for i in range(12):
    print(f"{i+1:2d}     {num_mcs_csam_era5[i]:6.1f}  {num_mcs_csam_cpm[i]:6.1f}")

print("\nEURR-3 Domain - Mean number of MCSs per month:")
print("Month    ERA5    RegCM5")
for i in range(12):
    print(f"{i+1:2d}     {num_mcs_eurr_era5[i]:6.1f}  {num_mcs_eurr_cpm[i]:6.1f}")

# plot
fig = plt.figure(figsize=(18, 6))
font_size = 10

width = 0.35
time = np.arange(0, 12)
xtick = ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D')

# CAR-4
ax = fig.add_subplot(1, 3, 1)
ax.bar(time - width/2, num_mcs_car_era5, width=width, color='red', alpha=0.75, edgecolor='red', label='ERA5')
ax.bar(time + width/2, num_mcs_car_cpm, width=width, color='blue', alpha=0.75, edgecolor='blue', label='RegCM5')
plot_domain_inset(ax, 'CAR-4')
plt.title('(a) CAR-4', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Number of MCSs', fontsize=font_size, fontweight='bold')
plt.xticks(time, xtick, fontsize=font_size)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.legend(loc=1, ncol=2, frameon=False, fontsize=font_size)

# CSAM-3
ax = fig.add_subplot(1, 3, 2)
ax.bar(time - width/2, num_mcs_csam_era5, width=width, color='red', alpha=0.75, edgecolor='red', label='ERA5')
ax.bar(time + width/2, num_mcs_csam_cpm, width=width, color='blue', alpha=0.75, edgecolor='blue', label='RegCM5')
plot_domain_inset(ax, 'CSAM-3')
plt.title('(b) CSAM-3', loc='left', fontsize=font_size, fontweight='bold')
plt.xticks(time, xtick, fontsize=font_size)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# EURR-3
ax = fig.add_subplot(1, 3, 3)
ax.bar(time - width/2, num_mcs_eurr_era5, width=width, color='red', alpha=0.75, edgecolor='red', label='ERA5')
ax.bar(time + width/2, num_mcs_eurr_cpm, width=width, color='blue', alpha=0.75, edgecolor='blue', label='RegCM5')
plot_domain_inset(ax, 'EURR-3')
plt.title('(c) EURR-3', loc='left', fontsize=font_size, fontweight='bold')
plt.xticks(time, xtick, fontsize=font_size)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Find peak months
car_peak_era5 = np.argmax(num_mcs_car_era5) + 1
car_peak_cpm = np.argmax(num_mcs_car_cpm) + 1
csam_peak_era5 = np.argmax(num_mcs_csam_era5) + 1
csam_peak_cpm = np.argmax(num_mcs_csam_cpm) + 1
eurr_peak_era5 = np.argmax(num_mcs_eurr_era5) + 1
eurr_peak_cpm = np.argmax(num_mcs_eurr_cpm) + 1

print(f"\nPeak months (mean number of MCSs):")
print(f"CAR-4:  ERA5={car_peak_era5}, RegCM5={car_peak_cpm}")
print(f"CSAM-3: ERA5={csam_peak_era5}, RegCM5={csam_peak_cpm}")
print(f"EURR-3: ERA5={eurr_peak_era5}, RegCM5={eurr_peak_cpm}")

# Path out to save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/figs'
name_out = f'pyplt_graphs_moaap_mcs_number_ac_domains_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
