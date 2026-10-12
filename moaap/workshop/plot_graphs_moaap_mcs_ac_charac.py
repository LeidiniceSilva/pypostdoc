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
from matplotlib.patches import FancyBboxPatch

warnings.filterwarnings("ignore")


def open_mcs(path_, start="2000-01", end="2009-12"):

    path = f"/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/dataset/{path_}"
    dates = pd.date_range(start=start, end=end, freq="MS")
    mcs = {}

    for d in dates:
        f = "MCSs_" + d.strftime("%Y%m") + "__dt-1h_MOAAP-masks.pkl"
        f = os.path.join(path, f)

        if os.path.exists(f):
            with open(f, "rb") as file:
                mcs[d.strftime("%Y-%m")] = pickle.load(file)

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

    # Calculate mean number per month 
    n_years = 10  # 2000-2009
    mean_by_month = count_by_month / n_years

    return mean_by_month


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


# Domain
domain = "EUR"

# Load datasets
mcs_eur_gpm = open_mcs("/GPM/EURR-3/output", start="2000-01", end="2009-12")
mcs_eur_cpm_eval = open_mcs("/CPMs/ICTP/EURR-3/evaluation/ERA5/output", start="2000-01", end="2009-12")
mcs_eur_cpm_hist = open_mcs("/CPMs/ICTP/EURR-3/historical/ECEarth/output", start="2000-01", end="2009-12")
mcs_eur_rcm_eval = open_mcs("/RCMs/ICTP/EUR-12/evaluation/ERA5/output", start="2000-01", end="2009-12")
mcs_eur_rcm_hist = open_mcs("/RCMs/ICTP/EUR-12/historical/ECEarth/output", start="2000-01", end="2009-12")

# Calculate metrics for each dataset
num_gpm = comp_annual_cycle_number(mcs_eur_gpm)
num_cpm_eval = comp_annual_cycle_number(mcs_eur_cpm_eval)
num_cpm_hist = comp_annual_cycle_number(mcs_eur_cpm_hist)
num_rcm_eval = comp_annual_cycle_number(mcs_eur_rcm_eval)
num_rcm_hist = comp_annual_cycle_number(mcs_eur_rcm_hist)

size_gpm, tot_gpm, max_gpm = comp_annual_cycle(mcs_eur_gpm)
size_cpm_eval, tot_cpm_eval, max_cpm_eval = comp_annual_cycle(mcs_eur_cpm_eval)
size_cpm_hist, tot_cpm_hist, max_cpm_hist = comp_annual_cycle(mcs_eur_cpm_hist)
size_rcm_eval, tot_rcm_eval, max_rcm_eval = comp_annual_cycle(mcs_eur_rcm_eval)
size_rcm_hist, tot_rcm_hist, max_rcm_hist = comp_annual_cycle(mcs_eur_rcm_hist)

# Plot parameters
fig = plt.figure(figsize=(18, 14))
font_size = 10
gs = gridspec.GridSpec(4, 1, figure=fig, hspace=0.35)

time = np.arange(0, 12)
xtick = ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D')
width = 0.15  # width of each bar

labels = ['GPM', 'CPM-3 Eval', 'CPM-3 Hist', 'RCM-12 Eval', 'RCM-12 Hist']
colors = ['black', 'red', 'blue', 'green', 'orange']

# 1) Count plot
ax1 = fig.add_subplot(gs[0, 0])
ax1.bar(time - 2*width, num_gpm/6, width, label=labels[0], color=colors[0], alpha=0.75, edgecolor='white', linewidth=1)
ax1.bar(time - width, num_cpm_eval, width, label=labels[1], color=colors[1], alpha=0.75, edgecolor='white', linewidth=1)
ax1.bar(time, num_cpm_hist, width, label=labels[2], color=colors[2], alpha=0.75, edgecolor='white', linewidth=1)
ax1.bar(time + width, num_rcm_eval, width, label=labels[3], color=colors[3], alpha=0.75, edgecolor='white', linewidth=1)
ax1.bar(time + 2*width, num_rcm_hist, width, label=labels[4], color=colors[4], alpha=0.75, edgecolor='white', linewidth=1)
ax1.set_title('(a) MCS count', fontsize=font_size+2, loc='left', fontweight='bold')
ax1.set_ylabel('Mean MCS count', fontsize=font_size)
ax1.set_ylim(0, 50)
ax1.set_xticks(time)
ax1.set_xticklabels(xtick, fontsize=font_size)
ax1.grid(True, linestyle='--', alpha=0.5)

# 2) Size plot
ax2 = fig.add_subplot(gs[1, 0])
ax2.bar(time - 2*width, size_gpm, width, label=labels[0], color=colors[0], alpha=0.75, edgecolor='white', linewidth=1)
ax2.bar(time - width, size_cpm_eval, width, label=labels[1], color=colors[1], alpha=0.75, edgecolor='white', linewidth=1)
ax2.bar(time, size_cpm_hist, width, label=labels[2], color=colors[2], alpha=0.75, edgecolor='white', linewidth=1)
ax2.bar(time + width, size_rcm_eval, width, label=labels[3], color=colors[3], alpha=0.75, edgecolor='white', linewidth=1)
ax2.bar(time + 2*width, size_rcm_hist, width, label=labels[4], color=colors[4], alpha=0.75, edgecolor='white', linewidth=1)
ax2.set_title('(b) MCS Size', fontsize=font_size+2, loc='left', fontweight='bold')
ax2.set_ylabel('Size (10³ km²)', fontsize=font_size)
ax2.set_ylim(0, 400000)
ax2.set_xticks(time)
ax2.set_xticklabels(xtick, fontsize=font_size)
ax2.grid(True, linestyle='--', alpha=0.5)

# 3) Total Precipitation plot
ax3 = fig.add_subplot(gs[2, 0])
ax3.bar(time - 2*width, tot_gpm, width, label=labels[0], color=colors[0], alpha=0.75, edgecolor='white', linewidth=1)
ax3.bar(time - width, tot_cpm_eval, width, label=labels[1], color=colors[1], alpha=0.75, edgecolor='white', linewidth=1)
ax3.bar(time, tot_cpm_hist, width, label=labels[2], color=colors[2], alpha=0.75, edgecolor='white', linewidth=1)
ax3.bar(time + width, tot_rcm_eval, width, label=labels[3], color=colors[3], alpha=0.75, edgecolor='white', linewidth=1)
ax3.bar(time + 2*width, tot_rcm_hist, width, label=labels[4], color=colors[4], alpha=0.75, edgecolor='white', linewidth=1)
ax3.set_title('(c) MCS total precipitation', fontsize=font_size+2, loc='left', fontweight='bold')
ax3.set_ylabel('Total precipitation (mm)', fontsize=font_size)
ax3.set_ylim(0, 1000)
ax3.set_xticks(time)
ax3.set_xticklabels(xtick, fontsize=font_size)
ax3.grid(True, linestyle='--', alpha=0.5)

# 4) Maximum Intensity plot
ax4 = fig.add_subplot(gs[3, 0])
ax4.bar(time - 2*width, max_gpm, width, label=labels[0], color=colors[0], alpha=0.75, edgecolor='white', linewidth=1)
ax4.bar(time - width, max_cpm_eval, width, label=labels[1], color=colors[1], alpha=0.75, edgecolor='white', linewidth=1)
ax4.bar(time, max_cpm_hist, width, label=labels[2], color=colors[2], alpha=0.75, edgecolor='white', linewidth=1)
ax4.bar(time + width, max_rcm_eval, width, label=labels[3], color=colors[3], alpha=0.75, edgecolor='white', linewidth=1)
ax4.bar(time + 2*width, max_rcm_hist, width, label=labels[4], color=colors[4], alpha=0.75, edgecolor='white', linewidth=1)
ax4.set_title('(d) MCS max precipitation intensity', fontsize=font_size+2, loc='left', fontweight='bold')
ax4.set_xlabel('Month', fontsize=font_size)
ax4.set_ylabel('Max Precip (mm/h)', fontsize=font_size)
ax4.set_ylim(0, 30)
ax4.set_xticks(time)
ax4.set_xticklabels(xtick, fontsize=font_size)
ax4.grid(True, linestyle='--', alpha=0.5)
ax4.legend(loc='upper right', fontsize=font_size, ncol=5)

# Path out to save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/figs/v2'
name_out = f'pyplt_graphs_moaap_mcs_ac_charac_{domain}_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
