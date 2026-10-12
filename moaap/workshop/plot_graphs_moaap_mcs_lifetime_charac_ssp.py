# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 03, 2026"
__description__ = "This script plot MCSs"

import os
import glob
import pickle
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch

warnings.filterwarnings("ignore")


def open_mcs(path_, start="2000-01", end="2009-12"):
    """Reads MCS pickle files for a given date range and path."""
    base_dir = "/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/dataset"
    path = os.path.join(base_dir, path_.lstrip("/"))
    dates = pd.date_range(start=start, end=end, freq="MS")
    mcs = {}

    for d in dates:
        pattern = os.path.join(path, f"MCSs_{d.strftime('%Y%m')}*__dt-1h_MOAAP-masks.pkl")
        matches = glob.glob(pattern)

        if matches:
            with open(matches[0], "rb") as file:
                mcs[d.strftime("%Y-%m")] = pickle.load(file)

    return mcs


def comp_lifetime_number(mcs_dict):
    """Calculates mean number of MCSs per year by lifetime hour (5-20h)."""
    count_by_hour = {h: 0 for h in range(5, 21)}

    for month in mcs_dict.keys():
        systems = mcs_dict[month]
        # Handle dict or list structure of systems
        if isinstance(systems, dict):
            system_list = systems.values()
        elif isinstance(systems, list):
            system_list = systems
        else:
            continue

        for mcs in system_list:
            times = mcs.get('times', [])
            lifetime = len(times)

            if 5 <= lifetime <= 20:
                for hour_idx in range(lifetime):
                    lifetime_hour = hour_idx + 1
                    if 5 <= lifetime_hour <= 20:
                        count_by_hour[lifetime_hour] += 1

    # Number of years in dataset (based on unique years loaded)
    years = set(m.split("-")[0] for m in mcs_dict.keys())
    n_years = len(years) if len(years) > 0 else 10

    hours = np.array(list(range(5, 21)))
    counts = np.array([count_by_hour[h] for h in hours]) / n_years

    return counts, hours


def comp_lifetime_cycle(mcs_dict):
    """Calculates average Size, Total Precip, and Max Intensity by lifetime hour (5-20h)."""
    size_by_hour = {h: [] for h in range(5, 21)}
    tot_by_hour = {h: [] for h in range(5, 21)}
    max_by_hour = {h: [] for h in range(5, 21)}

    for month in mcs_dict.keys():
        systems = mcs_dict[month]
        if isinstance(systems, dict):
            system_list = systems.values()
        elif isinstance(systems, list):
            system_list = systems
        else:
            continue

        for mcs in system_list:
            times = mcs.get('times', [])
            lifetime = len(times)

            if 5 <= lifetime <= 20:
                sizes = mcs.get('size', [])
                tots = mcs.get('tot', [])
                maxs = mcs.get('max', [])

                for hour_idx in range(lifetime):
                    lifetime_hour = hour_idx + 1

                    if 5 <= lifetime_hour <= 20 and hour_idx < len(sizes):
                        size_by_hour[lifetime_hour].append(sizes[hour_idx] / 1000**2)
                        tot_by_hour[lifetime_hour].append(tots[hour_idx])
                        max_by_hour[lifetime_hour].append(maxs[hour_idx])

    size_clim, tot_clim, max_clim, hours = [], [], [], []

    for h in range(5, 21):
        if len(size_by_hour[h]) > 0:
            size_clim.append(np.nanmean(size_by_hour[h]))
            tot_clim.append(np.nanmean(tot_by_hour[h]))
            max_clim.append(np.nanmean(max_by_hour[h]))
        else:
            size_clim.append(np.nan)
            tot_clim.append(np.nan)
            max_clim.append(np.nan)
        hours.append(h)

    return np.array(size_clim), np.array(tot_clim), np.array(max_clim), np.array(hours)


# Domain
domain = "EUR"

# Path definitions
rcm_hist_path = "RCMs/ICTP/EUR-12/historical/ECEarth/output"
rcm_ssp_path  = "RCMs/ICTP/EUR-12/ssp370/ECEarth/output"

# 1. Historical (2000-01 to 2009-12)
mcs_rcm_hist = open_mcs(rcm_hist_path, start="2000-01", end="2009-12")

# 2. GWL 1.5 (Historical 2000-01 to 2009-12 + SSP370 2015-01 to 2021-12)
mcs_rcm_gwl15_part1 = open_mcs(rcm_hist_path, start="2002-01", end="2009-12")
mcs_rcm_gwl15_part2 = open_mcs(rcm_ssp_path, start="2015-01", end="2021-12")
mcs_rcm_gwl15 = {**mcs_rcm_gwl15_part1, **mcs_rcm_gwl15_part2}

# 3. GWL 2.0 (SSP370 2023-01 to 2042-12)
mcs_rcm_gwl20 = open_mcs(rcm_ssp_path, start="2023-01", end="2042-12")

# Calculate Count metrics
num_hist, hours  = comp_lifetime_number(mcs_rcm_hist)
num_gwl15, _     = comp_lifetime_number(mcs_rcm_gwl15)
num_gwl20, _     = comp_lifetime_number(mcs_rcm_gwl20)

# Calculate Lifetime Cycle metrics
size_hist, tot_hist, max_hist, _    = comp_lifetime_cycle(mcs_rcm_hist)
size_gwl15, tot_gwl15, max_gwl15, _ = comp_lifetime_cycle(mcs_rcm_gwl15)
size_gwl20, tot_gwl20, max_gwl20, _ = comp_lifetime_cycle(mcs_rcm_gwl20)

# Plot parameters
fig = plt.figure(figsize=(18, 10))
font_size = 10
gs = gridspec.GridSpec(3, 1, figure=fig, hspace=0.35)

x = np.arange(len(hours))
width = 0.25

labels = ['Historical', 'GWL 1.5°C', 'GWL 2.0°C']
colors = ['#1f77b4', '#ff7f0e', '#d62728']  # Blue, Orange, Red

# 1) Size plot
ax1 = fig.add_subplot(gs[0, 0])
ax1.bar(x - width, size_hist, width, label=labels[0], color=colors[0], edgecolor='white', linewidth=1)
ax1.bar(x, size_gwl15, width, label=labels[1], color=colors[1], edgecolor='white', linewidth=1)
ax1.bar(x + width, size_gwl20, width, label=labels[2], color=colors[2], edgecolor='white', linewidth=1)
ax1.set_title('(a) MCS size', fontsize=font_size+2, loc='left', fontweight='bold')
ax1.set_ylabel('Size (10³ km²)', fontsize=font_size)
ax1.set_xticks(x)
ax1.set_xticklabels(hours, fontsize=font_size)
ax1.grid(True, linestyle='--', alpha=0.5)
ax1.legend(loc='upper right', fontsize=font_size, ncol=3)

# 3) Total Precipitation plot
ax2 = fig.add_subplot(gs[1, 0])
ax2.bar(x - width, tot_hist, width, label=labels[0], color=colors[0], edgecolor='white', linewidth=1)
ax2.bar(x, tot_gwl15, width, label=labels[1], color=colors[1], edgecolor='white', linewidth=1)
ax2.bar(x + width, tot_gwl20, width, label=labels[2], color=colors[2], edgecolor='white', linewidth=1)
ax2.set_title('(b) MCS total precipitation', fontsize=font_size+2, loc='left', fontweight='bold')
ax2.set_ylabel('Total precipitation (mm)', fontsize=font_size)
ax2.set_xticks(x)
ax2.set_xticklabels(hours, fontsize=font_size)
ax2.grid(True, linestyle='--', alpha=0.5)

# 4) Maximum Intensity plot
ax3 = fig.add_subplot(gs[2, 0])
ax3.bar(x - width, max_hist, width, label=labels[0], color=colors[0], edgecolor='white', linewidth=1)
ax3.bar(x, max_gwl15, width, label=labels[1], color=colors[1], edgecolor='white', linewidth=1)
ax3.bar(x + width, max_gwl20, width, label=labels[2], color=colors[2], edgecolor='white', linewidth=1)
ax3.set_title('(c) MCS max precipitation intensity', fontsize=font_size+2, loc='left', fontweight='bold')
ax3.set_xlabel('Lifetime hour (h)', fontsize=font_size)
ax3.set_ylabel('Max precipitation (mm/h)', fontsize=font_size)
ax3.set_xticks(x)
ax3.set_xticklabels(hours, fontsize=font_size)
ax3.grid(True, linestyle='--', alpha=0.5)

# Save output figure
path_out = "/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/figs/v2"
name_out = f"pyplt_graphs_moaap_mcs_lifetime_charac_{domain}_GWL_2000-2009.png"
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.show()
exit()
