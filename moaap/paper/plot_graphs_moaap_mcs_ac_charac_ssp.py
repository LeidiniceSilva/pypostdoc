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


def comp_annual_cycle_number(mcs_dict):
    """Calculates mean number of MCSs per year for each month (Jan-Dec)."""
    count_by_month = np.zeros(12)

    for month_key, systems in mcs_dict.items():
        if isinstance(systems, dict):
            system_list = systems.values()
        elif isinstance(systems, list):
            system_list = systems
        else:
            continue

        for mcs in system_list:
            times = mcs.get('times', [])
            if len(times) > 0:
                # Month index (0 to 11) derived from initiation time
                start_month = pd.DatetimeIndex(times)[0].month - 1
                count_by_month[start_month] += 1

    # Number of unique years loaded
    years = set(m.split("-")[0] for m in mcs_dict.keys())
    n_years = len(years) if len(years) > 0 else 10

    mean_by_month = count_by_month / n_years
    return mean_by_month


def comp_annual_cycle(mcs_dict):
    """Calculates average Size, Total Precip, and Max Intensity for each month (Jan-Dec)."""
    size_by_month = {m: [] for m in range(12)}
    tot_by_month = {m: [] for m in range(12)}
    max_by_month = {m: [] for m in range(12)}

    for month_key, systems in mcs_dict.items():
        if isinstance(systems, dict):
            system_list = systems.values()
        elif isinstance(systems, list):
            system_list = systems
        else:
            continue

        for mcs in system_list:
            times = mcs.get('times', [])
            if len(times) > 0:
                months = pd.DatetimeIndex(times).month - 1
                sizes = np.array(mcs.get('size', [])) / 1000**2
                tots = np.array(mcs.get('tot', []))
                maxs = np.array(mcs.get('max', []))

                for i, m_idx in enumerate(months):
                    if i < len(sizes):
                        size_by_month[m_idx].append(sizes[i])
                        tot_by_month[m_idx].append(tots[i])
                        max_by_month[m_idx].append(maxs[i])

    size_clim = np.full(12, np.nan)
    tot_clim  = np.full(12, np.nan)
    max_clim  = np.full(12, np.nan)

    for m_idx in range(12):
        if len(size_by_month[m_idx]) > 0:
            size_clim[m_idx] = np.nanmean(size_by_month[m_idx])
            tot_clim[m_idx]  = np.nanmean(tot_by_month[m_idx])
            max_clim[m_idx]  = np.nanmean(max_by_month[m_idx])

    return size_clim, tot_clim, max_clim


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
num_hist  = comp_annual_cycle_number(mcs_rcm_hist)
num_gwl15 = comp_annual_cycle_number(mcs_rcm_gwl15)
num_gwl20 = comp_annual_cycle_number(mcs_rcm_gwl20)

# Calculate Annual Cycle metrics
size_hist, tot_hist, max_hist    = comp_annual_cycle(mcs_rcm_hist)
size_gwl15, tot_gwl15, max_gwl15 = comp_annual_cycle(mcs_rcm_gwl15)
size_gwl20, tot_gwl20, max_gwl20 = comp_annual_cycle(mcs_rcm_gwl20)

# Plot parameters
fig = plt.figure(figsize=(18, 14))
font_size = 10
gs = gridspec.GridSpec(4, 1, figure=fig, hspace=0.35)

time = np.arange(12)
xtick = ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D')
width = 0.25

labels = ['Historical', 'GWL 1.5°C', 'GWL 2.0°C']
colors = ['#1f77b4', '#ff7f0e', '#d62728']  # Blue, Orange, Red

# 1) Count plot
ax1 = fig.add_subplot(gs[0, 0])
ax1.bar(time - width, num_hist, width, label=labels[0], color=colors[0], edgecolor='white', linewidth=1)
ax1.bar(time, num_gwl15, width, label=labels[1], color=colors[1], edgecolor='white', linewidth=1)
ax1.bar(time + width, num_gwl20, width, label=labels[2], color=colors[2], edgecolor='white', linewidth=1)
ax1.set_title('(a) MCS count', fontsize=font_size+2, loc='left', fontweight='bold')
ax1.set_ylabel('Mean MCS count', fontsize=font_size)
ax1.set_xticks(time)
ax1.set_xticklabels(xtick, fontsize=font_size)
ax1.grid(True, linestyle='--', alpha=0.5)
ax1.legend(loc='upper right', fontsize=font_size, ncol=3)

# 2) Size plot
ax2 = fig.add_subplot(gs[1, 0])
ax2.bar(time - width, size_hist, width, label=labels[0], color=colors[0], edgecolor='white', linewidth=1)
ax2.bar(time, size_gwl15, width, label=labels[1], color=colors[1], edgecolor='white', linewidth=1)
ax2.bar(time + width, size_gwl20, width, label=labels[2], color=colors[2], edgecolor='white', linewidth=1)
ax2.set_title('(b) MCS size', fontsize=font_size+2, loc='left', fontweight='bold')
ax2.set_ylabel('Size (10³ km²)', fontsize=font_size)
ax2.set_xticks(time)
ax2.set_xticklabels(xtick, fontsize=font_size)
ax2.grid(True, linestyle='--', alpha=0.5)

# 3) Total Precipitation plot
ax3 = fig.add_subplot(gs[2, 0])
ax3.bar(time - width, tot_hist, width, label=labels[0], color=colors[0], edgecolor='white', linewidth=1)
ax3.bar(time, tot_gwl15, width, label=labels[1], color=colors[1], edgecolor='white', linewidth=1)
ax3.bar(time + width, tot_gwl20, width, label=labels[2], color=colors[2], edgecolor='white', linewidth=1)
ax3.set_title('(c) MCS total precipitation', fontsize=font_size+2, loc='left', fontweight='bold')
ax3.set_ylabel('Total precipitation (mm)', fontsize=font_size)
ax3.set_xticks(time)
ax3.set_xticklabels(xtick, fontsize=font_size)
ax3.grid(True, linestyle='--', alpha=0.5)

# 4) Maximum Intensity plot
ax4 = fig.add_subplot(gs[3, 0])
ax4.bar(time - width, max_hist, width, label=labels[0], color=colors[0], edgecolor='white', linewidth=1)
ax4.bar(time, max_gwl15, width, label=labels[1], color=colors[1], edgecolor='white', linewidth=1)
ax4.bar(time + width, max_gwl20, width, label=labels[2], color=colors[2], edgecolor='white', linewidth=1)
ax4.set_title('(d) MCS max precipitation intensity', fontsize=font_size+2, loc='left', fontweight='bold')
ax4.set_xlabel('Month', fontsize=font_size)
ax4.set_ylabel('Max Precip (mm/h)', fontsize=font_size)
ax4.set_xticks(time)
ax4.set_xticklabels(xtick, fontsize=font_size)
ax4.grid(True, linestyle='--', alpha=0.5)

# Save output figure
path_out = "/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/figs/v2"
name_out = f"pyplt_graphs_moaap_mcs_ac_charac_{domain}_GWL_2000-2009.png"
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.show()
exit()
