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


def comp_lifetime_cycle(mcs_charac):

    # Dicionários para armazenar características por hora de vida
    size_by_hour = {h: [] for h in range(5, 21)}
    tot_by_hour = {h: [] for h in range(5, 21)}
    max_by_hour = {h: [] for h in range(5, 21)}

    for obj in mcs_charac.keys():

        for m in mcs_charac[obj].keys():

            times = mcs_charac[obj][m]['times']
            lifetime = len(times)
            
            # Selecionar MCSs com lifetime entre 5 e 20 horas
            if lifetime >= 5 and lifetime <= 20:
                
                # Para cada hora de vida do MCS (1 = primeira hora, 2 = segunda, etc.)
                for hour_idx in range(lifetime):
                    lifetime_hour = hour_idx + 1  # 1a hora, 2a hora, etc.
                    
                    # Verificar se esta hora de vida está no nosso range (5-20)
                    if lifetime_hour >= 5 and lifetime_hour <= 20:
                        
                        size = mcs_charac[obj][m]['size'][hour_idx] / 1000**2
                        tot = mcs_charac[obj][m]['tot'][hour_idx]
                        maxv = mcs_charac[obj][m]['max'][hour_idx]
                        
                        size_by_hour[lifetime_hour].append(size)
                        tot_by_hour[lifetime_hour].append(tot)
                        max_by_hour[lifetime_hour].append(maxv)

    # Calcular médias para cada hora de vida
    size_clim = []
    tot_clim = []
    max_clim = []
    hours = []
    
    for h in range(5, 21):
        if len(size_by_hour[h]) > 0:
            size_clim.append(np.nanmean(size_by_hour[h]))
            tot_clim.append(np.nanmean(tot_by_hour[h]))
            max_clim.append(np.nanmean(max_by_hour[h]))
            hours.append(h)
        else:
            size_clim.append(np.nan)
            tot_clim.append(np.nan)
            max_clim.append(np.nan)
            hours.append(h)

    return np.array(size_clim), np.array(tot_clim), np.array(max_clim), np.array(hours)


# Domain
domain = "EUR"

# Load datasets
mcs_eur_gpm = open_mcs("/GPM/EURR-3/output", start="2000-01", end="2009-12")
mcs_eur_cpm_eval = open_mcs("/CPMs/ICTP/EURR-3/evaluation/ERA5/output", start="2000-01", end="2009-12")
mcs_eur_cpm_hist = open_mcs("/CPMs/ICTP/EURR-3/historical/ECEarth/output", start="2000-01", end="2009-12")
mcs_eur_rcm_eval = open_mcs("/RCMs/ICTP/EUR-12/evaluation/ERA5/output", start="2000-01", end="2009-12")
mcs_eur_rcm_hist = open_mcs("/RCMs/ICTP/EUR-12/historical/ECEarth/output", start="2000-01", end="2009-12")

# Calculate metrics for each dataset
size_gpm, tot_gpm, max_gpm, hours = comp_lifetime_cycle(mcs_eur_gpm)
size_cpm_eval, tot_cpm_eval, max_cpm_eval, _ = comp_lifetime_cycle(mcs_eur_cpm_eval)
size_cpm_hist, tot_cpm_hist, max_cpm_hist, _ = comp_lifetime_cycle(mcs_eur_cpm_hist)
size_rcm_eval, tot_rcm_eval, max_rcm_eval, _ = comp_lifetime_cycle(mcs_eur_rcm_eval)
size_rcm_hist, tot_rcm_hist, max_rcm_hist, _ = comp_lifetime_cycle(mcs_eur_rcm_hist)

# Plot parameters
fig = plt.figure(figsize=(18, 10))
font_size = 10
gs = gridspec.GridSpec(3, 1, figure=fig, hspace=0.35)

x = np.arange(len(hours))  # bar locations along x-axis
width = 0.15  # width of each bar

labels = ['GPM', 'CPM-3 Eval', 'CPM-3 Hist', 'RCM-12 Eval', 'RCM-12 Hist']
colors = ['black', 'red', 'blue', 'green', 'orange']

# Size plot
ax1 = fig.add_subplot(gs[0, 0])
ax1.bar(x - 2*width, size_gpm, width, label=labels[0], color=colors[0], edgecolor='white', alpha=0.75)
ax1.bar(x - width, size_cpm_eval, width, label=labels[1], color=colors[1], edgecolor='white', alpha=0.75)
ax1.bar(x, size_cpm_hist, width, label=labels[2], color=colors[2], edgecolor='white', alpha=0.75)
ax1.bar(x + width, size_rcm_eval, width, label=labels[3], color=colors[3], edgecolor='white', alpha=0.75)
ax1.bar(x + 2*width, size_rcm_hist, width, label=labels[4], color=colors[4], edgecolor='white', alpha=0.75)
ax1.set_title('(a) MCS size', fontsize=font_size+2, loc='left', fontweight='bold', alpha=0.75)
ax1.set_ylabel('Size (10³ km²)', fontsize=font_size)
ax1.set_ylim(0, 200000)  
ax1.set_xticks(x)
ax1.set_xticklabels(hours, fontsize=font_size)
ax1.grid(True, linestyle='--', alpha=0.5)

# Total Precipitation plot
ax2 = fig.add_subplot(gs[1, 0])
ax2.bar(x - 2*width, tot_gpm, width, label=labels[0], color=colors[0], edgecolor='white', alpha=0.75)
ax2.bar(x - width, tot_cpm_eval, width, label=labels[1], color=colors[1], edgecolor='white', alpha=0.75)
ax2.bar(x, tot_cpm_hist, width, label=labels[2], color=colors[2], edgecolor='white', alpha=0.75)
ax2.bar(x + width, tot_rcm_eval, width, label=labels[3], color=colors[3], edgecolor='white', alpha=0.75)
ax2.bar(x + 2*width, tot_rcm_hist, width, label=labels[4], color=colors[4], edgecolor='white', alpha=0.75)
ax2.set_title('(b) MCS total precipitation', fontsize=font_size+2, loc='left', fontweight='bold', alpha=0.75)
ax2.set_ylabel('Total precipitation (mm)', fontsize=font_size)
ax2.set_ylim(0, 1000)  
ax2.set_xticks(x)
ax2.set_xticklabels(hours, fontsize=font_size)
ax2.grid(True, linestyle='--', alpha=0.5)

# Maximum Intensity plot
ax3 = fig.add_subplot(gs[2, 0])
ax3.bar(x - 2*width, max_gpm, width, label=labels[0], color=colors[0], edgecolor='white', alpha=0.75)
ax3.bar(x - width, max_cpm_eval, width, label=labels[1], color=colors[1], edgecolor='white', alpha=0.75)
ax3.bar(x, max_cpm_hist, width, label=labels[2], color=colors[2], edgecolor='white', alpha=0.75)
ax3.bar(x + width, max_rcm_eval, width, label=labels[3], color=colors[3], edgecolor='white', alpha=0.75)
ax3.bar(x + 2*width, max_rcm_hist, width, label=labels[4], color=colors[4], edgecolor='white', alpha=0.75)
ax3.set_title('(c) MCS max precipitation intensity', fontsize=font_size+2, loc='left', fontweight='bold', alpha=0.75)
ax3.set_xlabel('Lifetime Hour (h)', fontsize=font_size)
ax3.set_ylabel('Max precipitation (mm/h)', fontsize=font_size)
ax3.set_ylim(0, 25)  
ax3.set_xticks(x)
ax3.set_xticklabels(hours, fontsize=font_size)
ax3.grid(True, linestyle='--', alpha=0.5)
ax3.legend(loc='upper right', fontsize=font_size, ncol=5)

# Path out to save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/figs/v2'
name_out = f'pyplt_graphs_moaap_mcs_lifetime_charac_{domain}_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
