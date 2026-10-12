# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 03, 2026"
__description__ = "This script plot MCSs"

import os
import pickle
import warnings
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cartopy.mpl.geoaxes import GeoAxes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import gaussian_kde

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


def process_dataset(mcs_dict):
    length_scales = []
    durations = []
    for month_key, obj_dict in mcs_dict.items():
        for obj_id, data in obj_dict.items():
            if "times" in data and "size" in data:
                times = data["times"]
                sizes = data["size"]

                if len(times) < 2:
                    continue

                # Calculate lifetime in hours
                duration_hrs = (times[-1] - times[0]).total_seconds() / 3600.0

                # Compute characteristic horizontal scale (km)
                mean_area_m2 = np.mean(sizes)
                length_scale_km = np.sqrt(mean_area_m2) / 1000.0

                if duration_hrs >= 4 and length_scale_km >= 100:
                    length_scales.append(length_scale_km)
                    durations.append(duration_hrs)

    return np.array(length_scales), np.array(durations)


# Load Data
mcs_eur_gpm = open_mcs("/GPM/EURR-3/output", start="2000-01", end="2009-12")
mcs_eur_cpm_eval = open_mcs("/CPMs/ICTP/EURR-3/evaluation/ERA5/output", start="2000-01", end="2009-12")
mcs_eur_cpm_hist = open_mcs("/CPMs/ICTP/EURR-3/historical/ECEarth/output", start="2000-01", end="2009-12")
mcs_eur_rcm_eval = open_mcs("/RCMs/ICTP/EUR-12/evaluation/ERA5/output", start="2000-01", end="2009-12")
mcs_eur_rcm_hist = open_mcs("/RCMs/ICTP/EUR-12/historical/ECEarth/output", start="2000-01", end="2009-12")

# Extract metrics
ls_gpm, dur_gpm = process_dataset(mcs_eur_gpm)
ls_cpm_eval, dur_cpm_eval = process_dataset(mcs_eur_cpm_eval)
ls_cpm_hist, dur_cpm_hist = process_dataset(mcs_eur_cpm_hist)
ls_rcm_eval, dur_rcm_eval = process_dataset(mcs_eur_rcm_eval)
ls_rcm_hist, dur_rcm_hist = process_dataset(mcs_eur_rcm_hist)

if len(ls_gpm) == 0 or len(ls_cpm_eval) == 0 or len(ls_cpm_hist) == 0:
    raise ValueError("No valid MCS track data extracted.")

# Setup figure grid
fig = plt.figure(figsize=(10, 8), dpi=400, facecolor='white')
gs = gridspec.GridSpec(3, 2, width_ratios=[4, 1], height_ratios=[0.8, 4, 0.6], hspace=0.08, wspace=0.08)
ax_top = fig.add_subplot(gs[0, 0], facecolor='white')
ax_main = fig.add_subplot(gs[1, 0], sharex=ax_top, facecolor='white')
ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main, facecolor='white')
ax_bottom = fig.add_subplot(gs[2, 0], sharex=ax_main, facecolor='white')

color_gpm = "black"
color_cpm_eval = "red"
color_cpm_hist = "blue"
color_rcm_eval = "green"
color_rcm_hist = "orange"

# Main Panel: Log 2D Gaussian KDE
x_grid = np.linspace(1.5, 4.0, 100)  # ~30 km to 10000 km
y_grid = np.linspace(0.5, 2.5, 100)  # ~3 hours to 300 hours
X_log, Y_log = np.meshgrid(x_grid, y_grid)
X = 10**X_log
Y = 10**Y_log

def plot_kde(ax, ls, dur, color):
    if len(ls) < 2 or len(dur) < 2:
        return
        
    log_x = np.log10(ls)
    log_y = np.log10(dur)
    
    kde = gaussian_kde(np.vstack([log_x, log_y]), bw_method=0.4)
    Z = kde(np.vstack([X_log.ravel(), Y_log.ravel()])).reshape(X_log.shape)
    
    # Mask out region outside data range to remove edge boundary artifacts
    mask = (X_log >= np.min(log_x)) & (X_log <= np.max(log_x)) & \
           (Y_log >= np.min(log_y)) & (Y_log <= np.max(log_y))
    Z_masked = np.where(mask, Z, 0)
    
    if np.any(Z_masked > 0):
        z_min = np.percentile(Z_masked[Z_masked > 0], 80)
        levels = np.linspace(z_min, Z_masked.max(), 5)
        ax.contour(X, Y, Z_masked, levels=levels, colors=[color], linewidths=1.5)

plot_kde(ax_main, ls_gpm, dur_gpm, color_gpm)
plot_kde(ax_main, ls_cpm_eval, dur_cpm_eval, color_cpm_eval)
plot_kde(ax_main, ls_cpm_hist, dur_cpm_hist, color_cpm_hist)
plot_kde(ax_main, ls_rcm_eval, dur_rcm_eval, color_rcm_eval)
plot_kde(ax_main, ls_rcm_hist, dur_rcm_hist, color_rcm_hist)

ax_main.set_xscale("log")
ax_main.set_yscale("log")
ax_main.set_xlim(20, 10000)
ax_main.set_ylim(3, 300)

# Add Inset Domain Map on Top-Left
ax_inset = inset_axes(ax_main, width="30%", height="30%", loc="upper left", bbox_to_anchor=(0.03, -0.03, 1, 1), bbox_transform=ax_main.transAxes, axes_class=GeoAxes, axes_kwargs=dict(projection=ccrs.PlateCarree()))

# Set EURR-3 Domain: Lat (-17, 36), Lon (36, 58)
ax_inset.set_extent([-17, 36, 36, 58], crs=ccrs.PlateCarree())
ax_inset.add_feature(cfeature.OCEAN, facecolor="lightblue")
ax_inset.add_feature(cfeature.LAND, facecolor="tan", edgecolor="black", linewidth=0.5)
ax_inset.add_feature(cfeature.BORDERS, linestyle=":", linewidth=0.4, edgecolor="black")
ax_inset.add_feature(cfeature.COASTLINE, linewidth=0.5)

# Custom Y-axis Ticks
y_ticks = [6, 24, 240]
y_labels = ["6 hours", "1 day", "10 days"]
ax_main.set_yticks(y_ticks)
ax_main.set_yticklabels(y_labels)

# Background scale lines
ax_main.axhline(6, color="gray", linewidth=0.8, zorder=1)
ax_main.axhline(24, color="gray", linewidth=0.8, zorder=1)
ax_main.axhline(240, color="gray", linewidth=0.8, zorder=1)
ax_main.axvline(200, color="black", linewidth=1.0, zorder=1)
ax_main.axvline(2000, color="black", linewidth=1.0, zorder=1)

ax_main.set_xlabel("Horizontal length scale [km]", fontsize=12)
ax_main.set_ylabel("Time [hours]", fontsize=12)

# Top Boxplot (Length Scale) 
def plot_top_boxplot(ax, data, color, y_pos):
    if len(data) == 0:
        return
    q5, q25, q50, q75, q95 = np.percentile(data, [5, 25, 50, 75, 95])
    min_val, max_val = np.min(data), np.max(data)
    ax.hlines(y_pos, q5, q95, color=color, linewidth=1.2)
    ax.barh(y_pos, q75 - q25, left=q25, height=0.1, color=color, edgecolor=color, alpha=0.7)
    ax.plot([min_val, max_val], [y_pos, y_pos], "o", color=color, markersize=3)
    ax.plot(q50, y_pos, "o", color="white", markersize=3)

# datasets positioned evenly across vertical span (-0.6 to 0.6)
plot_top_boxplot(ax_top, ls_gpm, color_gpm, y_pos=0.6)
plot_top_boxplot(ax_top, ls_cpm_eval, color_cpm_eval, y_pos=0.3)
plot_top_boxplot(ax_top, ls_cpm_hist, color_cpm_hist, y_pos=0.0)
plot_top_boxplot(ax_top, ls_rcm_eval, color_rcm_eval, y_pos=-0.3)
plot_top_boxplot(ax_top, ls_rcm_hist, color_rcm_hist, y_pos=-0.6)

# Distinct label offsets matching the boxplot y-positions
ax_top.text(0.02, 0.90, "GPM", transform=ax_top.transAxes, color=color_gpm, fontsize=8, weight="bold", va="center")
ax_top.text(0.02, 0.70, "CPM-3 Eval", transform=ax_top.transAxes, color=color_cpm_eval, fontsize=8, weight="bold", va="center")
ax_top.text(0.02, 0.50, "CPM-3 Hist", transform=ax_top.transAxes, color=color_cpm_hist, fontsize=8, weight="bold", va="center")
ax_top.text(0.02, 0.30, "RCM-12 Eval", transform=ax_top.transAxes, color=color_rcm_eval, fontsize=8, weight="bold", va="center")
ax_top.text(0.02, 0.10, "RCM-12 Hist", transform=ax_top.transAxes, color=color_rcm_hist, fontsize=8, weight="bold", va="center")

ax_top.set_ylim(-0.8, 0.8)
ax_top.axis("off")

# Right Boxplot (Duration)
def plot_right_boxplot(ax, data, color, x_pos):
    if len(data) == 0:
        return
    q5, q25, q50, q75, q95 = np.percentile(data, [5, 25, 50, 75, 95])
    min_val, max_val = np.min(data), np.max(data)
    ax.vlines(x_pos, q5, q95, color=color, linewidth=1.2)
    ax.bar(x_pos, q75 - q25, bottom=q25, width=0.1, color=color, edgecolor=color, alpha=0.7)
    ax.plot([x_pos, x_pos], [min_val, max_val], "o", color=color, markersize=3)
    ax.plot(x_pos, q50, "o", color="white", markersize=3)

# Datasets positioned evenly across horizontal span (-0.6 to 0.6)
plot_right_boxplot(ax_right, dur_gpm, color_gpm, x_pos=-0.6)
plot_right_boxplot(ax_right, dur_cpm_eval, color_cpm_eval, x_pos=-0.3)
plot_right_boxplot(ax_right, dur_cpm_hist, color_cpm_hist, x_pos=0.0)
plot_right_boxplot(ax_right, dur_rcm_eval, color_rcm_eval, x_pos=0.3)
plot_right_boxplot(ax_right, dur_rcm_hist, color_rcm_hist, x_pos=0.6)

ax_right.set_xlim(-0.8, 0.8)
ax_right.axis("off")

# Bottom Annotations 
ax_bottom.axis("off")
ax_bottom.set_ylim(-0.8, 1.0)

ax_bottom.axvline(200, ymin=-0.5, ymax=0.4, color="black", linewidth=1.5, clip_on=False)
ax_bottom.axvline(2000, ymin=-0.5, ymax=0.4, color="black", linewidth=1.5, clip_on=False)

ax_bottom.text(60, -0.6, "Meso\nalpha", ha="center", va="center", fontsize=10, weight="bold")
ax_bottom.text(600, -0.6, "Macro\nbeta", ha="center", va="center", fontsize=10, weight="bold")
ax_bottom.text(4000, -0.6, "Macro\nalpha", ha="center", va="center", fontsize=10, weight="bold")

# Save figure 
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/figs/v2'
name_out = 'pyplt_graphs_moaap_mcs_scale_diagram_EUR_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.show()
exit()
