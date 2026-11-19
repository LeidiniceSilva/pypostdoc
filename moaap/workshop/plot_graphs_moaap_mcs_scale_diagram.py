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


def select_mcs_in_box(mcs, lon_min, lon_max, lat_min, lat_max, mode="any"):
    """
    Keep only the MCSs inside the lon/lat box, using mcs[month][id]['track']
    (column 0 = lat, column 1 = lon).
    """
    out = {}

    for month, systems in mcs.items():
        kept = {}

        for mid, mcs_i in systems.items():
            track = np.asarray(mcs_i['track'])

            if track.ndim != 2 or track.shape[0] == 0:
                continue

            lat, lon = track[:, 0], track[:, 1]

            inside = ((lon >= lon_min) & (lon <= lon_max) &
                      (lat >= lat_min) & (lat <= lat_max))

            if mode == "all":
                ok = inside.all()
            elif mode == "any":
                ok = inside.any()
            elif mode == "genesis":
                ok = inside[0]
            else:
                raise ValueError("mode must be 'all', 'any' or 'genesis'")

            if ok:
                kept[mid] = mcs_i

        out[month] = kept

    return out


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


# Domain e Caixa de Seleção Espacial
domain = "EUR"
box = dict(lon_min=-10, lon_max=20, lat_min=38, lat_max=53)

# Load Data
mcs_eur_era5 = open_mcs("/ERA5/EURR-3/output", start="2000-01", end="2009-12")
mcs_eur_gpm = open_mcs("/GPM/EURR-3/output", start="2000-01", end="2009-12")
mcs_eur_cnrm_cpm_eval = open_mcs("/CPMs/CNRM/ALPX-3/evaluation/ERA5/output", start="2000-01", end="2009-12")
mcs_eur_cnrm_cpm_hist = open_mcs("/CPMs/CNRM/ALPX-3/historical/CNRM-ESM2-1/output", start="2000-01", end="2009-12")
mcs_eur_cnrm_rcm_eval = open_mcs("/RCMs/CNRM/EUR-12/evaluation/ERA5/output", start="2000-01", end="2009-12")
mcs_eur_cnrm_rcm_hist = open_mcs("/RCMs/CNRM/EUR-12/historical/CNRM-ESM2-1/output", start="2000-01", end="2009-12")
mcs_eur_ictp_rcm_eval = open_mcs("/RCMs/ICTP/EUR-12/evaluation/ERA5/output", start="2000-01", end="2009-12")
mcs_eur_ictp_rcm_hist = open_mcs("/RCMs/ICTP/EUR-12/historical/ECEarth/output", start="2000-01", end="2009-12")

# Aplicação da seleção espacial por caixa (mode="any")
mcs_eur_era5 = select_mcs_in_box(mcs_eur_era5, **box, mode="any")
mcs_eur_gpm = select_mcs_in_box(mcs_eur_gpm, **box, mode="any")
mcs_eur_cnrm_cpm_eval = select_mcs_in_box(mcs_eur_cnrm_cpm_eval, **box, mode="any")
mcs_eur_cnrm_cpm_hist = select_mcs_in_box(mcs_eur_cnrm_cpm_hist, **box, mode="any")
mcs_eur_cnrm_rcm_eval = select_mcs_in_box(mcs_eur_cnrm_rcm_eval, **box, mode="any")
mcs_eur_cnrm_rcm_hist = select_mcs_in_box(mcs_eur_cnrm_rcm_hist, **box, mode="any")
mcs_eur_ictp_rcm_eval = select_mcs_in_box(mcs_eur_ictp_rcm_eval, **box, mode="any")
mcs_eur_ictp_rcm_hist = select_mcs_in_box(mcs_eur_ictp_rcm_hist, **box, mode="any")

# Extract metrics for all 8 datasets
ls_era5, dur_era5 = process_dataset(mcs_eur_era5)
ls_gpm, dur_gpm = process_dataset(mcs_eur_gpm)
ls_cnrm_cpm_eval, dur_cnrm_cpm_eval = process_dataset(mcs_eur_cnrm_cpm_eval)
ls_cnrm_cpm_hist, dur_cnrm_cpm_hist = process_dataset(mcs_eur_cnrm_cpm_hist)
ls_cnrm_rcm_eval, dur_cnrm_rcm_eval = process_dataset(mcs_eur_cnrm_rcm_eval)
ls_cnrm_rcm_hist, dur_cnrm_rcm_hist = process_dataset(mcs_eur_cnrm_rcm_hist)
ls_ictp_rcm_eval, dur_ictp_rcm_eval = process_dataset(mcs_eur_ictp_rcm_eval)
ls_ictp_rcm_hist, dur_ictp_rcm_hist = process_dataset(mcs_eur_ictp_rcm_hist)

if len(ls_gpm) == 0:
    raise ValueError("No valid MCS track data extracted for GPM.")

# Setup figure grid
fig = plt.figure(figsize=(10, 10), dpi=400, facecolor='white')
gs = gridspec.GridSpec(3, 2, width_ratios=[4, 1], height_ratios=[0.9, 4, 0.6], hspace=0.08, wspace=0.08)
ax_top = fig.add_subplot(gs[0, 0], facecolor='white')
ax_main = fig.add_subplot(gs[1, 0], sharex=ax_top, facecolor='white')
ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main, facecolor='white')
ax_bottom = fig.add_subplot(gs[2, 0], sharex=ax_main, facecolor='white')

# Cores
color_era5 = "#222222"
color_gpm = "#6c6c6c"
color_cnrm_cpm = "#d1495b"
color_cnrm_rcm = "#2e86ab"
color_ictp_rcm = "#edae49"

# Main Panel: Log 2D Gaussian KDE
x_grid = np.linspace(1.5, 4.0, 100)  # ~30 km to 10000 km
y_grid = np.linspace(0.5, 2.5, 100)  # ~3 hours to 300 hours
X_log, Y_log = np.meshgrid(x_grid, y_grid)
X = 10**X_log
Y = 10**Y_log

def plot_kde(ax, ls, dur, color, linestyle="-"):
    if len(ls) < 2 or len(dur) < 2:
        return
        
    log_x = np.log10(ls)
    log_y = np.log10(dur)
    
    kde = gaussian_kde(np.vstack([log_x, log_y]), bw_method=0.4)
    Z = kde(np.vstack([X_log.ravel(), Y_log.ravel()])).reshape(X_log.shape)
    
    mask = (X_log >= np.min(log_x)) & (X_log <= np.max(log_x)) & \
           (Y_log >= np.min(log_y)) & (Y_log <= np.max(log_y))
    Z_masked = np.where(mask, Z, 0)
    
    if np.any(Z_masked > 0):
        z_min = np.percentile(Z_masked[Z_masked > 0], 80)
        levels = np.linspace(z_min, Z_masked.max(), 5)
        ax.contour(X, Y, Z_masked, levels=levels, colors=[color], linewidths=1.5, linestyles=linestyle)

plot_kde(ax_main, ls_era5, dur_era5, color_era5, "-")
plot_kde(ax_main, ls_gpm, dur_gpm, color_gpm, "-")
plot_kde(ax_main, ls_cnrm_cpm_eval, dur_cnrm_cpm_eval, color_cnrm_cpm, "-")
plot_kde(ax_main, ls_cnrm_cpm_hist, dur_cnrm_cpm_hist, color_cnrm_cpm, "--")
plot_kde(ax_main, ls_cnrm_rcm_eval, dur_cnrm_rcm_eval, color_cnrm_rcm, "-")
plot_kde(ax_main, ls_cnrm_rcm_hist, dur_cnrm_rcm_hist, color_cnrm_rcm, "--")
plot_kde(ax_main, ls_ictp_rcm_eval, dur_ictp_rcm_eval, color_ictp_rcm, "-")
plot_kde(ax_main, ls_ictp_rcm_hist, dur_ictp_rcm_hist, color_ictp_rcm, "--")

ax_main.set_xscale("log")
ax_main.set_yscale("log")
ax_main.set_xlim(20, 10000)
ax_main.set_ylim(3, 300)

# Add Inset Domain Map on Top-Left
ax_inset = inset_axes(ax_main, width="30%", height="30%", loc="upper left", bbox_to_anchor=(0.03, -0.03, 1, 1), bbox_transform=ax_main.transAxes, axes_class=GeoAxes, axes_kwargs=dict(projection=ccrs.PlateCarree()))

ax_inset.set_extent([box["lon_min"], box["lon_max"], box["lat_min"], box["lat_max"]], crs=ccrs.PlateCarree())
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
def plot_top_boxplot(ax, data, color, y_pos, hatch=""):
    if len(data) == 0:
        return
    q5, q25, q50, q75, q95 = np.percentile(data, [5, 25, 50, 75, 95])
    min_val, max_val = np.min(data), np.max(data)
    ax.hlines(y_pos, q5, q95, color=color, linewidth=1.2)
    ax.barh(y_pos, q75 - q25, left=q25, height=0.08, color=color, edgecolor=color, alpha=0.7, hatch=hatch)
    ax.plot([min_val, max_val], [y_pos, y_pos], "o", color=color, markersize=2.5)
    ax.plot(q50, y_pos, "o", color="white", markersize=2.5)

y_positions = np.linspace(0.7, -0.7, 8)

plot_top_boxplot(ax_top, ls_era5, color_era5, y_positions[0], "")
plot_top_boxplot(ax_top, ls_gpm, color_gpm, y_positions[1], "")
plot_top_boxplot(ax_top, ls_cnrm_cpm_eval, color_cnrm_cpm, y_positions[2], "")
plot_top_boxplot(ax_top, ls_cnrm_cpm_hist, color_cnrm_cpm, y_positions[3], "///")
plot_top_boxplot(ax_top, ls_cnrm_rcm_eval, color_cnrm_rcm, y_positions[4], "")
plot_top_boxplot(ax_top, ls_cnrm_rcm_hist, color_cnrm_rcm, y_positions[5], "///")
plot_top_boxplot(ax_top, ls_ictp_rcm_eval, color_ictp_rcm, y_positions[6], "")
plot_top_boxplot(ax_top, ls_ictp_rcm_hist, color_ictp_rcm, y_positions[7], "///")

labels = [
    ("ERA5", color_era5),
    ("GPM", color_gpm),
    ("CNRM CPM-3 Eval", color_cnrm_cpm),
    ("CNRM CPM-3 Hist", color_cnrm_cpm),
    ("CNRM RCM-12 Eval", color_cnrm_rcm),
    ("CNRM RCM-12 Hist", color_cnrm_rcm),
    ("ICTP RCM-12 Eval", color_ictp_rcm),
    ("ICTP RCM-12 Hist", color_ictp_rcm),
]

text_y = np.linspace(0.92, 0.08, 8)
for idx, (lbl, col) in enumerate(labels):
    ax_top.text(0.02, text_y[idx], lbl, transform=ax_top.transAxes, color=col, fontsize=7.5, weight="bold", va="center")

ax_top.set_ylim(-0.9, 0.9)
ax_top.axis("off")

# Right Boxplot (Duration)
def plot_right_boxplot(ax, data, color, x_pos, hatch=""):
    if len(data) == 0:
        return
    q5, q25, q50, q75, q95 = np.percentile(data, [5, 25, 50, 75, 95])
    min_val, max_val = np.min(data), np.max(data)
    ax.vlines(x_pos, q5, q95, color=color, linewidth=1.2)
    ax.bar(x_pos, q75 - q25, bottom=q25, width=0.08, color=color, edgecolor=color, alpha=0.7, hatch=hatch)
    ax.plot([x_pos, x_pos], [min_val, max_val], "o", color=color, markersize=2.5)
    ax.plot(x_pos, q50, "o", color="white", markersize=2.5)

x_positions = np.linspace(-0.7, 0.7, 8)

plot_right_boxplot(ax_right, dur_era5, color_era5, x_positions[0], "")
plot_right_boxplot(ax_right, dur_gpm, color_gpm, x_positions[1], "")
plot_right_boxplot(ax_right, dur_cnrm_cpm_eval, color_cnrm_cpm, x_positions[2], "")
plot_right_boxplot(ax_right, dur_cnrm_cpm_hist, color_cnrm_cpm, x_positions[3], "///")
plot_right_boxplot(ax_right, dur_cnrm_rcm_eval, color_cnrm_rcm, x_positions[4], "")
plot_right_boxplot(ax_right, dur_cnrm_rcm_hist, color_cnrm_rcm, x_positions[5], "///")
plot_right_boxplot(ax_right, dur_ictp_rcm_eval, color_ictp_rcm, x_positions[6], "")
plot_right_boxplot(ax_right, dur_ictp_rcm_hist, color_ictp_rcm, x_positions[7], "///")

ax_right.set_xlim(-0.9, 0.9)
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
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/figs/workshop'
name_out = f'pyplt_graphs_moaap_mcs_scale_diagram_{domain}_2000-2009.png'

os.makedirs(path_out, exist_ok=True)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.show()
