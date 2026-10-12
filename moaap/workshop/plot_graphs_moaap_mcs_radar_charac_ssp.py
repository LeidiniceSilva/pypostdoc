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


def compute_mcs_metrics(mcs_dict):
    """Computes MCS median metrics. Returns (None, None, 0) if dictionary is empty or no MCSs exist."""
    if not mcs_dict:
        return None, None, 0

    # Flatten monthly dictionaries
    all_mcs = []
    for month, systems in mcs_dict.items():
        if isinstance(systems, dict):
            all_mcs.extend(systems.values())
        elif isinstance(systems, list):
            all_mcs.extend(systems)

    n_systems = len(all_mcs)
    if n_systems == 0:
        return None, None, 0

    # Data structures to collect MCS properties
    durations, distances, speeds = [], [], []
    mean_areas, volumes, hp_volumes, severities = [], [], [], []
    pr_mean, pr_p10, pr_p25, pr_median, pr_p75, pr_p90, pr_p99, pr_max = (
        [], [], [], [], [], [], [], []
    )

    for mcs in all_mcs:
        # Duration (hours)
        dur = len(mcs.get("times", []))
        durations.append(dur)

        # Distance Traveled (km)
        track = np.array(mcs.get("track", []))
        if len(track) > 1:
            dists = np.sqrt(np.diff(track[:, 0]) ** 2 + np.diff(track[:, 1]) ** 2) * 111.0
            dist_total = np.sum(dists)
        else:
            dist_total = 0.0
        distances.append(dist_total)

        # Speed (km/h)
        spd = mcs.get("speed", [])
        speed_mean = np.mean(spd) if len(spd) > 0 else 0.0
        speeds.append(speed_mean)

        # Spatial / Volumetric characteristics
        sizes = np.array(mcs.get("size", [0]))
        means = np.array(mcs.get("mean", [0]))
        maxs = np.array(mcs.get("max", [0]))
        tots = np.array(mcs.get("tot", [0]))

        m_area = np.mean(sizes) / 1e6  # Convert m^2 to km^2
        mean_areas.append(m_area)

        vol = np.sum(sizes * means) / 1e6  # Volume
        volumes.append(vol)

        hp_vol = np.sum(sizes * maxs)  # High-precipitation volume
        hp_volumes.append(hp_vol)

        sev = np.mean(tots)  # Severity
        severities.append(sev)

        # Intensity / Precipitation metrics
        pr_mean.append(np.mean(means))
        pr_p10.append(np.percentile(means, 10))
        pr_p25.append(np.percentile(means, 25))
        pr_median.append(np.median(means))
        pr_p75.append(np.percentile(means, 75))
        pr_p90.append(np.percentile(means, 90))
        pr_p99.append(np.percentile(means, 99))
        pr_max.append(np.max(maxs))

    medians = {
        "Duration": np.median(durations),
        "DistanceTraveled": np.median(distances),
        "Severity": np.median(severities),
        "Speed": np.median(speeds),
        "max(Pr)": np.median(pr_max),
        "P99(Pr)": np.median(pr_p99),
        "P90(Pr)": np.median(pr_p90),
        "P75(Pr)": np.median(pr_p75),
        "Median(Pr)": np.median(pr_median),
        "P25(Pr)": np.median(pr_p25),
        "P10(Pr)": np.median(pr_p10),
        "Mean(Pr)": np.median(pr_mean),
        "HPVolume": np.median(hp_volumes),
        "Volume": np.median(volumes),
        "MeanArea": np.median(mean_areas),
    }

    raw_metrics = {
        "durations": durations,
        "distances": distances,
        "severities": severities,
        "speeds": speeds,
        "pr_max": pr_max,
        "pr_p99": pr_p99,
        "pr_p90": pr_p90,
        "pr_p75": pr_p75,
        "pr_median": pr_median,
        "pr_p25": pr_p25,
        "pr_p10": pr_p10,
        "pr_mean": pr_mean,
        "hp_volumes": hp_volumes,
        "volumes": volumes,
        "mean_areas": mean_areas,
    }

    return medians, raw_metrics, n_systems


# Domain
domain = "EUR"

# Path definitions
rcm_hist_path = "RCMs/ICTP/EUR-12/historical/ECEarth/output"
rcm_ssp_path  = "RCMs/ICTP/EUR-12/ssp370/ECEarth/output"

cpm_hist_path = "CPMs/ICTP/EURR-3/historical/ECEarth/output"
cpm_ssp_path  = "CPMs/ICTP/EURR-3/ssp370/ECEarth/output"

# 1. Historical (2000-01 to 2009-12)
mcs_rcm_hist = open_mcs(rcm_hist_path, start="2000-01", end="2009-12")

# 2. GWL 1.5 (Part in historical 2000-01 to 2009-12 and part in ssp370 2015-01 to 2021-12)
mcs_rcm_gwl15_part1 = open_mcs(rcm_hist_path, start="2000-01", end="2009-12")
mcs_rcm_gwl15_part2 = open_mcs(rcm_ssp_path, start="2015-01", end="2021-12")
mcs_rcm_gwl15 = {**mcs_rcm_gwl15_part1, **mcs_rcm_gwl15_part2}

# 3. GWL 2.0 (2023-01 to 2042-12)
mcs_rcm_gwl20 = open_mcs(rcm_ssp_path, start="2023-01", end="2042-12")

# Check CPM-3 safely (returns 0 systems if missing)
mcs_cpm_hist = open_mcs(cpm_hist_path, start="2000-01", end="2009-12")
medians_cpm_hist, _, n_cpm_hist = compute_mcs_metrics(mcs_cpm_hist)

# Compute metrics for RCM datasets
medians_hist, _, n_hist = compute_mcs_metrics(mcs_rcm_hist)
medians_gwl15, _, n_gwl15 = compute_mcs_metrics(mcs_rcm_gwl15)
medians_gwl20, _, n_gwl20 = compute_mcs_metrics(mcs_rcm_gwl20)

# Unit formatting dictionary
units = {
    "Duration": "h",
    "DistanceTraveled": "km",
    "Severity": "E3 m^3",
    "Speed": "km h^{-1}",
    "max(Pr)": "mm h^{-1}",
    "P99(Pr)": "mm h^{-1}",
    "P90(Pr)": "mm h^{-1}",
    "P75(Pr)": "mm h^{-1}",
    "Median(Pr)": "mm h^{-1}",
    "P25(Pr)": "mm h^{-1}",
    "P10(Pr)": "mm h^{-1}",
    "Mean(Pr)": "mm h^{-1}",
    "HPVolume": "E8 m^3",
    "Volume": "E4 km^2 h",
    "MeanArea": "E3 km^2",
}

if medians_hist is None:
    raise RuntimeError("No historical MCS data available to compute radar limits.")

labels = list(medians_hist.keys())
num_vars = len(labels)

values_hist = np.array(list(medians_hist.values()))
values_gwl15 = np.array(list(medians_gwl15.values())) if medians_gwl15 else np.zeros(num_vars)
values_gwl20 = np.array(list(medians_gwl20.values())) if medians_gwl20 else np.zeros(num_vars)

# Set normalization axis baseline from historical
max_axis_limits = values_hist * 1.5

norm_hist = np.concatenate((values_hist / max_axis_limits, [(values_hist / max_axis_limits)[0]]))
norm_gwl15 = np.concatenate((values_gwl15 / max_axis_limits, [(values_gwl15 / max_axis_limits)[0]]))
norm_gwl20 = np.concatenate((values_gwl20 / max_axis_limits, [(values_gwl20 / max_axis_limits)[0]]))

angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
angles += angles[:1]

# Colors exactly as shown in the reference image
c_hist = "#1f77b4"   # Blue
c_gwl15 = "#ff7f0e"  # Orange / Nearfuture
c_gwl20 = "#d62728"  # Red / Farfuture

# Create Plot
fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))

if medians_gwl15:
    ax.plot(angles, norm_gwl15, color=c_gwl15, linewidth=1.0, marker="o", markersize=3, label="gwl15")
    ax.fill(angles, norm_gwl15, color=c_gwl15, alpha=0.25)

if medians_gwl20:
    ax.plot(angles, norm_gwl20, color=c_gwl20, linewidth=1.0, marker="o", markersize=3, label="gwl20")
    ax.fill(angles, norm_gwl20, color=c_gwl20, alpha=0.25)

# Single line plots and fills for each dataset/period
ax.plot(angles, norm_hist, color=c_hist, linewidth=1.0, marker="o", markersize=3, label="historical")
ax.fill(angles, norm_hist, color=c_hist, alpha=0.25)

# Radar polar axis settings
ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)
ax.set_xticks(angles[:-1])
ax.set_yticks(np.linspace(0.05, 1.0, 20))
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_ylim(0, 1)
ax.grid(True, color="black", linestyle="-", linewidth=0.3, alpha=0.4)

# Perimeter metric labels with M_hist values
for i, angle in enumerate(angles[:-1]):
    lbl = labels[i]
    val = values_hist[i]
    unit = units.get(lbl, "")
    val_str = f"{val/1000:.1f}" if val >= 1000 else f"{val:.1f}"

    text_label = f"${lbl}$\n$M_{{hist}} = {val_str}\\ {unit}$"
    ha = "left" if np.cos(angle) >= 0 else "right"
    ax.text(angle, 1.15, text_label, size=9, horizontalalignment=ha, verticalalignment="center")

# Annotations matched to the reference layout
ax.text(-0.20, 0.99, "Allyear\nAllover", transform=ax.transAxes, color="black", fontsize=12, fontweight="bold")
ax.text(-0.20, 0.93, f"historical: N={n_hist}", transform=ax.transAxes, color=c_hist, fontsize=10, fontweight="bold")
ax.text(-0.20, 0.89, f"GWL15: N={n_gwl15}", transform=ax.transAxes, color=c_gwl15, fontsize=10, fontweight="bold")
ax.text(-0.20, 0.85, f"GWL20: N={n_gwl20}", transform=ax.transAxes, color=c_gwl20, fontsize=10, fontweight="bold")

# Save figure
path_out = "/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/figs/v2"
name_out = f"pyplt_graphs_moaap_mcs_radar_charac_{domain}_GWL_2000-2009.png"
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches="tight", facecolor="white", edgecolor="none")
plt.show()
