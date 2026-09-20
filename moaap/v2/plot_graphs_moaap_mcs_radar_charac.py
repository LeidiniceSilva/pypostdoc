# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 03, 2026"
__description__ = "This script plot MCSs"

import os
import pickle
import warnings
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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


def compute_mcs_metrics(mcs_dict):

    # Flatten monthly dictionaries
    all_mcs = []
    for month, systems in mcs_dict.items():
        all_mcs.extend(systems.values())

    n_systems = len(all_mcs)
    if n_systems == 0:
        raise ValueError("No MCS found in the dict")

    # Data structures to collect MCS
    durations, distances, speeds = [], [], []
    mean_areas, volumes, hp_volumes, severities = [], [], [], []
    (pr_mean, pr_p10, pr_p25, pr_median, pr_p75, pr_p90, pr_p99, pr_max,) = ([], [], [], [], [], [], [], [],)

    for mcs in all_mcs:
        # Duration (hours)
        dur = len(mcs["times"])
        durations.append(dur)

        # Distance Traveled (km)
        track = mcs["track"]  # [[lat, lon], ...]
        if len(track) > 1:
            dists = ( np.sqrt(np.diff(track[:, 0]) ** 2 + np.diff(track[:, 1]) ** 2) * 111.0)
            dist_total = np.sum(dists)
        else:
            dist_total = 0.0
        distances.append(dist_total)

        # Speed (km/h)
        speed_mean = np.mean(mcs["speed"]) if len(mcs["speed"]) > 0 else 0.0
        speeds.append(speed_mean)

        # Spatial / Volumetric characteristics
        m_area = np.mean(mcs["size"]) / 1e6  # Convert m^2 to km^2
        mean_areas.append(m_area)

        vol = np.sum(mcs["size"] * mcs["mean"]) / 1e6  # Volume
        volumes.append(vol)

        hp_vol = np.sum(mcs["size"] * mcs["max"])  # High-precipitation volume
        hp_volumes.append(hp_vol)

        sev = np.mean(mcs["tot"])  # Severity
        severities.append(sev)

        # Intensity / Precipitation metrics (pr_HPE)
        pr_mean.append(np.mean(mcs["mean"]))
        pr_p10.append(np.percentile(mcs["mean"], 10))
        pr_p25.append(np.percentile(mcs["mean"], 25))
        pr_median.append(np.median(mcs["mean"]))
        pr_p75.append(np.percentile(mcs["mean"], 75))
        pr_p90.append(np.percentile(mcs["mean"], 90))
        pr_p99.append(np.percentile(mcs["mean"], 99))
        pr_max.append(np.max(mcs["max"]))

    # Dictionary for each metric
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

# Load datasets
mcs_eur_gpm = open_mcs("/GPM/EURR-3/output", start="2000-01", end="2009-12")
mcs_eur_cpm_eval = open_mcs("/CPMs/ICTP/EURR-3/evaluation/ERA5/output", start="2000-01", end="2009-12",)
mcs_eur_cpm_hist = open_mcs("/CPMs/ICTP/EURR-3/historical/ECEarth/output", start="2000-01", end="2009-12",)
mcs_eur_rcm_eval = open_mcs("/RCMs/ICTP/EUR-12/evaluation/ERA5/output", start="2000-01", end="2009-12",)
mcs_eur_rcm_hist = open_mcs("/RCMs/ICTP/EUR-12/historical/ECEarth/output", start="2000-01", end="2009-12",)

# Compute metrics 
gpm_medians, gpm_raw, n_gpm = compute_mcs_metrics(mcs_eur_gpm)
cpm_eval_medians, cpm_eval_raw, n_cpm_eval = compute_mcs_metrics(mcs_eur_cpm_eval)
cpm_hist_medians, cpm_hist_raw, n_cpm_hist = compute_mcs_metrics(mcs_eur_cpm_hist)
rcm_eval_medians, rcm_eval_raw, n_rcm_eval = compute_mcs_metrics(mcs_eur_rcm_eval)
rcm_hist_medians, rcm_hist_raw, n_rcm_hist = compute_mcs_metrics(mcs_eur_rcm_hist)

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

labels = list(gpm_medians.keys())
values_gpm = np.array(list(gpm_medians.values()))
values_cpm_eval = np.array(list(cpm_eval_medians.values()))
values_cpm_hist = np.array(list(cpm_hist_medians.values()))
values_rcm_eval = np.array(list(rcm_eval_medians.values()))
values_rcm_hist = np.array(list(rcm_hist_medians.values()))
num_vars = len(labels)

# Normalize each variable 
max_axis_limits = values_gpm * 1.5

norm_gpm = values_gpm / max_axis_limits
norm_gpm = np.concatenate((norm_gpm, [norm_gpm[0]]))

norm_cpm_eval = values_cpm_eval / max_axis_limits
norm_cpm_eval = np.concatenate((norm_cpm_eval, [norm_cpm_eval[0]]))

norm_cpm_hist = values_cpm_hist / max_axis_limits
norm_cpm_hist = np.concatenate((norm_cpm_hist, [norm_cpm_hist[0]]))

norm_rcm_eval = values_rcm_eval / max_axis_limits
norm_rcm_eval = np.concatenate((norm_rcm_eval, [norm_rcm_eval[0]]))

norm_rcm_hist = values_rcm_hist / max_axis_limits
norm_rcm_hist = np.concatenate((norm_rcm_hist, [norm_rcm_hist[0]]))

angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
angles += angles[:1]

# Create Plot
fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))

ax.plot(angles, norm_gpm, color="#4daf4a", linewidth=0.8, marker="o", markersize=3, label="GPM",)
ax.fill(angles, norm_gpm, color="#4daf4a", alpha=0.35)
ax.plot(angles, norm_cpm_eval, color="#1f77b4", linewidth=0.8, marker="s", markersize=3, label="CPM-3 Eval",)
ax.plot(angles, norm_cpm_hist, color="#41b6c4", linewidth=0.8, marker="^", markersize=3, label="CPM-3 Hist",)
ax.plot(angles, norm_rcm_eval, color="#225ea8", linewidth=0.8, marker="d", markersize=3, label="RCM-12 Eval",)
ax.plot(angles, norm_rcm_hist, color="#081d58", linewidth=0.8, marker="x", markersize=3, label="RCM-12 Hist",)

ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)
ax.set_xticks(angles[:-1])
ax.set_yticks(np.linspace(0.05, 1.0, 20))
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_ylim(0, 1)
ax.grid(True, color="black", linestyle="-", linewidth=0.3, alpha=0.4)

# Perimeter labels 
for i, angle in enumerate(angles[:-1]):
    lbl = labels[i]
    val = values_gpm[i]
    unit = units[lbl]
    val_str = f"{val/1000:.1f}" if val >= 1000 else f"{val:.1f}"

    text_label = f"${lbl}$\n$M_{{obs}} = {val_str}\\ {unit}$"
    ha = "left" if np.cos(angle) >= 0 else "right"
    ax.text(angle, 1.15, text_label, size=10, horizontalalignment=ha, verticalalignment="center",)

# Title and annotations
plt.title("Allyear\nAllover", loc="left", fontsize=12, fontweight="bold", pad=30)
ax.text(-0.20, 0.95, f"GPM: N={n_gpm}", transform=ax.transAxes, color="#4daf4a", fontsize=10, fontweight="bold",)
ax.text(-0.20, 0.91, f"CPM-3 Eval: N={n_cpm_eval}", transform=ax.transAxes, color="#1f77b4", fontsize=10, fontweight="bold",)
ax.text(-0.20, 0.87, f"CPM-3 Hist: N={n_cpm_hist}", transform=ax.transAxes, color="#41b6c4", fontsize=10, fontweight="bold", )
ax.text(-0.20, 0.83, f"RCM-12 Eval: N={n_rcm_eval}", transform=ax.transAxes, color="#225ea8", fontsize=10, fontweight="bold",)
ax.text(-0.20, 0.79, f"RCM-12 Hist: N={n_rcm_hist}", transform=ax.transAxes, color="#081d58", fontsize=10, fontweight="bold",)

# Save figure
path_out = ("/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/figs/v2")
name_out = f"pyplt_graphs_moaap_mcs_radar_charac_{domain}_2000-2009.png"
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches="tight", facecolor="white", edgecolor="none",)
plt.show()

