# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 02, 2026"
__description__ = "This script plot trend heatmap"

import glob
import os
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.colors import BoundaryNorm
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import linregress

# Indices
idx = "RX1day" # RX1day RR1 TXx TN20
idx_ = "RX1DAY" # RX1DAY RR1 TXx TN20

# Path setup
obs_dir = f"/home/netapp-clima-scratch/mda_silv/CORDEX-CORE2/OBS/CPC/{idx}_ts"

# CMIP5 GCMs
cmip5_dir = f"/home/mda_silv/scratch/CORDEX-CORE2/CMIP5/{idx}_ts"
cmip5_m1 = "ICHEC-EC-EARTH_r1"
cmip5_m2 = "MPI-M-MPI-ESM-LR_r1"
cmip5_m3 = "NCC-NorESM1-M_r1"

# CMIP6 GCMs
cmip6_m1_dir = f"/home/mda_silv/scratch/CORDEX-CORE2/CMIP6/EC-Earth3-Veg/{idx}_ts"
cmip6_m2_dir = f"/home/mda_silv/scratch/CORDEX-CORE2/CMIP6/MPI-ESM1-2-HR/{idx}_ts"
cmip6_m3_dir = f"/home/mda_silv/scratch/CORDEX-CORE2/CMIP6/NorESM2-MM/{idx}_ts"

# CORDEX-CORE 1 base directory
cordex_core1_base = "/home/clima-archive4-b/nzazulie/trendsCORDEX"

# CORDEX-CORE 2 base directory
cordex_core2_base = "/home/clima-archive4-b/nzazulie/trendsCORDEX/CORDEX-CORE2"

# Base model tag prefixes for CORDEX-CORE 1
rcm4_tag_m1 = "MOHC-HadGEM2-ES_r1_ICTP-RegCM4"
rcm4_tag_m2 = "MPI-M-MPI-ESM-MR_r1_ICTP-RegCM4"
rcm4_tag_m3 = "NCC-NorESM1-M_r1_ICTP-RegCM4"
remo_tag_m1 = "NCC-NorESM1-M_r1_GERICS-REMO2015"
remo_tag_m2 = "MPI-M-MPI-ESM-LR_r1_GERICS-REMO2015"

# CORDEX-CORE 2 GCM tags for RegCM5
rcm5_gcm1 = "EC-Earth3-Veg"
rcm5_gcm2 = "MPI-ESM1-2-HR"
rcm5_gcm3 = "NorESM2-MM"

# Mapping to CORDEX domains
domain_dict = {
    "NAM-22": ["ALA", "ANA", "NWN", "NEN", "CWN", "WSN", "CNA", "CEN", "SEN", "ESN", "WSN"],
    "CAM-22": ["CAM", "CAR"],
    "SAM-22": ["NWS", "NSA", "SAM", "NES", "SES", "SWS", "SSA"],
    "AFR-22": ["NWSH", "NCSH", "NESH", "SWSH", "SCSH", "SESH", "WAF", "NCAF", "NEAF", "SEAF", "ARP", "SCAF", "WSAF", "ESAF", "MDG"],
    "WAS-22": ["WCA", "SAS"],
    "EAS-22": ["ESB", "RFE", "ECA", "TIB", "EAS"],
    "SEA-22": ["NWSEA", "NESEA", "SSEA"],
    "AUS-22": ["NAU", "CAU", "EAU", "WSAU", "ESAU", "NZE"],
    "EUR-11": ["SMED", "NMED", "NEU", "WCE"],
}

# Domain mapping lookup table from 22/11 km to 12 km (CORDEX-CORE 2)
dom_22_to_12 = {
    "NAM-22": "NAM-12",
    "CAM-22": "CAM-12",
    "SAM-22": "SAM-12",
    "AFR-22": "AFR-12",
    "AUS-22": "AUS-12",
    "EUR-11": "EUR-12",
}

# Domain grouping order & subregion-to-domain lookup map
ordered_subregs = []
subreg_to_domain = {}

for dom, subregs in domain_dict.items():
  for reg in subregs:
    if reg not in ordered_subregs:
      ordered_subregs.append(reg)
      subreg_to_domain[reg] = dom


# Helper function to find dynamic RegCM file paths using wildcard globbing
def find_rcm_file(folder, search_pattern):
  pattern = os.path.join(folder, search_pattern)
  matches = glob.glob(pattern)
  if matches:
    return matches[0]
  return None


# Function to automatically resolve RegCM4 version (-7 or -6)
def resolve_rcm4_path(base_dir, idx_str, tag_prefix, subreg):
  path_7 = os.path.join(base_dir, f"{idx_str}_{tag_prefix}-7_{subreg}_ts.nc")
  path_6 = os.path.join(base_dir, f"{idx_str}_{tag_prefix}-6_{subreg}_ts.nc")

  if os.path.exists(path_7):
    return path_7
  elif os.path.exists(path_6):
    return path_6
  return path_7


# Function to compute spatial mean and linear trend (sliced 1970\u20132024)
def get_stats(file_path, start_year=1970, end_year=2024):
  if file_path is None or not os.path.exists(file_path):
    return np.nan, np.nan
  try:
    # Explicitly set decode_timedelta=False to prevent xarray from treating TN20 as a timedelta object
    ds = xr.open_dataset(
        file_path, decode_times=True, decode_timedelta=False, use_cftime=True
    )
    var_name = list(ds.data_vars)[0]
    da = ds[var_name].squeeze()

    # Select target time period (1970-2024)
    if "time" in da.coords:
      try:
        da = da.sel(time=slice(f"{start_year}-01-01", f"{end_year}-12-31"))
      except Exception:
        # Fallback for integer/float time coordinates or unparsed cftime years
        years = da.time.dt.year if hasattr(da.time, "dt") else da.time.values
        da = da.isel(time=(years >= start_year) & (years <= end_year))

    # Spatial/Regional Mean
    mean_val = float(da.mean().values)

    # Linear Trend
    y = da.values
    x = np.arange(len(y))
    valid = ~np.isnan(y)

    if np.sum(valid) < 2:
      return mean_val, np.nan

    slope, _, _, _, _ = linregress(x[valid], y[valid])
    return mean_val, slope

  except Exception as e:
    print(f"Error reading {file_path}: {e}")
    return np.nan, np.nan


# Data extraction containers
obs_means, obs_trends = [], []

cmip5_m1_means, cmip5_m1_trends = [], []
cmip5_m2_means, cmip5_m2_trends = [], []
cmip5_m3_means, cmip5_m3_trends = [], []

cmip6_m1_means, cmip6_m1_trends = [], []
cmip6_m2_means, cmip6_m2_trends = [], []
cmip6_m3_means, cmip6_m3_trends = [], []

rcm4_m1_means, rcm4_m1_trends = [], []
rcm4_m2_means, rcm4_m2_trends = [], []
rcm4_m3_means, rcm4_m3_trends = [], []

remo_m1_means, remo_m1_trends = [], []
remo_m2_means, remo_m2_trends = [], []

rcm5_m1_means, rcm5_m1_trends = [], []
rcm5_m2_means, rcm5_m2_trends = [], []
rcm5_m3_means, rcm5_m3_trends = [], []

for reg in ordered_subregs:
  dom = subreg_to_domain[reg]

  # OBS
  m_obs, t_obs = get_stats(os.path.join(obs_dir, f"{idx}_{reg}_ts.nc"))
  obs_means.append(m_obs)
  obs_trends.append(t_obs)

  # CMIP5
  m_c5_1, t_c5_1 = get_stats(os.path.join(cmip5_dir, f"{idx}_{cmip5_m1}_{reg}_ts.nc"))
  cmip5_m1_means.append(m_c5_1)
  cmip5_m1_trends.append(t_c5_1)
  
  m_c5_2, t_c5_2 = get_stats(os.path.join(cmip5_dir, f"{idx}_{cmip5_m2}_{reg}_ts.nc"))
  cmip5_m2_means.append(m_c5_2)
  cmip5_m2_trends.append(t_c5_2)

  m_c5_3, t_c5_3 = get_stats(os.path.join(cmip5_dir, f"{idx}_{cmip5_m3}_{reg}_ts.nc"))
  cmip5_m3_means.append(m_c5_3)
  cmip5_m3_trends.append(t_c5_3)

  # CMIP6
  m_c6_1, t_c6_1 = get_stats(os.path.join(cmip6_m1_dir, f"{idx}_{reg}_ts.nc"))
  cmip6_m1_means.append(m_c6_1)
  cmip6_m1_trends.append(t_c6_1)

  m_c6_2, t_c6_2 = get_stats(os.path.join(cmip6_m2_dir, f"{idx}_{reg}_ts.nc"))
  cmip6_m2_means.append(m_c6_2)
  cmip6_m2_trends.append(t_c6_2)

  m_c6_3, t_c6_3 = get_stats(os.path.join(cmip6_m3_dir, f"{idx}_{reg}_ts.nc"))
  cmip6_m3_means.append(m_c6_3)
  cmip6_m3_trends.append(t_c6_3)

  # CORDEX-CORE 1 (RegCM4 & REMO2015)
  core1_dir = os.path.join(cordex_core1_base, dom, f"{idx_}_ts_AR7")

  rcm4_m1_path = resolve_rcm4_path(core1_dir, idx_, rcm4_tag_m1, reg)
  m_rcm4_1, t_rcm4_1 = get_stats(rcm4_m1_path)
  rcm4_m1_means.append(m_rcm4_1)
  rcm4_m1_trends.append(t_rcm4_1)

  rcm4_m2_path = resolve_rcm4_path(core1_dir, idx_, rcm4_tag_m2, reg)
  m_rcm4_2, t_rcm4_2 = get_stats(rcm4_m2_path)
  rcm4_m2_means.append(m_rcm4_2)
  rcm4_m2_trends.append(t_rcm4_2)

  rcm4_m3_path = resolve_rcm4_path(core1_dir, idx_, rcm4_tag_m3, reg)
  m_rcm4_3, t_rcm4_3 = get_stats(rcm4_m3_path)
  rcm4_m3_means.append(m_rcm4_3)
  rcm4_m3_trends.append(t_rcm4_3)

  p_remo_1 = find_rcm_file(core1_dir, f"{idx_}_{remo_tag_m1}_{reg}_ts.nc")
  m_remo_1, t_remo_1 = get_stats(p_remo_1, start_year=1970, end_year=2024)
  remo_m1_means.append(m_remo_1)
  remo_m1_trends.append(t_remo_1)

  p_remo_2 = find_rcm_file(core1_dir, f"{idx_}_{remo_tag_m2}_{reg}_ts.nc")
  m_remo_2, t_remo_2 = get_stats(p_remo_2, start_year=1970, end_year=2024)
  remo_m2_means.append(m_remo_2)
  remo_m2_trends.append(t_remo_2)

  # CORDEX-CORE 2 (RegCM5 - 12km)
  dom12 = dom_22_to_12.get(dom, None)
  if dom12:
    core2_dir = os.path.join(cordex_core2_base, dom12, f"{idx_}_ts")

    # Locate exact RegCM5 files dynamically using wildcard matching
    p_rcm5_1 = find_rcm_file(core2_dir, f"{idx_}_{dom12}_{rcm5_gcm1}_*_{reg}_ts.nc")
    p_rcm5_2 = find_rcm_file(core2_dir, f"{idx_}_{dom12}_{rcm5_gcm2}_*_{reg}_ts.nc")
    p_rcm5_3 = find_rcm_file(core2_dir, f"{idx_}_{dom12}_{rcm5_gcm3}_*_{reg}_ts.nc")

    m_rcm5_1, t_rcm5_1 = get_stats(p_rcm5_1, start_year=1970, end_year=2024)
    m_rcm5_2, t_rcm5_2 = get_stats(p_rcm5_2, start_year=1970, end_year=2024)
    m_rcm5_3, t_rcm5_3 = get_stats(p_rcm5_3, start_year=1970, end_year=2024)
  else:
    m_rcm5_1, t_rcm5_1 = np.nan, np.nan
    m_rcm5_2, t_rcm5_2 = np.nan, np.nan
    m_rcm5_3, t_rcm5_3 = np.nan, np.nan

  rcm5_m1_means.append(m_rcm5_1)
  rcm5_m1_trends.append(t_rcm5_1)
  rcm5_m2_means.append(m_rcm5_2)
  rcm5_m2_trends.append(t_rcm5_2)
  rcm5_m3_means.append(m_rcm5_3)
  rcm5_m3_trends.append(t_rcm5_3)

# Combine matrices
means_matrix = np.column_stack([
    obs_means,
    cmip5_m1_means,
    cmip5_m2_means,
    cmip5_m3_means,
    cmip6_m1_means,
    cmip6_m2_means,
    cmip6_m3_means,
    rcm4_m1_means,
    rcm4_m2_means,
    rcm4_m3_means,
    remo_m1_means,
    remo_m2_means,
    rcm5_m1_means,
    rcm5_m2_means,
    rcm5_m3_means,
])

trends_matrix = np.column_stack([
    obs_trends,
    cmip5_m1_trends,
    cmip5_m2_trends,
    cmip5_m3_trends,
    cmip6_m1_trends,
    cmip6_m2_trends,
    cmip6_m3_trends,
    rcm4_m1_trends,
    rcm4_m2_trends,
    rcm4_m3_trends,
    remo_m1_trends,
    remo_m2_trends,
    rcm5_m1_trends,
    rcm5_m2_trends,
    rcm5_m3_trends,
])

n_domains = len(ordered_subregs)
n_cols = 15
cols = [
    "OBS",
    "EC-Earth",
    "MPI-ESM-MR",
    "NorESM1-M",
    "EC-Earth3-Veg",
    "MPI-ESM1-2-HR",
    "NorESM2-MM",
    "RegCM4-HadGEM2-ES",
    "RegCM4-MPI",
    "RegCM4-Nor",
    "REMO2015-Nor",
    "REMO2015-MPI",
    "RegCM5-EC-Earth3",
    "RegCM5-MPI",
    "RegCM5-Nor",
]

# Discrete levels for Mean and Trend
levels_mean = np.linspace(0, 100, 21)
levels_trend = np.linspace(-2, 2, 21)

if idx == "RX1day":
  cmap_mean = plt.colormaps["Blues"].resampled(len(levels_mean) - 1)
  cmap_trend = plt.colormaps["BrBG"].resampled(len(levels_trend) - 1)
else:
  cmap_mean = plt.colormaps["Reds"].resampled(len(levels_mean) - 1)
  cmap_trend = plt.colormaps["bwr"].resampled(len(levels_trend) - 1)

norm_mean = BoundaryNorm(levels_mean, ncolors=cmap_mean.N, clip=True)
norm_trend = BoundaryNorm(levels_trend, ncolors=cmap_trend.N, clip=True)

# Plot heatmap
fig, ax = plt.subplots(figsize=(15, 14))
ax.set_facecolor("white")
pad = 0.08

for row in range(n_domains):
  for col in range(n_cols):
    m_val = means_matrix[row, col]
    t_val = trends_matrix[row, col]

    x_min, x_max = col + pad, col + 1 - pad
    y_min, y_max = row + pad, row + 1 - pad

    # Mean triangle (upper-left split)
    if not np.isnan(m_val):
      p1 = Path([(x_min, y_min), (x_max, y_min), (x_min, y_max), (x_min, y_min)])
      patch1 = PathPatch(p1, facecolor=cmap_mean(norm_mean(m_val)), edgecolor="none")
      ax.add_patch(patch1)

    # Trend triangle (lower-right split)
    if not np.isnan(t_val):
      p2 = Path([(x_max, y_min), (x_max, y_max), (x_min, y_max), (x_max, y_min)])
      patch2 = PathPatch(p2, facecolor=cmap_trend(norm_trend(t_val)), edgecolor="none")
      ax.add_patch(patch2)

# Configure axes
ax.set_aspect("equal")
ax.set_xlim(0, n_cols)
ax.set_ylim(0, n_domains)
ax.set_xticks(np.arange(n_cols) + 0.5)
ax.set_xticklabels(cols, fontsize=10, rotation=30, ha="right")
ax.set_yticks(np.arange(n_domains) + 0.5)
ax.set_yticklabels(ordered_subregs, fontsize=10)
ax.invert_yaxis()

for spine in ax.spines.values():
  spine.set_visible(False)

# Colorbars setup
sm_mean = plt.cm.ScalarMappable(cmap=cmap_mean, norm=norm_mean)
sm_trend = plt.cm.ScalarMappable(cmap=cmap_trend, norm=norm_trend)

cax_mean = inset_axes(ax, width="4%", height="42%", loc="upper right", bbox_to_anchor=(0.22, 0.0, 1.0, 1.0), bbox_transform=ax.transAxes,)
cbar_mean = plt.colorbar(sm_mean, cax=cax_mean, orientation="vertical", ticks=levels_mean)
cbar_mean.set_label("Mean", fontsize=10)

cax_trend = inset_axes(ax, width="4%", height="42%", loc="lower right", bbox_to_anchor=(0.22, 0.0, 1.0, 1.0),bbox_transform=ax.transAxes,)
cbar_trend = plt.colorbar(sm_trend, cax=cax_trend, orientation="vertical", ticks=levels_trend)
cbar_trend.set_label("Trend", fontsize=10)

# Save figure
plt.tight_layout()
out_dir = "/home/mda_silv/scratch/CORDEX-CORE2/figs"
os.makedirs(out_dir, exist_ok=True)
plt.savefig(f"{out_dir}/{idx}_CORDEX-CORE2_heatmap_AR7.png", dpi=400, bbox_inches='tight')
plt.show()
