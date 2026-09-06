# -*- coding: utf-8 -*-

# author      = "Leidinice Silva"
# email       = "leidinicesilva@gmail.com"
# date        = "Jul 28, 2026"
# description = "This script plots Taylor diagrams"

import os
import gc
import glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.floating_axes as FA
import mpl_toolkits.axisartist.grid_finder as GF


def setup_taylor_diagram(fig, subplot_spec, title="", max_std_val=1.5, ref_label='CPC'):
    """Sets up the Taylor Diagram coordinate system, grid, reference point, and RMSD contours."""
    ref_std = 1.0  # Reference normalized std for observations
    
    tr = PolarAxes.PolarTransform(apply_theta_transforms=False)

    # Correlation locator and format
    rlocs = np.array([0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0])
    tlocs = np.arccos(rlocs)
    gl1 = GF.FixedLocator(tlocs)
    tf1 = GF.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

    # Define standard deviation scale range dynamically
    smax = 1.25 * max(ref_std, max_std_val)

    ghelper = FA.GridHelperCurveLinear(
        tr,
        extremes=(0, np.pi / 2, 0, smax),
        grid_locator1=gl1,
        tick_formatter1=tf1
    )

    ax = FA.FloatingSubplot(fig, subplot_spec, grid_helper=ghelper)
    fig.add_subplot(ax)

    # Configure axes
    ax.axis["top"].set_axis_direction("bottom")
    ax.axis["top"].toggle(ticklabels=True, label=True)
    ax.axis["top"].major_ticklabels.set_axis_direction("top")
    ax.axis["top"].label.set_axis_direction("top")
    ax.axis["top"].label.set_text("Correlation")

    ax.axis["left"].set_axis_direction("bottom")
    ax.axis["left"].label.set_text("Norm. STD")

    ax.axis["right"].set_axis_direction("top")
    ax.axis["right"].toggle(ticklabels=True)
    ax.axis["right"].major_ticklabels.set_axis_direction("left")

    ax.axis["bottom"].set_visible(False)

    polar_ax = ax.get_aux_axes(tr)

    # Enable grid lines
    ax.grid(True, linestyle=':', color='gray', alpha=0.5)

    # Plot reference observation point
    polar_ax.plot([0], [ref_std], 'k*', markersize=10, label='OBS')

    # Reference STD arc
    t = np.linspace(0, np.pi / 2, 100)
    polar_ax.plot(t, np.full_like(t, ref_std), 'k--', linewidth=0.8)

    # RMSD contours (Centered RMSE)
    rs, ts = np.meshgrid(np.linspace(0, smax, 200), np.linspace(0, np.pi / 2, 200))
    rms = np.sqrt(ref_std**2 + rs**2 - 2 * ref_std * rs * np.cos(ts))

    contours = polar_ax.contour(ts, rs, rms, levels=4, colors='gray', linestyles='--', linewidths=0.6)
    polar_ax.clabel(contours, inline=1, fontsize=6, fmt='%.1f')

    if title:
        ax.set_title(title, loc='left', pad=18, fontweight='bold', fontsize=10)

    return ax, polar_ax


def add_experiment_points(polar_ax, std, ccoef, label, exp_color='tab:blue'):
    """Plots model marker for a specific experiment on an existing Taylor diagram axis."""
    theta = np.arccos(np.clip(ccoef, -1.0, 1.0))
    polar_ax.plot(theta, std, marker='o', color=exp_color, markersize=8, alpha=0.80, ls='', label=label)


# Input & Output paths
path_data = "/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11/postproc/paper/txt_files"
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11/figs/paper'
os.makedirs(path_out, exist_ok=True)

# Grid layout setup: 3 Variables (rows) x 5 Seasons (columns)
variables = ["pr", "tas", "clt"]
var_titles = {"pr": "PR", "tas": "TAS", "clt": "CLT"}
seasons = ['DJF', 'MAM', 'JJA', 'SON', 'ANN']

experiments = [
    {"name": "NoTo", "prefix_mid": "NoTo-EUR", "color": "red"},
    {"name": "WSM5", "prefix_mid": "WSM5-EUR", "color": "blue"},
    {"name": "WSM7", "prefix_mid": "WSM7-EUR", "color": "green"},
    {"name": "WDM7", "prefix_mid": "WDM7-EUR", "color": "orange"}
]

# 1. Load data for all variables
data_store = {}
all_stds = []

for var in variables:
    obs_name = "CPC" if var == "pr" else "ERA5"
    data_store[var] = {}

    for exp in experiments:
        prefix = f"{var}_RegCM5_{exp['prefix_mid']}_{obs_name}_2000-2009"
        ccoef = np.atleast_1d(np.loadtxt(os.path.join(path_data, f"{prefix}_cc.txt"))).flatten()
        std_norm = np.atleast_1d(np.loadtxt(os.path.join(path_data, f"{prefix}_ratio.txt"))).flatten()

        data_store[var][exp["name"]] = {
            "ccoef": ccoef,
            "std_norm": std_norm,
            "color": exp["color"]
        }
        all_stds.extend(std_norm)

max_std_val = max(all_stds)

# 2. Setup Figure and GridSpec (3 rows x 5 columns)
fig = plt.figure(figsize=(16, 10))
fig.patch.set_facecolor('#E0E0E0')
fig.patch.set_alpha(0.75)

gs = GridSpec(3, 5, figure=fig, hspace=0.35, wspace=0.25)

handles_list, labels_list = [], []
panel_idx = 0

# Outer loop: Rows (Variables) | Inner loop: Columns (Seasons)
for i_row, var in enumerate(variables):
    for j_col, (isea, sea) in enumerate(enumerate(seasons)):
        
        letter = chr(97 + panel_idx)  # (a), (b), (c)...
        subplot_title = f"({letter}) {var_titles[var]} - {sea}"

        ax, polar_ax = setup_taylor_diagram(
            fig=fig,
            subplot_spec=gs[i_row, j_col],
            title=subplot_title,
            max_std_val=max_std_val,
            ref_label=("CPC" if var == "pr" else "ERA5")
        )

        # Plot experiments for current variable and season
        for exp in experiments:
            add_experiment_points(
                polar_ax=polar_ax,
                std=data_store[var][exp["name"]]["std_norm"][isea],
                ccoef=data_store[var][exp["name"]]["ccoef"][isea],
                label=exp["name"],
                exp_color=data_store[var][exp["name"]]["color"]
            )

        # Grab handles from the first panel for global legend
        if i_row == 0 and j_col == 0:
            handles_list, labels_list = polar_ax.get_legend_handles_labels()

        panel_idx += 1

# 3. Global legend at the bottom
fig.legend(handles_list,labels_list,loc='lower center',bbox_to_anchor=(0.5, 0.01),ncol=len(labels_list),frameon=True,fontsize=10)

plt.subplots_adjust(left=0.04, right=0.96, top=0.94, bottom=0.08)

# Save figure
name_out = 'pyplt_taylor_diagram_RegCM5_EUR-11_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.close()
