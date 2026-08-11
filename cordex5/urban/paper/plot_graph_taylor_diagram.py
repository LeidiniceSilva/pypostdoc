# -*- coding: utf-8 -*-

# author      = "Leidinice Silva"
# email       = "leidinicesilva@gmail.com"
# date        = "Jul 28, 2026"
# description = "This script plots Taylor diagrams with multiple experiments"

import os
import gc
import glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.floating_axes as FA
import mpl_toolkits.axisartist.grid_finder as GF


def setup_taylor_diagram(fig, rect, title="", max_std_val=1.5, ref_label='CPC'):
    """Sets up the Taylor Diagram coordinate system, grid, reference point, and RMSD contours."""
    ref_std = 1.0  # Reference normalized std for observations
    tr = PolarAxes.PolarTransform()

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

    ax = FA.FloatingSubplot(fig, rect, grid_helper=ghelper)
    fig.add_subplot(ax)

    # Configure axes
    ax.axis["top"].set_axis_direction("bottom")
    ax.axis["top"].toggle(ticklabels=True, label=True)
    ax.axis["top"].major_ticklabels.set_axis_direction("top")
    ax.axis["top"].label.set_axis_direction("top")
    ax.axis["top"].label.set_text("Correlation")

    ax.axis["left"].set_axis_direction("bottom")
    ax.axis["left"].label.set_text("Normalized Standard Deviation")

    ax.axis["right"].set_axis_direction("top")
    ax.axis["right"].toggle(ticklabels=True)
    ax.axis["right"].major_ticklabels.set_axis_direction("left")

    ax.axis["bottom"].set_visible(False)

    polar_ax = ax.get_aux_axes(tr)

    # Enable grid lines
    ax.grid(True, linestyle=':', color='gray', alpha=0.5)

    # Plot reference observation point
    polar_ax.plot([0], [ref_std], 'k*', markersize=12, label=ref_label)

    # Reference STD arc
    t = np.linspace(0, np.pi / 2, 100)
    polar_ax.plot(t, np.full_like(t, ref_std), 'k--', linewidth=1.0)

    # RMSD contours (Centered RMSE)
    rs, ts = np.meshgrid(
        np.linspace(0, smax, 200),
        np.linspace(0, np.pi / 2, 200)
    )

    rms = np.sqrt(ref_std**2 + rs**2 - 2 * ref_std * rs * np.cos(ts))

    contours = polar_ax.contour(
        ts,
        rs,
        rms,
        levels=5,
        colors='gray',
        linestyles='--',
        linewidths=0.8
    )

    polar_ax.clabel(contours, inline=1, fontsize=8, fmt='%.1f')

    if title:
        ax.set_title(title, loc='left', pad=25, fontweight='bold')

    return ax, polar_ax


def add_experiment_points(polar_ax, stds, ccoefs, labels, exp_color='tab:blue'):
    """Plots model markers for a specific experiment on an existing Taylor diagram axis."""
    markers = ['+', 'o', 's', '^', 'D']

    for i, (std, cc, label) in enumerate(zip(stds, ccoefs, labels)):
        theta = np.arccos(np.clip(cc, -1.0, 1.0))

        polar_ax.plot(
            theta,
            std,
            marker=markers[i % len(markers)],
            color=exp_color,
            markersize=8,
            ls='',
            label=label
        )


# Input path
path_data = "/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/postproc/urban/paper/txt_files"

# Variables and seasons
variables = ["pr", "tasmax", "tasmin"]
seasons = ['ANN', 'DJF', 'MAM', 'JJA', 'SON']

# Create figure
fig = plt.figure(figsize=(18, 6))

handles_list = []
labels_list = []

for ivar, var in enumerate(variables):
    print(f"Taylor statistics: {var}")

    # Read URBAN txt files
    prefix_urban = f"{var}_RegCM5-ERA5_URBAN_CPC_CSAM-3_2000-2009"
    ccoef_urban = np.atleast_1d(np.loadtxt(os.path.join(path_data, f"{prefix_urban}_cc.txt"))).flatten()
    std_norm_urban = np.atleast_1d(np.loadtxt(os.path.join(path_data, f"{prefix_urban}_ratio.txt"))).flatten()

    # Read CTRL txt files
    prefix_ctrl = f"{var}_RegCM5-ERA5_CTRL_CPC_CSAM-3_2000-2009"
    ccoef_ctrl = np.atleast_1d(np.loadtxt(os.path.join(path_data, f"{prefix_ctrl}_cc.txt"))).flatten()
    std_norm_ctrl = np.atleast_1d(np.loadtxt(os.path.join(path_data, f"{prefix_ctrl}_ratio.txt"))).flatten()

    rect = 131 + ivar
    subplot_title = f"({chr(97 + ivar)}) {var.upper()}"

    # Calculate overall max std to ensure outer limits fit both URBAN and CTRL
    max_std_val = max(np.max(std_norm_urban), np.max(std_norm_ctrl))

    # 1. Setup Taylor base plot ONCE for this variable
    ax, polar_ax = setup_taylor_diagram(
        fig=fig,
        rect=rect,
        title=subplot_title,
        max_std_val=max_std_val,
        ref_label='CPC'
    )

    # 2. Add URBAN experiment points
    add_experiment_points(
        polar_ax=polar_ax,
        stds=std_norm_urban,
        ccoefs=ccoef_urban,
        labels=[f"URBAN_{s}" for s in seasons],
        exp_color='tab:red'
    )

    # 3. Add CTRL experiment points onto the SAME polar_ax
    add_experiment_points(
        polar_ax=polar_ax,
        stds=std_norm_ctrl,
        ccoefs=ccoef_ctrl,
        labels=[f"CTRL_{s}" for s in seasons],
        exp_color='blue'
    )

    # Collect legend entries from the first subplot only
    if ivar == 0:
        handles_list, labels_list = polar_ax.get_legend_handles_labels()

# Single shared legend box below subplots
fig.legend(
    handles_list,
    labels_list,
    loc='lower center',
    bbox_to_anchor=(0.5, -0.1),
    ncol=len(labels_list),
    frameon=True,
    fontsize=10
)

# Save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/figs/urban/paper'
os.makedirs(path_out, exist_ok=True)
name_out = 'pyplt_taylor_diagram_RegCM5_CSAM-3_2000-2009.png'

plt.tight_layout()
plt.savefig(
    os.path.join(path_out, name_out),
    dpi=400,
    bbox_inches='tight'
)

plt.show()
