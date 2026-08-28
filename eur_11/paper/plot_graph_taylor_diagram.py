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
    polar_ax.plot([0], [ref_std], 'k*', markersize=14, label='OBS')

    # Reference STD arc
    t = np.linspace(0, np.pi / 2, 100)
    polar_ax.plot(t, np.full_like(t, ref_std), 'k--', linewidth=1.0)

    # RMSD contours (Centered RMSE)
    rs, ts = np.meshgrid(np.linspace(0, smax, 200), np.linspace(0, np.pi / 2, 200))
    rms = np.sqrt(ref_std**2 + rs**2 - 2 * ref_std * rs * np.cos(ts))

    contours = polar_ax.contour(ts,rs,rms,levels=5,colors='gray',linestyles='--',linewidths=0.8)
    polar_ax.clabel(contours, inline=1, fontsize=8, fmt='%.1f')

    if title:
        ax.set_title(title, loc='left', pad=25, fontweight='bold')

    return ax, polar_ax


def add_experiment_points(polar_ax, std, ccoef, label, exp_color='tab:blue'):
    """Plots model marker for a specific experiment on an existing Taylor diagram axis."""
    theta = np.arccos(np.clip(ccoef, -1.0, 1.0))
    polar_ax.plot(theta,std,marker='o',color=exp_color,markersize=14,alpha=0.70,ls='',label=label)


# Input path
path_data = "/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11/postproc/paper/txt_files"
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11/figs/paper'
os.makedirs(path_out, exist_ok=True)

# Variables, seasons, and experiments setup
variables = ["pr", "tas", "clt"]
seasons = ['DJF', 'MAM', 'JJA', 'SON', 'ANN']

experiments = [
    {"name": "NoTo", "prefix_mid": "NoTo-EUR", "color": "red"},
    {"name": "WSM5", "prefix_mid": "WSM5-EUR", "color": "blue"},
    {"name": "WSM7", "prefix_mid": "WSM7-EUR", "color": "green"},
    {"name": "WDM7", "prefix_mid": "WDM7-EUR", "color": "orange"}
]

for var in variables:
    print(f"Taylor statistics: {var}")

    # Define observational dataset
    if var in ["pr"]:
        obs_name = "CPC"
    elif var in ["tas", "clt"]:
        obs_name = "ERA5"

    # Read data for each experiment
    exp_data = {}
    all_stds = []

    for exp in experiments:
        prefix = f"{var}_RegCM5_{exp['prefix_mid']}_{obs_name}_2000-2009"
        ccoef = np.atleast_1d(np.loadtxt(os.path.join(path_data, f"{prefix}_cc.txt"))).flatten()
        std_norm = np.atleast_1d(np.loadtxt(os.path.join(path_data, f"{prefix}_ratio.txt"))).flatten()
        
        exp_data[exp["name"]] = {
            "ccoef": ccoef,
            "std_norm": std_norm,
            "color": exp["color"]
        }
        all_stds.extend(std_norm)

    # Calculate overall max std to ensure proper axis scaling
    max_std_val = max(all_stds)

    # Create figure with 2 rows and 3 columns
    fig = plt.figure(figsize=(15, 10))
    handles_list, labels_list = [], []

    for isea, sea in enumerate(seasons):
        rect = 231 + isea
        subplot_title = f"({chr(97 + isea)}) {sea}"

        # Setup Taylor diagram for this season
        ax, polar_ax = setup_taylor_diagram(fig=fig,rect=rect,title=subplot_title,max_std_val=max_std_val,ref_label=obs_name)

        # Plot each experiment's point for the current season
        for exp in experiments:
            add_experiment_points(polar_ax=polar_ax,
                std=exp_data[exp["name"]]["std_norm"][isea],
                ccoef=exp_data[exp["name"]]["ccoef"][isea],
                label=exp["name"],
                exp_color=exp_data[exp["name"]]["color"])

        # Collect legend handles from the first panel
        if isea == 0:
            handles_list, labels_list = polar_ax.get_legend_handles_labels()

    # Create empty 6th subplot axis specifically for the legend
    ax_legend = fig.add_subplot(2, 3, 6)
    ax_legend.axis('off')  # Hide background grid and ticks

    # Place legend centered in the 6th panel
    ax_legend.legend(handles_list,labels_list,loc='center',frameon=True,fontsize=14,title="Experiments",title_fontsize=14)

    # Save figure for the variable
    name_out = f'pyplt_taylor_diagram_{var}_RegCM5_EUR-11_2000-2009.png'
    plt.tight_layout()
    plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
    plt.close()

