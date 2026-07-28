# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 28, 2026"
__description__ = "This script plot Taylor diagram"

import os
import gc
import glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.floating_axes as FA
import mpl_toolkits.axisartist.grid_finder as GF


def plot_taylor_diagram(stds, ccoefs, labels, fig=None, rect=111, title=""):
    """
    Plots a Taylor Diagram using native Matplotlib floating axes.
    stds: list of normalized standard deviations [mod_1, mod_2, ...]
    ccoefs: list of correlation coefficients [mod_1, mod_2, ...]
    labels: list of labels ['CPC', 'DJF', 'MAM', 'JJA', 'SON']
    """
    if fig is None:
        fig = plt.gcf()

    ref_std = 1.0  # Reference normalized std for observations
    tr = PolarAxes.PolarTransform()

    # Correlation locator and format
    rlocs = np.array([0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0])
    tlocs = np.arccos(rlocs)
    gl1 = GF.FixedLocator(tlocs)
    tf1 = GF.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

    # Define standard deviation scale range
    max_std = max(ref_std, max(stds)) if len(stds) > 0 else ref_std
    smax = 1.5 * max_std

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

    # Plot reference observation point
    polar_ax.plot([0], [ref_std], 'k*', markersize=12, label=labels[0])

    # Reference STD arc
    t = np.linspace(0, np.pi / 2, 100)
    polar_ax.plot(t, np.full_like(t, ref_std), 'k--', linewidth=1.0)

    # RMSD contours
    rs, ts = np.meshgrid(
        np.linspace(0, smax, 100),
        np.linspace(0, np.pi / 2, 100)
    )
    rms = np.sqrt(ref_std**2 + rs**2 - 2 * ref_std * rs * np.cos(ts))
    contours = polar_ax.contour(ts, rs, rms, colors='0.5', linestyles=':', linewidths=0.8)
    polar_ax.clabel(contours, inline=1, fontsize=7, fmt='%.1f')

    # Plot model markers per season
    colors = ['tab:blue', 'tab:blue', 'tab:blue', 'tab:blue']
    markers = ['o', 's', '^', 'D']

    for i, (std, cc, label) in enumerate(zip(stds, ccoefs, labels[1:])):
        theta = np.arccos(np.clip(cc, -1.0, 1.0))
        polar_ax.plot(
            theta, std,
            marker=markers[i % len(markers)],
            color=colors[i % len(colors)],
            markersize=8,
            ls='',
            label=label
        )

    if title:
        ax.set_title(title, pad=25, fontweight='bold')

    polar_ax.legend(loc='upper right', bbox_to_anchor=(1.30, 1.15), frameon=True)
    return ax


# Input path
path_data = "/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/postproc/urban/paper"

# Variables dictionary
variables = {
    'pr': {
        'obs': 'precip',
        'mod': 'pr'
    },
    'tasmax': {
        'obs': 'tmax',
        'mod': 'tasmax'
    },
    'tasmin': {
        'obs': 'tmin',
        'mod': 'tasmin'
    }
}

seasons = ['DJF', 'MAM', 'JJA', 'SON']

# Create figure 
fig = plt.figure(figsize=(18, 6))

for ivar, var in enumerate(variables):

    std_norm = []
    ccoef = []

    for season in seasons:

        # File names
        file_mod = os.path.join(path_data, f'{var}_CSAM-3_RegCM5-ERA5_URBAN_{season}_2000-2009_0.25_box.nc')
        file_obs = os.path.join(path_data, f"{variables[var]['obs']}_CSAM-3_CPC_{season}_2000-2009_0.25_box.nc")
        print(f"Mod: {file_mod}")
        print(f"Obs: {file_obs}")

        # Read data 
        with xr.open_dataset(file_obs) as ds_obs:
            obs = ds_obs[variables[var]['obs']].squeeze().values.astype(np.float32)

        with xr.open_dataset(file_mod) as ds_mod:
            mod = ds_mod[variables[var]['mod']].squeeze().values.astype(np.float32)

        # Spatial pattern 
        obs = obs.ravel()
        mod = mod.ravel()

        mask = np.isfinite(obs) & np.isfinite(mod)

        obs = obs[mask]
        mod = mod[mask]

        # Calculate statistics
        std_obs = np.std(obs, ddof=1)
        std_mod = np.std(mod, ddof=1)

        # Normalized standard deviation relative to observations
        std_norm.append(std_mod / std_obs)
        ccoef.append(np.corrcoef(obs, mod)[0, 1])

        del obs, mod, mask
        gc.collect()

    # Subplot placement 
    rect = 131 + ivar

    # Plot Taylor Diagram for the current variable
    plot_taylor_diagram(
        stds=std_norm,
        ccoefs=ccoef,
        labels=['CPC'] + seasons,
        fig=fig,
        rect=rect,
        title=var.upper()
    )

# Save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/figs/urb/paper'
os.makedirs(path_out, exist_ok=True)
name_out = 'pyplt_taylor_diagram_RegCM5_CSAM-3_2000-2009.png'

plt.tight_layout()
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
