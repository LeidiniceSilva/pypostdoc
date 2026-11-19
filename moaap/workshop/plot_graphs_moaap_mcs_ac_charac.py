# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 03, 2026"
__description__ = "This script plot MCSs (fancy grouped bars)"

import os
import pickle
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

warnings.filterwarnings("ignore")


def open_mcs(path_, start="2000-01", end="2009-12"):
    path = f"/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/dataset/{path_}"
    dates = pd.date_range(start=start, end=end, freq="MS")
    mcs = {}

    for d in dates:
        fname = "MCSs_" + d.strftime("%Y%m") + "__dt-1h_MOAAP-masks.pkl"
        f = os.path.join(path, fname)

        if os.path.exists(f):
            with open(f, "rb") as file:
                mcs[d.strftime("%Y-%m")] = pickle.load(file)

    return mcs


def comp_annual_cycle_number(mcs_charac):
    count_by_month = np.zeros(12)
    years_by_month = [set() for i in range(12)]

    for obj in mcs_charac.keys():
        for m in mcs_charac[obj].keys():
            times = pd.DatetimeIndex(mcs_charac[obj][m]['times'])

            if len(times) == 0:
                continue

            month = times[0].month - 1
            year = times[0].year

            count_by_month[month] += 1
            years_by_month[month].add(year)

    num_clim = np.full(12, np.nan)
    for k in range(12):
        if len(years_by_month[k]) > 0:
            num_clim[k] = count_by_month[k] / len(years_by_month[k])

    return num_clim


def diagnostics(label, mcs):
    """Print numbers that reveal grid / time-step / duplication problems."""
    n_sys, n_steps = 0, 0
    dur, dt, size_km2, cells, cell_area = [], [], [], [], []
    keys = set()
    dup = 0
    years = set()
    months = set()

    for obj in mcs:
        for m in mcs[obj]:
            d = mcs[obj][m]
            t = pd.DatetimeIndex(d['times'])

            if len(t) == 0:
                continue

            n_sys += 1
            n_steps += len(t)

            years.update(t.year)
            months.add(obj)
            dur.append(len(t))

            if len(t) > 1:
                dt.append(np.median(np.diff(t).astype('timedelta64[m]').astype(float)))

            s = np.asarray(d['size'], dtype=float)
            if np.nanmedian(s) > 1e4:
                s = s / 1e6
            size_km2.append(np.nanmedian(s))

            mean = np.asarray(d['mean'], dtype=float)
            tot = np.asarray(d['tot'], dtype=float)

            ok = (mean > 0) & np.isfinite(mean) & np.isfinite(tot)

            if ok.any():
                c = tot[ok] / mean[ok]
                valid = c > 0

                if valid.any():
                    cells.append(np.nanmedian(c[valid]))
                    cell_area.append(np.nanmedian(s[ok][valid] / c[valid]))

            key = (t[0], t[-1], round(float(np.asarray(d['size'])[0])))
            dup += key in keys
            keys.add(key)

    print(f"{label:18s} systems={n_sys:7d} hourly-steps={n_steps:8d} "
          f"years={len(years):2d} months={len(months):3d} "
          f"median_dur={np.nanmedian(dur) if dur else np.nan:5.1f} "
          f"step={np.nanmedian(dt) if dt else np.nan:5.0f}min "
          f"median_size={np.nanmedian(size_km2) if size_km2 else np.nan:9.0f}km2 "
          f"pixels/MCS={np.nanmedian(cells) if cells else np.nan:8.0f} "
          f"cell_area={np.nanmedian(cell_area) if cell_area else np.nan:7.1f}km2 "
          f"duplicates={dup}")


def comp_annual_cycle_max_precip(mcs_charac):
    monthly_max = {k: [] for k in range(1, 13)}

    for obj in mcs_charac.keys():
        for m in mcs_charac[obj].keys():
            times = pd.DatetimeIndex(mcs_charac[obj][m]['times'])
            if len(times) == 0:
                continue

            month_idx = times[0].month
            max_pr = np.asarray(mcs_charac[obj][m]['max'], dtype=float)

            # PASSO 1: Média ao longo do ciclo de vida DESTE objeto específico (ou use np.nanmax se preferir o pico absoluto)
            mean_max_obj = np.nanmean(max_pr)

            # PASSO 2: Agrupa na lista correspondente ao mês do ano
            if np.isfinite(mean_max_obj):
                monthly_max[month_idx].append(mean_max_obj)

    max_clim = np.full(12, np.nan)

    # PASSO 3: Média de todos os objetos de cada mês ao longo dos anos
    for k in range(1, 13):
        if monthly_max[k]:
            max_clim[k - 1] = np.nanmean(monthly_max[k])

    return max_clim


def select_mcs_in_box(mcs, lon_min, lon_max, lat_min, lat_max, mode="any"):
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


# Domain
domain = "EUR"

# Selection box (lon_min, lon_max, lat_min, lat_max)
box = dict(lon_min=-10, lon_max=20, lat_min=38, lat_max=53)

# Datasets: label -> (path, color, hatch)
datasets = {
    'ERA5':             ("/ERA5/EURR-3/output",                              '#222222', ''),
    'GPM':              ("/GPM/EURR-3/output",                               '#6c6c6c', ''),
    'CNRM CPM-3 Eval':  ("/CPMs/CNRM/ALPX-3/evaluation/ERA5/output",         '#d1495b', ''),
    'CNRM CPM-3 Hist':  ("/CPMs/CNRM/ALPX-3/historical/CNRM-ESM2-1/output",  '#d1495b', '///'),
    'CNRM RCM-12 Eval': ("/RCMs/CNRM/EUR-12/evaluation/ERA5/output",         '#2e86ab', ''),
    'CNRM RCM-12 Hist': ("/RCMs/CNRM/EUR-12/historical/CNRM-ESM2-1/output",  '#2e86ab', '///'),
    'ICTP RCM-12 Eval': ("/RCMs/ICTP/EUR-12/evaluation/ERA5/output",         '#edae49', ''),
    'ICTP RCM-12 Hist': ("/RCMs/ICTP/EUR-12/historical/ECEarth/output",      '#edae49', '///'),
}

# Load datasets and compute metrics
results = {}

for label, (p, color, hatch) in datasets.items():
    print(f"Processando: {label}...")
    mcs = open_mcs(p, start="2000-01", end="2009-12")
    mcs = select_mcs_in_box(mcs, **box, mode="any")

    diagnostics(label, mcs)

    num = comp_annual_cycle_number(mcs)
    max_pr = comp_annual_cycle_max_precip(mcs)

    results[label] = dict(num=num, max_pr=max_pr)


# Plot parameters
plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.edgecolor': '#444444',
    'axes.linewidth': 1.0,
    'xtick.color': '#333333',
    'ytick.color': '#333333',
})

# Layout 1 linha e 2 colunas
fig = plt.figure(figsize=(16, 6), facecolor='white')

gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.18,
                       left=0.06, right=0.98, top=0.85, bottom=0.18)

font_size = 11
n = len(datasets)
width = 0.8 / n

time = np.arange(12)

xtick = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
         'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

panels = [
    ('num',    '(a) MCS count',                 'MCS per month'),
    ('max_pr', '(b) Max precipitation intensity', 'Max precipitation (mm h⁻¹)'),
]

for i, (key, title, ylabel) in enumerate(panels):
    ax = fig.add_subplot(gs[0, i])
    ax.set_facecolor('#fafafa')

    # Highlight summer (JJA) background
    ax.axvspan(4.5, 7.5, color='#fff1d6', alpha=0.6, zorder=0, lw=0)

    for j, (label, (p, color, hatch)) in enumerate(datasets.items()):
        offset = (j - (n - 1) / 2) * width

        ax.bar(time + offset, results[label][key], width * 0.95,
               color=color, hatch=hatch, alpha=0.9,
               edgecolor='white', linewidth=0.6, zorder=3)

    ax.set_title(title, fontsize=font_size + 3, loc='left',
                 fontweight='bold', color='#222222', pad=10)

    ax.set_ylabel(ylabel, fontsize=font_size)
    ax.set_xticks(time)
    ax.set_xticklabels(xtick, fontsize=font_size - 1)
    ax.set_xlim(-0.6, 11.6)
    ax.tick_params(axis='both', length=3)
    ax.grid(True, axis='y', linestyle='--', alpha=0.4, zorder=1)
    ax.set_xlabel('Month', fontsize=font_size)


# Shared legend: colours = model, hatch = historical
handles = [Patch(facecolor=c, edgecolor='white', hatch=h, label=l)
           for l, (_, c, h) in datasets.items()]

fig.legend(handles=handles, loc='lower center', ncol=4, fontsize=font_size,
           frameon=True, fancybox=True, framealpha=0.95, edgecolor='#cccccc',
           bbox_to_anchor=(0.52, 0.01))


# Path out to save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/paper/figs/workshop'
name_out = f'pyplt_graphs_moaap_mcs_ac_charac_{domain}_2000-2009_maxpr.png'

os.makedirs(path_out, exist_ok=True)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
