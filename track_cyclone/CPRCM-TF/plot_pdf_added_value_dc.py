# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot pdf"

import os
import netCDF4
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt

dt = '2000-2009'
var = 'pr'
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/TRACK-CYCLONE/CPRCM-TF/data/Dataset'

# Box: -60W to -40W, -35S to -15S
lon_min, lon_max = -60.0, -40.0
lat_min, lat_max = -30.0, -20.0

# Hourly pdf bins (mm/h) and daily pdf bins (mm/d)
step_h, rain_min_h, rain_max_h = 0.1, 0.1, 200.0
step_d, rain_min_d, rain_max_d = 1.0, 1.0, 600.0
bins_h_edges = np.arange(rain_min_h, rain_max_h + step_h, step_h)
bins_d_edges = np.arange(rain_min_d, rain_max_d + step_d, step_d)


def get_latlon(data):
    """Try multiple common lat/lon variable name conventions."""
    lat_names = ['lat', 'latitude', 'rlat', 'xlat']
    lon_names = ['lon', 'longitude', 'rlon', 'xlon']

    lat, lon = None, None
    for name in lat_names:
        if name in data.variables:
            lat = data.variables[name][:]
            break
    for name in lon_names:
        if name in data.variables:
            lon = data.variables[name][:]
            break

    if lat is None or lon is None:
        raise KeyError(f"Could not find lat/lon variables. Available: {list(data.variables.keys())}")

    return lat, lon


def get_box(lat, lon):
    """Return row/col slices and a 2D mask of the lat/lon box (1D and 2D lat/lon)."""
    lat = np.ma.filled(lat, np.nan)
    lon = np.ma.filled(lon, np.nan)
    lon = np.where(lon > 180, lon - 360, lon)  # 0-360 -> -180-180

    if lat.ndim == 1 and lon.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon, lat)
    else:
        lat2d, lon2d = lat, lon

    box = (lat2d >= lat_min) & (lat2d <= lat_max) & (lon2d >= lon_min) & (lon2d <= lon_max)
    if not box.any():
        raise ValueError('No grid points found inside the selected box.')

    rows = np.where(box.any(axis=1))[0]
    cols = np.where(box.any(axis=0))[0]
    rs = slice(rows.min(), rows.max() + 1)
    cs = slice(cols.min(), cols.max() + 1)

    return rs, cs, box[rs, cs]


def compute_pdfs(param, domain, dataset):
    """Read one day at a time inside the box and build hourly and daily histograms."""

    arq  = '{0}/{1}_{2}_{3}_mmhr_{4}_mask.nc'.format(path, param, domain, dataset, dt)
    data = netCDF4.Dataset(arq)
    vals = data.variables[param]
    lat, lon = get_latlon(data)
    rs, cs, box = get_box(lat, lon)

    hist_h = np.zeros(len(bins_h_edges) - 1)
    hist_d = np.zeros(len(bins_d_edges) - 1)

    ndays = vals.shape[0] // 24
    for d in range(ndays):
        chunk = np.ma.filled(vals[d*24:(d+1)*24, rs, cs], np.nan).astype('float32')
        chunk[:, ~box] = np.nan

        # hourly (mm/h)
        wet_h = chunk[chunk >= rain_min_h]
        h, _ = np.histogram(wet_h, bins=bins_h_edges)
        hist_h += h

        # daily accumulated (mm/d)
        daily = chunk.sum(axis=0)
        wet_d = daily[daily >= rain_min_d]
        h, _ = np.histogram(wet_d, bins=bins_d_edges)
        hist_d += h

    data.close()

    def to_pdf(hist, step):
        total = np.sum(hist)
        hist = hist.copy()
        hist[hist < 1] = np.nan
        return hist / (total * step)

    return to_pdf(hist_d, step_d), to_pdf(hist_h, step_h)


bins_d = (bins_d_edges[:-1] + bins_d_edges[1:]) / 2
bins_h = (bins_h_edges[:-1] + bins_h_edges[1:]) / 2

# Import model and obs dataset
dict_var = {'pr': ['pr', 'tp']}

imerg_pdf_d, imerg_pdf_h = compute_pdfs(dict_var[var][0], 'CSAM-3', 'GPM_IMERG')
era5_pdf_d, era5_pdf_h = compute_pdfs(dict_var[var][1], 'CSAM-3', 'ERA5')
sam12_pdf_d, sam12_pdf_h = compute_pdfs(dict_var[var][0], 'CSAM-3', 'ERA5_evaluation_RCM_RegCM5')
csam3_pdf_d, csam3_pdf_h = compute_pdfs(dict_var[var][0], 'CSAM-3', 'ERA5_evaluation_CPRCM_RegCM5')

# Datasets to plot: (label, color, daily pdf, hourly pdf)
plot_list = [
    ('Reanalysis (25 km)',           'blue',   era5_pdf_d,  era5_pdf_h),
    ('Satellite estimation (10 km)', 'black',  imerg_pdf_d, imerg_pdf_h),
    ('RCM (12 km)',  'red',    sam12_pdf_d, sam12_pdf_h),
    ('CPRCM (3 km)', 'green',  csam3_pdf_d, csam3_pdf_h),
]

# Plot figure
fig = plt.figure()
font_size = 10

ax = fig.add_subplot(1, 1, 1)
for label, color, pdf_d, pdf_h in plot_list:
    ax.plot(bins_d, pdf_d, marker='o', markersize=4, mfc=color, mec=color, alpha=0.70, linestyle='None', label=label)
ax.set_yscale('log')
ax.set_ylabel('Frequency (Log)', fontsize=font_size, fontweight='bold')
ax.set_xlabel('Precipitation (mm d$^{-1}$)', fontsize=font_size, fontweight='bold')
ax.grid(linestyle='--')
ax.legend(loc=3, ncol=1, fontsize=font_size-2, shadow=True)

# Hourly pdf inside the daily pdf plot
ax_in = ax.inset_axes([0.40, 0.50, 0.55, 0.45])
for label, color, pdf_d, pdf_h in plot_list:
    ax_in.plot(bins_h, pdf_h, marker='o', markersize=2, mfc=color, mec=color, alpha=0.70, linestyle='None', label=label)
ax_in.set_yscale('log')
ax_in.set_xlabel('Precipitation (mm h$^{-1}$)', fontsize=font_size-2, fontweight='bold')
ax_in.tick_params(labelsize=font_size-2)
ax_in.grid(linestyle='--')

# Path out to save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/TRACK-CYCLONE/CPRCM-TF/figs'
name_out = 'pyplt_test_graph_pdf_{0}_box_GPM_ERA5_RegCM5_{1}.png'.format(var, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
exit()
