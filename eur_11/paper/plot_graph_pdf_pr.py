# -*- coding: utf-8 -*-

# author       = "Leidinice Silva"
# email        = "leidinicesilva@gmail.com"
# date         = "Jul 28, 2026"
# description  = "This script plots PDFs"

import os
import glob
import netCDF4
import numpy as np
import matplotlib.pyplot as plt


def compute_pdf(data):

    step = 1
    rain_min = 1.0
    rain_max = 500.0
    bins = np.arange(rain_min, rain_max + step, step)
    rain_hist = np.zeros(len(bins) - 1)
    for t in range(data.shape[0]):
        values = data[t, :, :]
        wet = values[values >= rain_min]
        hist, _ = np.histogram(wet, bins=bins)
        rain_hist += hist

    total = np.sum(rain_hist)
    rain_hist[rain_hist < 1] = np.nan  
    pdf = rain_hist / (total * step)

    bin_centers = (bins[:-1] + bins[1:] + step) / 2

    return bin_centers, pdf


def import_obs(path, param, obs_name, dt, tag):

    arq = f"{path}/{param}_{obs_name}_{dt}_{tag}.nc"
    files = glob.glob(arq)

    data = netCDF4.Dataset(files[0])
    var = np.squeeze(data.variables[param][:])

    if np.ma.isMaskedArray(var):
        var = var.filled(np.nan)

    lat_key = "latitude" if "latitude" in data.variables else "lat"
    lon_key = "longitude" if "longitude" in data.variables else "lon"

    lat = np.squeeze(data.variables[lat_key][:])
    lon = np.squeeze(data.variables[lon_key][:])

    data.close()

    return var, lat, lon


def import_rcm(path, param, exp_name, dt, tag):

    arq = f"{path}/{param}_RegCM5_{exp_name}_{dt}_{tag}.nc"
    files = glob.glob(arq)

    data = netCDF4.Dataset(files[0])
    var = np.squeeze(data.variables[param][:])

    if np.ma.isMaskedArray(var):
        var = var.filled(np.nan)

    lat_key = "latitude" if "latitude" in data.variables else "xlat"
    lon_key = "longitude" if "longitude" in data.variables else "xlon"

    lat = np.squeeze(data.variables[lat_key][:])
    lon = np.squeeze(data.variables[lon_key][:])

    data.close()

    return var, lat, lon


# Paths and Setup
path_data = "/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11/postproc/paper/pdf"
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11/figs/paper'
os.makedirs(path_out, exist_ok=True)

yr = "2000-2009"
var = "pr"

# Import
data_cpc, _, _ = import_obs(path_data, "precip", "CPC", f"day_{yr}", "land_box")
x_cpc, pdf_cpc = compute_pdf(data_cpc)

data_era5, _, _ = import_obs(path_data, "tp", "ERA5", f"day_{yr}", "land_box")
x_era5, pdf_era5 = compute_pdf(data_era5)

data_exp1, _, _ = import_rcm(path_data, "pr", "NoTo-EUR", f"day_{yr}", "land_box")
x_exp1, pdf_exp1 = compute_pdf(data_exp1)

data_exp2, _, _ = import_rcm(path_data, "pr", "WSM5-EUR", f"day_{yr}", "land_box")
x_exp2, pdf_exp2 = compute_pdf(data_exp2)

data_exp3, _, _ = import_rcm(path_data, "pr", "WSM7-EUR", f"day_{yr}", "land_box")
x_exp3, pdf_exp3 = compute_pdf(data_exp3)

data_exp4, _, _ = import_rcm(path_data, "pr", "WDM7-EUR", f"day_{yr}", "land_box")
x_exp4, pdf_exp4 = compute_pdf(data_exp4)

# Ploting
plt.figure()

plt.plot(x_cpc, pdf_cpc,  marker='o', markersize=4, mfc='black', mec='black', alpha=0.70, linestyle='None', label='CPC')
plt.plot(x_era5, pdf_era5, marker='o', markersize=4, mfc='gray', mec='gray', alpha=0.70, linestyle='None', label='ERA5')
plt.plot(x_exp1, pdf_exp1, marker='o', markersize=4, mfc='red', mec='red', alpha=0.70, linestyle='None', label='NoTo')
plt.plot(x_exp2, pdf_exp2, marker='o', markersize=4, mfc='blue', mec='blue', alpha=0.70, linestyle='None', label='WSM5')
plt.plot(x_exp3, pdf_exp3, marker='o', markersize=4, mfc='green', mec='green', alpha=0.70, linestyle='None', label='WSM7')
plt.plot(x_exp4, pdf_exp4, marker='o', markersize=4, mfc='orange', mec='orange', alpha=0.70, linestyle='None', label='WDM7')

plt.yscale('log')
plt.xlabel("Precipitation (mm/day)")
plt.ylabel("PDF")
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend()

# Save
name_out = f'pyplt_pdf_{var}_RegCM5_EUR-11_{yr}.png'
plt.tight_layout()
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.close()
exit()

