# -*- coding: utf-8 -*-

# author       = "Leidinice Silva"
# email        = "leidinicesilva@gmail.com"
# date         = "Jul 28, 2026"
# description  = "This script plots 2x2 PDFs including hourly and daily precipitation"

import os
import glob
import netCDF4
import numpy as np
import matplotlib.pyplot as plt


def compute_pdf(data):
    step = 1
    rain_min = 1
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


def compute_pdf_(data):
    step = 0.1
    rain_min = 0
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

# Daily data 
data_cpc_d, _, _ = import_obs(path_data, "precip", "CPC", f"day_{yr}", "land_box")
x_cpc_d, pdf_cpc_d = compute_pdf(data_cpc_d)

data_era5_d, _, _ = import_obs(path_data, "tp", "ERA5", f"day_{yr}", "land_box")
x_era5_d, pdf_era5_d = compute_pdf(data_era5_d)

data_exp1_d, _, _ = import_rcm(path_data, "pr", "NoTo-EUR", f"day_{yr}", "land_box")
x_exp1_d, pdf_exp1_d = compute_pdf(data_exp1_d)

data_exp2_d, _, _ = import_rcm(path_data, "pr", "WSM5-EUR", f"day_{yr}", "land_box")
x_exp2_d, pdf_exp2_d = compute_pdf(data_exp2_d)

data_exp3_d, _, _ = import_rcm(path_data, "pr", "WSM7-EUR", f"day_{yr}", "land_box")
x_exp3_d, pdf_exp3_d = compute_pdf(data_exp3_d)

data_exp4_d, _, _ = import_rcm(path_data, "pr", "WDM7-EUR", f"day_{yr}", "land_box")
x_exp4_d, pdf_exp4_d = compute_pdf(data_exp4_d)

# Hourly data 
data_era5_h, _, _ = import_obs(path_data, "tp", "ERA5", f"1hr_{yr}", "land_box")
x_era5_h, pdf_era5_h = compute_pdf_(data_era5_h)

# Plotting 
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
font_size = 8

# Daily Precipitation
ax1.plot(x_cpc_d, pdf_cpc_d, marker='o', markersize=4, mfc='black', mec='black', alpha=0.70, linestyle='None', label='CPC')
ax1.plot(x_era5_d, pdf_era5_d, marker='o', markersize=4, mfc='gray', mec='gray', alpha=0.70, linestyle='None', label='ERA5')
ax1.plot(x_exp1_d, pdf_exp1_d, marker='o', markersize=4, mfc='red', mec='red', alpha=0.70, linestyle='None', label='NoTo')
ax1.plot(x_exp2_d, pdf_exp2_d, marker='o', markersize=4, mfc='blue', mec='blue', alpha=0.70, linestyle='None', label='WSM5')
ax1.plot(x_exp3_d, pdf_exp3_d, marker='o', markersize=4, mfc='green', mec='green', alpha=0.70, linestyle='None', label='WSM7')
ax1.plot(x_exp4_d, pdf_exp4_d, marker='o', markersize=4, mfc='orange', mec='orange', alpha=0.70, linestyle='None', label='WDM7')
ax1.set_yscale('log')
ax1.set_xlabel("Precipitation (mm/day)", fontsize=font_size)
ax1.set_ylabel("PDF", fontsize=font_size)
ax1.set_title("(a)", loc='left', fontsize=font_size, fontweight='bold')
ax1.grid(True, which="both", ls="--", alpha=0.5)
ax1.legend()

# Hourly Precipitation
ax2.plot(x_era5_h, pdf_era5_h, marker='o', markersize=4, mfc='gray', mec='gray', alpha=0.70, linestyle='None', label='ERA5')
ax2.set_yscale('log')
ax2.set_xlabel("Precipitation (mm/hr)", fontsize=font_size)
ax2.set_ylabel("PDF", fontsize=font_size)
ax2.set_title("(b)", loc='left', fontsize=font_size, fontweight='bold')
ax2.grid(True, which="both", ls="--", alpha=0.5)
ax2.legend()

# Save
name_out = f'pyplt_pdf_{var}_RegCM5_EUR-11_{yr}.png'
plt.tight_layout()
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.close()
exit()
