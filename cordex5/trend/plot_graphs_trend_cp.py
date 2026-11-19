# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 02, 2026"
__description__ = "This script plot trend maps"

import os
import netCDF4
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

parser = argparse.ArgumentParser(description='Plot trend graphs')
parser.add_argument('--var', required=True, help='Variable name')
parser.add_argument('--domain', required=True, help='Domain name')
parser.add_argument('--idt', required=True, help='Initial year')
parser.add_argument('--fdt', required=True, help='Final year')
args = parser.parse_args()

# Configuration
var = args.var
domain = args.domain
idt = int(args.idt)
fdt = int(args.fdt)
years = np.arange(idt, fdt + 1)
dt = f"{idt}-{fdt}"
font_size = 10

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/postproc/trend'
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/figs/trend'

VAR_MAP = {
    'pr': 'tp', 
    'tas': 'tas',
    'tasmax': 'tasmax',
    'tasmin': 'tasmin'
}

UNITS = {
    'pr': 'Precipitation (mm yr$^{-1}$)',
    'tas': 'Air temperature (°C yr$^{-1}$)',
    'tasmax': 'Maximum air temperature (°C yr$^{-1}$)',
    'tasmin': 'Minimum air temperature (°C yr$^{-1}$)'
}

def import_simulation(param):

    filename = f"{path}/{domain}/{param}_{domain}_ERA5_RegCM5_year_{dt}.nc"
    with netCDF4.Dataset(filename) as nc:
        data = nc.variables[param][:]
        lat = nc.variables['lat'][:]
        lon = nc.variables['lon'][:]

    if data.shape[0] == 11:
        data = data[1:]

    print(f"  Data shape: {data.shape}")
    print(f"  Lon range: [{lon.min():.1f}, {lon.max():.1f}]")
    print(f"  Lat range: [{lat.min():.1f}, {lat.max():.1f}]")
    print()

    if lat.ndim == 1 and lon.ndim == 1:
        lon, lat = np.meshgrid(lon, lat)

    return lat, lon, data


def import_observation(param, lat_regcm, lon_regcm):

    var_name = VAR_MAP[param]
    filename = f"{path}/global/{var_name}_ERA5_year_{dt}.nc"

    with netCDF4.Dataset(filename) as nc:
        data = nc.variables[var_name][:]
        lat = nc.variables['latitude'][:]
        lon = nc.variables['longitude'][:]

    print(f"  Data shape: {data.shape}")
    print(f"  Lon range before: [{lon.min():.1f}, {lon.max():.1f}]")
    print(f"  Lat range: [{lat.min():.1f}, {lat.max():.1f}]")

    # fix lon
    if lon.max() > 180:
        lon = ((lon + 180) % 360) - 180
        idx = np.argsort(lon)
        lon = lon[idx]
        data = data[..., idx]

    if lat.ndim == 1 and lon.ndim == 1:
        lon, lat = np.meshgrid(lon, lat)

    print(f"  Lon range after: [{lon.min():.1f}, {lon.max():.1f}]")
    print()

    # crop
    lat_min, lat_max = np.nanmin(lat_regcm), np.nanmax(lat_regcm)
    lon_min, lon_max = np.nanmin(lon_regcm), np.nanmax(lon_regcm)

    mask = (
        (lat >= lat_min) & (lat <= lat_max) &
        (lon >= lon_min) & (lon <= lon_max)
    )

    idx = np.where(mask)
    i_min, i_max = idx[0].min(), idx[0].max()
    j_min, j_max = idx[1].min(), idx[1].max()

    data = data[:, i_min:i_max+1, j_min:j_max+1]

    return data


def compute_mean(data):
    return np.nanmean(data, axis=(1, 2))


def compute_trend(ts):
    slope, intercept, r, p, std = linregress(years, ts)
    trend = intercept + slope * years
    return slope, intercept, p, trend


def main():

    print(f"{var} {dt}")

    # importa datasets
    lat_regcm, lon_regcm, data_regcm = import_simulation(var)
    data_era5 = import_observation(var, lat_regcm, lon_regcm)

    # time series
    ts_regcm = compute_mean(data_regcm)
    ts_era5  = compute_mean(data_era5)

    # trends
    s_r, i_r, p_r, t_r = compute_trend(ts_regcm)
    s_e, i_e, p_e, t_e = compute_trend(ts_era5)

    print(f"RegCM: {s_r:.4f}/yr (p={p_r:.3f})")
    print(f"ERA5 : {s_e:.4f}/yr (p={p_e:.3f})")

    # plot
    plt.figure(figsize=(8,4))
    ax = plt.gca()
    ax.set_facecolor('#e6e6e6')

    plt.plot(years, ts_regcm, marker='.', linewidth=1, color='blue')
    plt.plot(years, t_r, '--', color='blue', label=f'RegCM5: {s_r:.4f}/yr (p={p_r:.4f})')
    plt.plot(years, ts_era5, marker='.', linewidth=1, color='red')
    plt.plot(years, t_e, '--', color='red', label=f'ERA5: {s_e:.4f}/yr (p={p_e:.4f})')
    plt.xlabel('Year', fontsize=font_size, fontweight='bold')
    plt.ylabel(f'{UNITS[var]}', fontsize=font_size, fontweight='bold')
    plt.grid(True, color='white', linestyle='--', linewidth=1)
    plt.legend(fontsize=font_size)

    # save figure
    os.makedirs(path_out, exist_ok=True)
    outfile = f"{path_out}/pyplt_graph_trend_{var}_{domain}_RegCM5_{dt}.png"
    plt.savefig(outfile, dpi=400, bbox_inches='tight')
    print(f"Saved: {outfile}")

if __name__ == "__main__":
    main()

