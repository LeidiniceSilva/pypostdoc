# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jan 25, 2026"
__description__ = "This script plot trend maps"

import os
import netCDF4
import argparse
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from scipy.stats import linregress

from cartopy import config
from matplotlib.colors import LinearSegmentedColormap
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

parser = argparse.ArgumentParser()
parser.add_argument('--var', required=True, help='Variable name')
parser.add_argument('--domain', required=True, help='Domain name')
parser.add_argument('--idt', required=True, help='Initial year')
parser.add_argument('--fdt', required=True, help='Final year')
args = parser.parse_args()

var = args.var
domain = args.domain
idt = args.idt
fdt = args.fdt
dt = '{0}-{1}'.format(idt, fdt)
font_size = 8

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5'


def import_dataset(param, dataset):

	folder = 'obs' if dataset in ['CRU', 'ERA5'] else 'rcm'
	arq = "{0}/postproc/evaluate/{1}/{2}_{3}_year_{4}.nc".format(path, folder, param, dataset, dt)

	data = netCDF4.Dataset(arq)
	lat_name = "latitude" if dataset == "ERA5" else "lat"
	lon_name = "longitude" if dataset == "ERA5" else "lon"
	var = data.variables[param][:]
	lat = data.variables[lat_name][:]
	lon = data.variables[lon_name][:]
	
	return lat, lon, var


def comp_trend_sig(data, min_valid=10, per_decade=False):
	
	nt, ny, nx = data.shape
	t = np.arange(nt)
	trend = np.full((ny, nx), np.nan)
	pval  = np.full((ny, nx), np.nan)

	for i in range(ny):
		for j in range(nx):
			y = data[:, i, j]
			mask = np.isfinite(y)

			if mask.sum() >= min_valid:
				res = linregress(t[mask], y[mask])
				trend[i, j] = res.slope
				pval[i, j]  = res.pvalue

	if per_decade:
		trend = trend * 10.0

	sig = np.ma.masked_where(pval >= 0.05, trend)

	return trend, sig


# Import model and obs dataset
dict_var = {'pr': ['pre', 'tp'],
'tas': ['tmp', 'tas'],
'tasmax': ['tmx', 'tasmax'],
'tasmin': ['tmn', 'tasmin']}

lat_cru, lon_cru, cru_yr = import_dataset(dict_var[var][0], 'CRU')
lat_era5, lon_era5, era5_yr = import_dataset(dict_var[var][1], 'ERA5')
lat_regcm, lon_regcm, regcm_yr = import_dataset(var, 'CSAM-3_RegCM5')

# Compute trend and significance
cru_trend, cru_sig = comp_trend_sig(cru_yr, min_valid=10, per_decade=False)
era5_trend, era5_sig = comp_trend_sig(era5_yr, min_valid=10, per_decade=False)
regcm_trend, regcm_sig = comp_trend_sig(regcm_yr, min_valid=10, per_decade=False)

# Plot figure
def configure_subplot(ax):

	lon_min = np.round(np.min(lon_regcm), 1)
	lon_max = np.round(np.max(lon_regcm), 1)
	lat_min = np.round(np.min(lat_regcm), 1)
	lat_max = np.round(np.max(lat_regcm), 1)
	ax.set_extent([np.min(lon_regcm), np.max(lon_regcm), np.min(lat_regcm), np.max(lat_regcm)], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(lon_min,lon_max,10), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(lat_min,lat_max,5), crs=ccrs.PlateCarree())

	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(6)

	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.add_feature(cfeat.BORDERS, linewidth=0.5)
	ax.coastlines(linewidth=0.5)

fig, axes = plt.subplots(1, 3, figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree()})

dict_plot = {'pr': ['Precipitation trend (mm yr$^-$$^1$) p ≥ 0.05', np.arange(-60, 62, 2), cm.BrBG],
'tas': ['Air temperature trend (°C yr$^-$$^1$) Dots: p ≥ 0.05', np.arange(-0.3, 0.31, 0.01), cm.bwr],
'tasmax': ['Maximum air temperature trend (°C yr$^-$$^1$) Dots: p ≥ 0.05', np.arange(-0.3, 0.31, 0.01), cm.bwr],
'tasmin': ['Minimum air temperature trend (°C yr$^-$$^1$) Dots: p ≥ 0.05', np.arange(-0.3, 0.31, 0.01), cm.bwr]}

datasets = [('CRU', lon_cru, lat_cru, cru_trend, cru_sig),
('ERA5', lon_era5, lat_era5, era5_trend, era5_sig),
('RegCM5', lon_regcm, lat_regcm, regcm_trend, regcm_sig),]

for ax, (name, lon, lat, trend, sig), lab in zip(
        axes, datasets, ['(a)', '(b)', '(c)']):

    plt_map = ax.contourf(lon, lat, trend, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both', transform=ccrs.PlateCarree())
    ax.contourf(lon, lat, sig, hatches=['...'], colors='none', transform=ccrs.PlateCarree())
    ax.set_title('{} {}'.format(lab, name), loc='left', fontsize=font_size, fontweight='bold')
    configure_subplot(ax)

# Set colobar
cbar = fig.colorbar(plt_map, ax=fig.axes, orientation='horizontal', pad=0.1, aspect=50)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/figs/evaluate'.format(path)
name_out = 'pyplt_maps_trend_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


