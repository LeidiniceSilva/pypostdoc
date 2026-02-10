# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot bias maps"

import os
import string
import netCDF4
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat

from cartopy import config
from matplotlib.colors import LinearSegmentedColormap
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

var = 'pr'
obs = 'ERA5'
dt = '197001-02'
domain = 'SAM-12'
latlon = [-105, -16, -57, 18]

exp_i = 'RegCM5-ERA5_ICTP'
exp_ii = 'RegCM5-ERA5_USP'
exp_iii = 'RegCM5-ERA5_USP-dt60'
exp_iv = 'RegCM5-ERA5_USP-garoa'

dict_var = {'pr': ['tp']}

subdomains = {'AMZ': {'lat': (-15, -5), 'lon': (-68, -48)}, 
'LPB': {'lat': (-33, -20), 'lon': (-62, -49)},
'AND': {'lat': (-45, -35), 'lon': (-73, -69)},
'NEB': {'lat': (-15, -2), 'lon': (-45, -35)}}

subplot_labels = [f"({letter})" for letter in string.ascii_lowercase[:16]]

font_size = 8
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/{0}'.format(domain)


def pdf_rain_subdomain(value, lat, lon, subdomain):

	if subdomain in subdomains:
		bounds = subdomains[subdomain]
		lat_inds = np.where((lat >= bounds['lat'][0]) & (lat <= bounds['lat'][1]))[0]
		lon_inds = np.where((lon >= bounds['lon'][0]) & (lon <= bounds['lon'][1]))[0]
		value_ = value[:, lat_inds[:, None], lon_inds]

	step = 1
	rain_min = 0.1
	rain_max = 500.0
	bins = np.arange(rain_min, rain_max + step, step)

	rain_hist = np.zeros(len(bins) - 1)

	for t in range(value_.shape[0]):
		values = value_[t, :, :]
		wet = values[values >= rain_min]
		hist, _ = np.histogram(wet, bins=bins)
		rain_hist += hist

	total = np.sum(rain_hist)
	rain_hist[rain_hist < 1] = np.nan
	pdf = rain_hist / (total * step)

	return pdf


def import_obs(param, dataset, subdomain):

	arq   = '{0}/postproc/obs/{1}_{2}_{3}_day_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	data.close()

	pdf = pdf_rain_subdomain(value, lat, lon, subdomain)

	return pdf


def import_rcm(param, dataset, subdomain):

	arq   = '{0}/postproc/rcm/{1}_{2}_{3}_day_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	data.close()

	pdf = pdf_rain_subdomain(value, lat, lon, subdomain)

	return pdf


# Import model and obs dataset
regions = ['AMZ', 'LPB', 'AND', 'NEB']

obs_profiles = {f'obs_{region.lower()}': import_obs(dict_var[var][0], obs, region) for region in regions}
exp_i_profiles = {f'exp_i_{region.lower()}': import_rcm(var, exp_i, region) for region in regions}
exp_ii_profiles = {f'exp_ii_{region.lower()}': import_rcm(var, exp_ii, region) for region in regions}
exp_iii_profiles = {f'exp_iii_{region.lower()}': import_rcm(var, exp_iii, region) for region in regions}
exp_iv_profiles = {f'exp_iv_{region.lower()}': import_rcm(var, exp_iv, region) for region in regions}

# Plot figure  
fig, axes = plt.subplots(2, 2, figsize=(12, 12), sharex=True, sharey=True)

dict_plot = {'pr': ['Precipitation (mm d$^{-1}$)']}
subplot_labels = ['(a)', '(b)', '(c)', '(d)']

for idx, region in enumerate(regions):
	row = idx // 2
	col = idx % 2
	ax = axes[row][col]

	key = f'obs_{region.lower()}'
	key_i = f'exp_i_{region.lower()}'
	key_ii = f'exp_ii_{region.lower()}'
	key_iii = f'exp_iii_{region.lower()}'
	key_iv = f'exp_iv_{region.lower()}'

	obs_profile = obs_profiles[key]
	exp_i_profile = exp_i_profiles[key_i]
	exp_ii_profile = exp_ii_profiles[key_ii]
	exp_iii_profile = exp_iii_profiles[key_iii]
	exp_iv_profile = exp_iv_profiles[key_iv]

	ax.plot(obs_profile, marker='o', markersize=4, mfc='black', mec='black', alpha=0.70, linestyle='None', label='{0}'.format(obs))
	ax.plot(exp_i_profile, marker='o', markersize=4, mfc='blue', mec='black', alpha=0.70, linestyle='None', label='{0}'.format(exp_i))
	ax.plot(exp_ii_profile, marker='o', markersize=4, mfc='red', mec='black', alpha=0.70, linestyle='None', label='{0}'.format(exp_ii))
	ax.plot(exp_iii_profile, marker='o', markersize=4, mfc='green', mec='black', alpha=0.70, linestyle='None', label='{0}'.format(exp_iii))
	ax.plot(exp_iv_profile, marker='o', markersize=4, mfc='orange', mec='black', alpha=0.70, linestyle='None', label='{0}'.format(exp_iv))
	ax.set_title(f"{subplot_labels[idx]} {region}", loc='left', fontsize=font_size, fontweight='bold')

	if col == 0:
		ax.set_ylabel('Frequency (Log)', fontsize=font_size, fontweight='bold')
	if row == 1:
		ax.set_xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')

	ax.set_yscale('log')

	ax.grid(True, linestyle='--')

plt.legend(loc=1, ncol=1, fontsize=font_size)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_graph_pdf_{0}_{1}_RegCM5_{2}_v2.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
