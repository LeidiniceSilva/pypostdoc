# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot diurnal cycle"

import os
import string
import netCDF4
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt

var = 'pr'
obs = 'ERA5'
dt = '1970-1971'
domain = 'SAM-22'
latlon = [-105, -16, -57, 18]

exp_i = 'ctrl_RegCM5'
exp_i_tg = exp_i.split('_RegCM5')[0]
exp_i_up = exp_i_tg.upper()

exp_ii = 'pbl_RegCM5'
exp_ii_tg = exp_ii.split('_RegCM5')[0]
exp_ii_up = exp_ii_tg.upper()

font_size = 8
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/{0}'.format(domain)

dict_var = {'pr': ['tp']}

subdomains = {'AMZ': {'lat': (-15, -5), 'lon': (-68, -48)}, 
'LPB': {'lat': (-33, -20), 'lon': (-62, -49)},
'AND': {'lat': (-45, -35), 'lon': (-73, -69)},
'NEB': {'lat': (-15, -2), 'lon': (-45, -35)}}

subplot_labels = [f"({letter})" for letter in string.ascii_lowercase[:16]]


def import_obs(param, dataset, subdomain):

	arq  = '{0}/postproc/obs/{1}_{2}_{3}_diurnal_cycle_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	if subdomain in subdomains:
		bounds = subdomains[subdomain]
		lat_inds = np.where((lat >= bounds['lat'][0]) & (lat <= bounds['lat'][1]))[0]
		lon_inds = np.where((lon >= bounds['lon'][0]) & (lon <= bounds['lon'][1]))[0]
		value = value[:, lat_inds[:, None], lon_inds]
	clim = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return clim


def import_rcm(param, dataset, subdomain):

	arq  = '{0}/postproc/rcm/{1}_{2}_{3}_diurnal_cycle_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	if subdomain in subdomains:
		bounds = subdomains[subdomain]
		lat_inds = np.where((lat >= bounds['lat'][0]) & (lat <= bounds['lat'][1]))[0]
		lon_inds = np.where((lon >= bounds['lon'][0]) & (lon <= bounds['lon'][1]))[0]
		value = value[:, lat_inds[:, None], lon_inds]
	mean = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return mean
	
	
# Import model and obs dataset
regions = ['AMZ', 'LPB', 'AND', 'NEB']

obs_profiles = {f'obs_{region.lower()}': import_obs(dict_var[var][0], obs, region) for region in regions}
exp_i_profiles = {f'exp_i_{region.lower()}': import_rcm(var, exp_i, region) for region in regions}
exp_ii_profiles = {f'exp_ii_{region.lower()}': import_rcm(var, exp_ii, region) for region in regions}

# Plot figure  
fig, axes = plt.subplots(2, 2, figsize=(12, 6), sharex=True, sharey=True)

dict_plot = {'pr': ['Precipitation (mm h$^{-1}$)', np.arange(0, 0.7, 0.01)]}
month_labels = ['00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', '']
subplot_labels = ['(a)', '(b)', '(c)', '(d)']
time = np.arange(0, 24) 

for idx, region in enumerate(regions):
	row = idx // 2
	col = idx % 2
	ax = axes[row][col]

	key = f'obs_{region.lower()}'
	key_i = f'exp_i_{region.lower()}'
	key_ii = f'exp_ii_{region.lower()}'

	obs_profile = obs_profiles[key]
	exp_i_profile = exp_i_profiles[key_i]
	exp_ii_profile = exp_ii_profiles[key_ii]

	ax.plot(time, obs_profile, color='black', label='{0}'.format(obs), linewidth=1., marker='.', linestyle='-')
	ax.plot(time, exp_i_profile, color='blue', label='{0}'.format(exp_i_up), linewidth=1., marker='.', linestyle='-')
	ax.plot(time, exp_ii_profile, color='red', label='{0}'.format(exp_ii_up), linewidth=1., marker='.', linestyle='-')
	#ax.set_yticks(dict_plot[var][1])
	ax.set_xticks(time)
	ax.set_xticklabels(month_labels, fontsize=font_size)
	ax.set_title(f"{subplot_labels[idx]} {region}", loc='left', fontsize=font_size, fontweight='bold')

	if col == 0:
		ax.set_ylabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')
	if row == 1:
		ax.set_xlabel('Hour', fontsize=font_size, fontweight='bold')

	ax.grid(True, linestyle='--')

plt.legend(loc=1, ncol=1, fontsize=font_size)

# Path out to save figure
path_out = '{0}/figs/{1}'.format(path, exp_ii_tg)
name_out = 'pyplt_graph_diurnal_cycle_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()



