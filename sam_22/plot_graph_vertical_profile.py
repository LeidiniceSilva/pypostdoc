# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot vertical profile"

import os
import string
import netCDF4
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt

var = 'rh'
obs = 'ERA5'
dt = '1970-1971'
domain = 'SAM-22'
latlon = [-105, -16, -57, 18]

exp_i = 'ctrl_RegCM5'
exp_ii = 'vfqi_RegCM5'

font_size = 8
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/{0}'.format(domain)

dict_var = {'cl': ['cc'], 
'clw': ['clwc'],
'cli': ['ciwc'],
'rh': ['r'],
'hus': ['q']}

subdomains = {'AMZ': {'lat': (-15, -5), 'lon': (-68, -48)}, 
'LPB': {'lat': (-33, -20), 'lon': (-62, -49)},
'AND': {'lat': (-45, -35), 'lon': (-73, -69)},
'NEB': {'lat': (-15, -2), 'lon': (-45, -35)}}

subplot_labels = [f"({letter})" for letter in string.ascii_lowercase[:16]]

def import_obs(param, dataset, season, subdomain):

	arq  = '{0}/postproc/obs/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][0,:,:,:]
	
	if subdomain in subdomains:
		bounds = subdomains[subdomain]
		lat_inds = np.where((lat >= bounds['lat'][0]) & (lat <= bounds['lat'][1]))[0]
		lon_inds = np.where((lon >= bounds['lon'][0]) & (lon <= bounds['lon'][1]))[0]
		value = value[:, lat_inds[:, None], lon_inds]

	mean = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return mean


def import_rcm(param, dataset, season, subdomain):

	arq  = '{0}/postproc/rcm/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][0,:,:,:]
	
	if subdomain in subdomains:
		bounds = subdomains[subdomain]
		lat_inds = np.where((lat >= bounds['lat'][0]) & (lat <= bounds['lat'][1]))[0]
		lon_inds = np.where((lon >= bounds['lon'][0]) & (lon <= bounds['lon'][1]))[0]
		value = value[:, lat_inds[:, None], lon_inds]

	mean = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return mean
	
	
# Import model and obs dataset
seasons = ['DJF', 'MAM', 'JJA', 'SON']
regions = ['AMZ', 'LPB', 'AND', 'NEB']

obs_profiles = {f'obs_{season.lower()}_{region.lower()}': (import_obs(dict_var[var][0], obs, season, region) * (1 if var == 'rh' else 100 if var == 'cl' else 1000000))
    for season in seasons
    for region in regions}

exp_i_profiles = {f'exp_i_{season.lower()}_{region.lower()}': (import_rcm(var, exp_i, season, region) * (1 if var == 'rh' else 100 if var == 'cl' else 1000000))
    for season in seasons
    for region in regions}

exp_ii_profiles = {f'exp_ii_{season.lower()}_{region.lower()}': (import_rcm(var, exp_ii, season, region) * (1 if var == 'rh' else 100 if var == 'cl' else 1000000))
    for season in seasons
    for region in regions}

# Plot figure  
fig, axes = plt.subplots(4, 4, figsize=(16, 16), sharex=True, sharey=True)

dict_plot = {
'cl': ['Cloud fraction (%)', 0, 40, np.arange(0, 44, 4)],
'clw': ['Cloud liquid water (mg kg$^-$$^1$)', 0, 50, np.arange(0, 55, 5)],
'cli': ['Cloud liquid ice (mg kg$^-$$^1$)', 0, 20, np.arange(0, 22, 2)],
'rh': ['Relative humidity (%)', 0, 100, np.arange(0, 110, 10)]
}

levels_i = (1000,975,950,925,900,875,850,825,800,775,750,700,650,600,550,500,450,400,350,300,250,225,200,175,150,125,100,70,50,30,20,10,7,5,3,2,1)
levels_ii = (1000,925,850,700,600,500,400,300,250,200,150,100)

for i, season in enumerate(seasons):
	for j, region in enumerate(regions):
		key = f'obs_{season.lower()}_{region.lower()}'
		key_i = f'exp_i_{season.lower()}_{region.lower()}'
		key_ii = f'exp_ii_{season.lower()}_{region.lower()}'

		obs_profile = obs_profiles[key]
		exp_i_profile = exp_i_profiles[key_i]
		exp_ii_profile = exp_ii_profiles[key_ii]
	
		ax = axes[i, j]
		ax.plot(obs_profile, levels_i, color='black', label='ERA5', linewidth=1)
		ax.plot(exp_i_profile, levels_ii, color='blue', label='CTRL', linewidth=1)
		ax.plot(exp_ii_profile, levels_ii, color='red', label='VFQI', linewidth=1)
		ax.set_xlim(dict_plot[var][1], dict_plot[var][2])
		ax.set_ylim(0,1000)
		ax.set_xticks(dict_plot[var][3])
		plt.gca().invert_yaxis()

		subplot_idx = i * 4 + j
		index_label = subplot_labels[subplot_idx]
		ax.set_title(f"{index_label} {region} {season}", loc='left', fontsize=font_size, fontweight='bold')

		if j == 0:
			ax.set_ylabel('Pressure Level (hPa)', fontsize=font_size, fontweight='bold')
		if i == 3:
			ax.set_xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')

		ax.grid(True, linestyle='--')

plt.legend(loc=1, ncol=1, fontsize=font_size)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_graph_vertical_profile_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

