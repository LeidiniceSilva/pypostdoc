# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot bias maps"

import os
import netCDF4
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat

from cartopy import config
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from import_climate_tools import compute_mbe

var = 'clt'
obs = 'ERA5'
dt = '1970-1971'
domain = 'SAM-22'
latlon = [-105, -16, -57, 18]

exp_i = 'ctrl_RegCM5'
exp_ii = 'vfqi_RegCM5'

font_size = 8
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/{0}'.format(domain)

dict_var = {'pr': ['tp'],
'tas': ['t2m'],
'clt': ['tcc']}


def import_obs(param, dataset, season):

	arq   = '{0}/postproc/obs/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:]
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

	return lat, lon, mean


def import_rcm(param, dataset, season):

	arq   = '{0}/postproc/rcm/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:]
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	
	if param == 'tas':
		mean = var[:][0,0,:,:]
	else:
		mean = var[:][0,:,:]
	
	return lat, lon, mean


def configure_subplot(ax):

	ax.set_extent(latlon, crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(latlon[0], latlon[1], 20), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(latlon[2], latlon[3], 20), crs=ccrs.PlateCarree())

	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(6)

	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.add_feature(cfeat.BORDERS, linewidth=0.5)
	ax.coastlines(linewidth=0.5)	


# Import model and obs dataset
lat, lon, obs_djf = import_obs(dict_var[var][0], obs, 'DJF')
lat, lon, obs_mam = import_obs(dict_var[var][0], obs, 'MAM')
lat, lon, obs_jja = import_obs(dict_var[var][0], obs, 'JJA')
lat, lon, obs_son = import_obs(dict_var[var][0], obs, 'SON')

lat, lon, exp_i_djf = import_rcm(var, exp_i, 'DJF')
lat, lon, exp_i_mam = import_rcm(var, exp_i, 'MAM')
lat, lon, exp_i_jja = import_rcm(var, exp_i, 'JJA')
lat, lon, exp_i_son = import_rcm(var, exp_i, 'SON')

lat, lon, exp_ii_djf = import_rcm(var, exp_ii, 'DJF')
lat, lon, exp_ii_mam = import_rcm(var, exp_ii, 'MAM')
lat, lon, exp_ii_jja = import_rcm(var, exp_ii, 'JJA')
lat, lon, exp_ii_son = import_rcm(var, exp_ii, 'SON')

mbe_djf_exp_i_obs = compute_mbe(exp_i_djf, obs_djf)
mbe_mam_exp_i_obs = compute_mbe(exp_i_mam, obs_mam)
mbe_jja_exp_i_obs = compute_mbe(exp_i_jja, obs_jja)
mbe_son_exp_i_obs = compute_mbe(exp_i_son, obs_son)

mbe_djf_exp_ii_obs = compute_mbe(exp_ii_djf, obs_djf)
mbe_mam_exp_ii_obs = compute_mbe(exp_ii_mam, obs_mam)
mbe_jja_exp_ii_obs = compute_mbe(exp_ii_jja, obs_jja)
mbe_son_exp_ii_obs = compute_mbe(exp_ii_son, obs_son)

# Plot figure
fig, axes = plt.subplots(2, 4, figsize=(12, 5), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()

dict_plot = {'pr': ['Bias of precipitation (mm d$^-$$^1$)', np.arange(-8, 8.5, 0.5), cm.BrBG],
'tas': ['Bias of air temperature (°C)', np.arange(-8, 8.5, 0.5), cm.RdBu_r],
'clt': ['Bias of total cloud cover (%)', np.arange(-60, 65, 5), cm.RdGy]}

plot_data = {'Plot 1': {'data': mbe_djf_exp_i_obs, 'title': '(a) CTRL-{0} DJF'.format(obs)},
'Plot 2': {'data': mbe_mam_exp_i_obs, 'title': '(b) CTRL-{0} MAM'.format(obs)},
'Plot 3': {'data': mbe_jja_exp_i_obs, 'title': '(c) CTRL-{0} JJA'.format(obs)},
'Plot 4': {'data': mbe_son_exp_i_obs, 'title': '(d) CTRL-{0} SON'.format(obs)},
'Plot 5': {'data': mbe_djf_exp_ii_obs, 'title': '(e) VFQI-{0} DJF'.format(obs)},
'Plot 6': {'data': mbe_mam_exp_ii_obs, 'title': '(f) VFQI-{0} MAM'.format(obs)},
'Plot 7': {'data': mbe_jja_exp_ii_obs, 'title': '(g) VFQI-{0} JJA'.format(obs)},
'Plot 8': {'data': mbe_son_exp_ii_obs, 'title': '(h) VFQI-{0} SON'.format(obs)}}

for ax, (key, value) in zip(axes, plot_data.items()):
	data = value['data']
	title = value['title']
    
	contour = ax.contourf(lon, lat, data, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	ax.set_title(title, loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax)

# Set colobar
cbar = fig.colorbar(contour, ax=fig.axes, orientation='vertical', pad=0.025, aspect=50)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

