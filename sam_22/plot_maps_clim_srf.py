# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot clim maps"

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

var = 'clt'
obs = 'ERA5'
dt = '1970-1971'
domain = 'SAM-22'
latlon = [-105, -16, -57, 18]

exp_i = 'ctrl_RegCM5'
exp_ii = 'vfqr_RegCM5'

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

# Plot figure
fig, axes = plt.subplots(3, 4, figsize=(12, 8), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()
font_size = 6 

color=['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']

dict_plot={
'pr': ['Precipitation (mm d$^-$$^1$)', np.arange(0, 18.5, 0.5), matplotlib.colors.ListedColormap(color)],
'tas': ['Air temperature (°C)', np.arange(-18, 39, 3), cm.jet],
'clt': ['Total cloud cover (%)', np.arange(0, 105, 5), cm.Greys]}


plot_data = {'Plot 1': {'data': obs_djf, 'title': '(a) {0} DJF'.format(obs)},
'Plot 2': {'data': obs_mam, 'title': '(b) {0} MAM'.format(obs)},
'Plot 3': {'data': obs_jja, 'title': '(c) {0} JJA'.format(obs)},
'Plot 4': {'data': obs_son, 'title': '(d) {0} SON'.format(obs)},
'Plot 5': {'data': exp_i_djf, 'title': '(e) CTRL DJF'},
'Plot 6': {'data': exp_i_mam, 'title': '(f) CTRL MAM'},
'Plot 7': {'data': exp_i_jja, 'title': '(g) CTRL JJA'},
'Plot 8': {'data': exp_i_son, 'title': '(h) CTRL SON'},
'Plot 9': {'data': exp_ii_djf, 'title': '(i) VFQR DJF'},
'Plot 10': {'data': exp_ii_mam, 'title': '(j) VFQR MAM'},
'Plot 11': {'data': exp_ii_jja, 'title': '(k) VFQR JJA'},
'Plot 12': {'data': exp_ii_son, 'title': '(l) VFQR SON'}}

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
path_out = '{0}/figs/vfqr'.format(path)
name_out = 'pyplt_maps_clim_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


