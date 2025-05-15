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

var = 'q'
level = '850hPa'
obs = 'ERA5'
dt = '1970-1971'
domain = 'SAM-22'
latlon = [-105, -16, -57, 18]

exp_i = 'ctrl_RegCM5'
exp_ii = 'vfqr_RegCM5'

font_size = 8
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/{0}'.format(domain)

dict_var = {var: ['u', 'v', 'q', 'ua', 'va', 'hus']}


def import_obs(param, dataset, season):

	arq   = '{0}/postproc/obs/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]

	if level == '200hPa':
		mean = var[:][0,22,:,:]
	elif level == '500hPa':
		mean = var[:][0,15,:,:]
	else:
		mean = var[:][0,6,:,:]

	return lat, lon, mean


def import_rcm(param, dataset, season):

	arq   = '{0}/postproc/rcm/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]

	if level == '200hPa':
		mean = var[:][0,9,:,:]
	elif level == '500hPa':
		mean = var[:][0,5,:,:]
	else:
		mean = var[:][0,2,:,:]

	return lat, lon, mean


def compute_ws(dataset_u, dataset_v):

	ws = np.sqrt(dataset_u**2+dataset_v**2)

	return ws


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
lat, lon, u_obs_djf = import_obs(dict_var[var][0], obs, 'DJF')
lat, lon, u_obs_mam = import_obs(dict_var[var][0], obs, 'MAM')
lat, lon, u_obs_jja = import_obs(dict_var[var][0], obs, 'JJA')
lat, lon, u_obs_son = import_obs(dict_var[var][0], obs, 'SON')

lat, lon, v_obs_djf = import_obs(dict_var[var][1], obs, 'DJF')
lat, lon, v_obs_mam = import_obs(dict_var[var][1], obs, 'MAM')
lat, lon, v_obs_jja = import_obs(dict_var[var][1], obs, 'JJA')
lat, lon, v_obs_son = import_obs(dict_var[var][1], obs, 'SON')

lat, lon, q_obs_djf = import_obs(dict_var[var][2], obs, 'DJF')
lat, lon, q_obs_mam = import_obs(dict_var[var][2], obs, 'MAM')
lat, lon, q_obs_jja = import_obs(dict_var[var][2], obs, 'JJA')
lat, lon, q_obs_son = import_obs(dict_var[var][2], obs, 'SON')

lat, lon, u_exp_i_djf = import_rcm(dict_var[var][3], exp_i, 'DJF')
lat, lon, u_exp_i_mam = import_rcm(dict_var[var][3], exp_i, 'MAM')
lat, lon, u_exp_i_jja = import_rcm(dict_var[var][3], exp_i, 'JJA')
lat, lon, u_exp_i_son = import_rcm(dict_var[var][3], exp_i, 'SON')

lat, lon, v_exp_i_djf = import_rcm(dict_var[var][4], exp_i, 'DJF')
lat, lon, v_exp_i_mam = import_rcm(dict_var[var][4], exp_i, 'MAM')
lat, lon, v_exp_i_jja = import_rcm(dict_var[var][4], exp_i, 'JJA')
lat, lon, v_exp_i_son = import_rcm(dict_var[var][4], exp_i, 'SON')

lat, lon, u_exp_ii_djf = import_rcm(dict_var[var][3], exp_ii, 'DJF')
lat, lon, u_exp_ii_mam = import_rcm(dict_var[var][3], exp_ii, 'MAM')
lat, lon, u_exp_ii_jja = import_rcm(dict_var[var][3], exp_ii, 'JJA')
lat, lon, u_exp_ii_son = import_rcm(dict_var[var][3], exp_ii, 'SON')

lat, lon, v_exp_ii_djf = import_rcm(dict_var[var][4], exp_ii, 'DJF')
lat, lon, v_exp_ii_mam = import_rcm(dict_var[var][4], exp_ii, 'MAM')
lat, lon, v_exp_ii_jja = import_rcm(dict_var[var][4], exp_ii, 'JJA')
lat, lon, v_exp_ii_son = import_rcm(dict_var[var][4], exp_ii, 'SON')

lat, lon, q_exp_i_djf = import_rcm(dict_var[var][5], exp_i, 'DJF')
lat, lon, q_exp_i_mam = import_rcm(dict_var[var][5], exp_i, 'MAM')
lat, lon, q_exp_i_jja = import_rcm(dict_var[var][5], exp_i, 'JJA')
lat, lon, q_exp_i_son = import_rcm(dict_var[var][5], exp_i, 'SON')

lat, lon, q_exp_ii_djf = import_rcm(dict_var[var][5], exp_ii, 'DJF')
lat, lon, q_exp_ii_mam = import_rcm(dict_var[var][5], exp_ii, 'MAM')
lat, lon, q_exp_ii_jja = import_rcm(dict_var[var][5], exp_ii, 'JJA')
lat, lon, q_exp_ii_son = import_rcm(dict_var[var][5], exp_ii, 'SON')

uv_obs_djf = compute_ws(u_obs_djf, v_obs_djf)
uv_obs_mam = compute_ws(u_obs_mam, v_obs_mam)
uv_obs_jja = compute_ws(u_obs_jja, v_obs_jja)
uv_obs_son = compute_ws(u_obs_son, v_obs_son)

uv_exp_i_djf = compute_ws(u_exp_i_djf, v_exp_i_djf)
uv_exp_i_mam = compute_ws(u_exp_i_mam, v_exp_i_mam)
uv_exp_i_jja = compute_ws(u_exp_i_jja, v_exp_i_jja)
uv_exp_i_son = compute_ws(u_exp_i_son, v_exp_i_son)

uv_exp_ii_djf = compute_ws(u_exp_ii_djf, v_exp_ii_djf)
uv_exp_ii_mam = compute_ws(u_exp_ii_mam, v_exp_ii_mam)
uv_exp_ii_jja = compute_ws(u_exp_ii_jja, v_exp_ii_jja)
uv_exp_ii_son = compute_ws(u_exp_ii_son, v_exp_ii_son)

# Plot figure
fig, axes = plt.subplots(3, 4, figsize=(12, 8), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()
font_size = 6 
vector = 25

if level == '200hPa':
	dict_plot={'uv': ['Wind speed {0} (m s$^-$$^1$)'.format(level), np.arange(0, 62, 2), cm.viridis_r],
	'q': ['Specific humidity {0} (g kg$^-$$^1$)'.format(level), np.arange(0, 1.1, 0.01), cm.magma_r]}
elif level == '500hPa':
	dict_plot={'uv': ['Wind speed {0} (m s$^-$$^1$)'.format(level), np.arange(0, 31, 1), cm.viridis_r],
	'q': ['Specific humidity {0} (g kg$^-$$^1$)'.format(level), np.arange(0, 11.25, 0.25), cm.magma_r]}
else:
	dict_plot={'uv': ['Wind speed {0} (m s$^-$$^1$)'.format(level), np.arange(0, 15.5, 0.5), cm.viridis_r],
	'q': ['Specific humidity {0} (g kg$^-$$^1$)'.format(level), np.arange(0, 22.5, 0.5), cm.hsv]}

plot_data = {'Plot 1': {'data 0': uv_obs_djf, 'data 1': u_obs_djf, 'data 2': v_obs_djf, 'data 3': q_obs_djf, 'title': '(a) {0} DJF'.format(obs)},
	'Plot 2': {'data 0': uv_obs_mam, 'data 1': u_obs_mam, 'data 2': v_obs_mam, 'data 3': q_obs_mam, 'title': '(b) {0} MAM'.format(obs)},
	'Plot 3': {'data 0': uv_obs_jja, 'data 1': u_obs_jja, 'data 2': v_obs_jja, 'data 3': q_obs_jja, 'title': '(c) {0} JJA'.format(obs)},
	'Plot 4': {'data 0': uv_obs_son, 'data 1': u_obs_son, 'data 2': v_obs_son, 'data 3': q_obs_son, 'title': '(d) {0} SON'.format(obs)},
	'Plot 5': {'data 0': uv_exp_i_djf, 'data 1': u_exp_i_djf, 'data 2': v_exp_i_djf, 'data 3': q_exp_i_djf, 'title': '(e) CTRL DJF'},
	'Plot 6': {'data 0': uv_exp_i_mam, 'data 1': u_exp_i_mam, 'data 2': v_exp_i_mam, 'data 3': q_exp_i_mam, 'title': '(f) CTRL MAM'},
	'Plot 7': {'data 0': uv_exp_i_jja, 'data 1': u_exp_i_jja, 'data 2': v_exp_i_jja, 'data 3': q_exp_i_jja, 'title': '(g) CTRL JJA'},
	'Plot 8': {'data 0': uv_exp_i_son, 'data 1': u_exp_i_son, 'data 2': v_exp_i_son, 'data 3': q_exp_i_son, 'title': '(h) CTRL SON'},
	'Plot 9': {'data 0': uv_exp_ii_djf, 'data 1': u_exp_ii_djf, 'data 2': v_exp_ii_djf, 'data 3': q_exp_ii_djf, 'title': '(c) VFQR DJF'},
	'Plot 10': {'data 0': uv_exp_ii_mam, 'data 1': u_exp_ii_mam, 'data 2': v_exp_ii_mam, 'data 3': q_exp_ii_mam, 'title': '(f) VFQR MAM'},
	'Plot 11': {'data 0': uv_exp_ii_jja, 'data 1': u_exp_ii_jja, 'data 2': v_exp_ii_jja, 'data 3': q_exp_ii_jja, 'title': '(i) VFQR JJA'},
	'Plot 12': {'data 0': uv_exp_ii_son, 'data 1': u_exp_ii_son, 'data 2': v_exp_ii_son, 'data 3': q_exp_ii_son, 'title': '(l) VFQR SON'}}

for ax, (key, value) in zip(axes, plot_data.items()):
	data_ = value['data 0']
	data_i = value['data 1']
	data_ii = value['data 2']
	data_iii = value['data 3']
	title = value['title']

	if var == 'ua':
		lons, lats = np.meshgrid(lon, lat)
		contourf = ax.contourf(lons, lats, data_, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
		quiv = ax.quiver(lons[::vector, ::vector], lats[::vector, ::vector], data_i[::vector, ::vector], data_ii[::vector, ::vector], color='black') 
		ax.set_title(title, loc='left', fontsize=font_size, fontweight='bold')
		configure_subplot(ax)     
	else:
		lons, lats = np.meshgrid(lon, lat)
		contourf = ax.contourf(lons, lats, data_iii*1000, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
		ax.set_title(title, loc='left', fontsize=font_size, fontweight='bold')
		configure_subplot(ax) 

# Set colobar
cbar = fig.colorbar(contourf, ax=fig.axes, orientation='vertical', pad=0.025, aspect=50)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/figs/vfqr'.format(path)
name_out = 'pyplt_maps_clim_{0}_{1}_{2}_RegCM5_{3}.png'.format(var, level, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


