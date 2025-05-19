# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot clim maps"

import os
import netCDF4
import argparse
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat

from cartopy import config
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from import_climate_tools import compute_mbe

obs = 'ERA5'
dt = '1970-1971'
domain = 'SAM-22'
latlon = [-105, -16, -57, 18]

font_size = 8
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/{0}'.format(domain)


def main():

	parser = argparse.ArgumentParser(description='Process experiment configuration.')
	parser.add_argument('--var', type=str, required=True, help='Variable name (e.g., uv)')
	parser.add_argument('--level', type=str, required=True, help='Pressure level (e.g., 200hPa)')
	parser.add_argument('--exp_i', type=str, required=True, help='First experiment name (e.g., ctrl_RegCM5)')
	parser.add_argument('--exp_ii', type=str, required=True, help='Second experiment name (e.g., rclcrit_RegCM5)')

	args = parser.parse_args()

	var = args.var
	level = args.level
	exp_i = args.exp_i
	exp_ii = args.exp_ii

	exp_i_tg = exp_i.split('_RegCM5')[0]
	exp_i_up = exp_i_tg.upper()

	exp_ii_tg = exp_ii.split('_RegCM5')[0]
	exp_ii_up = exp_ii_tg.upper()

	dict_var = {var: ['u', 'v', 'q', 'ua', 'va', 'hus']}

	print(f"Variable: {var}")
	print(f"Level: {level}")
	print(f"exp_i_up: {exp_i_up}")
	print(f"exp_ii_up: {exp_ii_up}")

	return dict_var, var, level, exp_i, exp_i_tg, exp_i_up, exp_ii, exp_ii_tg, exp_ii_up  


if __name__ == '__main__':
	dict_var, var, level, exp_i, exp_i_tg, exp_i_up, exp_ii, exp_ii_tg, exp_ii_up = main()


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

lat, lon, q_exp_i_djf = import_rcm(dict_var[var][5], exp_i, 'DJF')
lat, lon, q_exp_i_mam = import_rcm(dict_var[var][5], exp_i, 'MAM')
lat, lon, q_exp_i_jja = import_rcm(dict_var[var][5], exp_i, 'JJA')
lat, lon, q_exp_i_son = import_rcm(dict_var[var][5], exp_i, 'SON')

lat, lon, u_exp_ii_djf = import_rcm(dict_var[var][3], exp_ii, 'DJF')
lat, lon, u_exp_ii_mam = import_rcm(dict_var[var][3], exp_ii, 'MAM')
lat, lon, u_exp_ii_jja = import_rcm(dict_var[var][3], exp_ii, 'JJA')
lat, lon, u_exp_ii_son = import_rcm(dict_var[var][3], exp_ii, 'SON')

lat, lon, v_exp_ii_djf = import_rcm(dict_var[var][4], exp_ii, 'DJF')
lat, lon, v_exp_ii_mam = import_rcm(dict_var[var][4], exp_ii, 'MAM')
lat, lon, v_exp_ii_jja = import_rcm(dict_var[var][4], exp_ii, 'JJA')
lat, lon, v_exp_ii_son = import_rcm(dict_var[var][4], exp_ii, 'SON')

lat, lon, q_exp_ii_djf = import_rcm(dict_var[var][5], exp_ii, 'DJF')
lat, lon, q_exp_ii_mam = import_rcm(dict_var[var][5], exp_ii, 'MAM')
lat, lon, q_exp_ii_jja = import_rcm(dict_var[var][5], exp_ii, 'JJA')
lat, lon, q_exp_ii_son = import_rcm(dict_var[var][5], exp_ii, 'SON')

mbe_djf_exp_i_obs_u = compute_mbe(u_exp_i_djf, u_obs_djf)
mbe_mam_exp_i_obs_u = compute_mbe(u_exp_i_mam, u_obs_mam)
mbe_jja_exp_i_obs_u = compute_mbe(u_exp_i_jja, u_obs_jja)
mbe_son_exp_i_obs_u = compute_mbe(u_exp_i_son, u_obs_son)

mbe_djf_exp_i_obs_v = compute_mbe(v_exp_i_djf, v_obs_djf)
mbe_mam_exp_i_obs_v = compute_mbe(v_exp_i_mam, v_obs_mam)
mbe_jja_exp_i_obs_v = compute_mbe(v_exp_i_jja, v_obs_jja)
mbe_son_exp_i_obs_v = compute_mbe(v_exp_i_son, v_obs_son)

mbe_djf_exp_ii_obs_u = compute_mbe(u_exp_ii_djf, u_obs_djf)
mbe_mam_exp_ii_obs_u = compute_mbe(u_exp_ii_mam, u_obs_mam)
mbe_jja_exp_ii_obs_u = compute_mbe(u_exp_ii_jja, u_obs_jja)
mbe_son_exp_ii_obs_u = compute_mbe(u_exp_ii_son, u_obs_son)

mbe_djf_exp_ii_obs_v = compute_mbe(v_exp_ii_djf, v_obs_djf)
mbe_mam_exp_ii_obs_v = compute_mbe(v_exp_ii_mam, v_obs_mam)
mbe_jja_exp_ii_obs_v = compute_mbe(v_exp_ii_jja, v_obs_jja)
mbe_son_exp_ii_obs_v = compute_mbe(v_exp_ii_son, v_obs_son)

mbe_djf_exp_i_obs_q = compute_mbe(q_exp_i_djf, q_obs_djf)
mbe_mam_exp_i_obs_q = compute_mbe(q_exp_i_mam, q_obs_mam)
mbe_jja_exp_i_obs_q = compute_mbe(q_exp_i_jja, q_obs_jja)
mbe_son_exp_i_obs_q = compute_mbe(q_exp_i_son, q_obs_son)

mbe_djf_exp_ii_obs_q = compute_mbe(q_exp_ii_djf, q_obs_djf)
mbe_mam_exp_ii_obs_q = compute_mbe(q_exp_ii_mam, q_obs_mam)
mbe_jja_exp_ii_obs_q = compute_mbe(q_exp_ii_jja, q_obs_jja)
mbe_son_exp_ii_obs_q = compute_mbe(q_exp_ii_son, q_obs_son)

uv_exp_i_djf = compute_ws(mbe_djf_exp_i_obs_u, mbe_djf_exp_i_obs_v)
uv_exp_i_mam = compute_ws(mbe_mam_exp_i_obs_u, mbe_mam_exp_i_obs_v)
uv_exp_i_jja = compute_ws(mbe_jja_exp_i_obs_u, mbe_jja_exp_i_obs_v)
uv_exp_i_son = compute_ws(mbe_son_exp_i_obs_u, mbe_son_exp_i_obs_v)

uv_exp_ii_djf = compute_ws(mbe_djf_exp_ii_obs_u, mbe_djf_exp_ii_obs_v)
uv_exp_ii_mam = compute_ws(mbe_mam_exp_ii_obs_u, mbe_mam_exp_ii_obs_v)
uv_exp_ii_jja = compute_ws(mbe_jja_exp_ii_obs_u, mbe_jja_exp_ii_obs_v)
uv_exp_ii_son = compute_ws(mbe_son_exp_ii_obs_u, mbe_son_exp_ii_obs_v)

# Plot figure
fig, axes = plt.subplots(2, 4, figsize=(12, 5), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()
vector = 25

if level == '200hPa':
	dict_plot={'uv': ['Bias of wind speed {0} (m s$^-$$^1$)'.format(level), np.arange(0, 31, 1), cm.viridis_r],
	'q': ['Bias of specific humidity {0} (g kg$^-$$^1$)'.format(level), np.arange(-0.02, 0.021, 0.001), cm.Spectral]}
elif level == '500hPa':
	dict_plot={'uv': ['Bias of wind speed {0} (m s$^-$$^1$)'.format(level), np.arange(0, 15.5, 0.5), cm.viridis_r],
	'q': ['Bias of specific humidity {0} (g kg$^-$$^1$)'.format(level), np.arange(-0.9, 0.91, 0.01), cm.Spectral]}
else:
	dict_plot={'uv': ['Bias of wind speed {0} (m s$^-$$^1$)'.format(level), np.arange(0, 7.5, 0.5), cm.viridis_r],
	'q': ['Bias of specific humidity {0} (g kg$^-$$^1$)'.format(level), np.arange(-3, 3.1, 0.1), cm.Spectral]}

plot_data = {'Plot 1': {'data 0': mbe_djf_exp_i_obs_u, 'data 1': mbe_djf_exp_i_obs_v, 'data 2': uv_exp_i_djf, 'data 3': mbe_djf_exp_i_obs_q, 'title': '(a) CTRL-{0} DJF'.format(obs)},
'Plot 2': {'data 0': mbe_mam_exp_i_obs_u, 'data 1': mbe_mam_exp_i_obs_v, 'data 2': uv_exp_i_mam, 'data 3': mbe_mam_exp_i_obs_q, 'title': '(b) {0}-{1} MAM'.format(exp_i_up, obs)},
'Plot 3': {'data 0': mbe_jja_exp_i_obs_u, 'data 1': mbe_jja_exp_i_obs_v, 'data 2': uv_exp_i_jja, 'data 3': mbe_jja_exp_i_obs_q, 'title': '(c) {0}-{1} JJA'.format(exp_i_up, obs)},
'Plot 4': {'data 0': mbe_son_exp_i_obs_u, 'data 1': mbe_son_exp_i_obs_v, 'data 2': uv_exp_i_son, 'data 3': mbe_son_exp_i_obs_q, 'title': '(d) {0}-{1} SON'.format(exp_i_up, obs)},
'Plot 5': {'data 0': mbe_djf_exp_ii_obs_u, 'data 1': mbe_djf_exp_ii_obs_v, 'data 2': uv_exp_ii_djf, 'data 3': mbe_djf_exp_ii_obs_q, 'title': '(e) {0}-{1} DJF'.format(exp_ii_up, obs)},
'Plot 6': {'data 0': mbe_mam_exp_ii_obs_u, 'data 1': mbe_mam_exp_ii_obs_v, 'data 2': uv_exp_ii_mam, 'data 3': mbe_mam_exp_ii_obs_q, 'title': '(f) {0}-{1} MAM'.format(exp_ii_up, obs)},
'Plot 7': {'data 0': mbe_jja_exp_ii_obs_u, 'data 1': mbe_jja_exp_ii_obs_v, 'data 2': uv_exp_ii_jja, 'data 3': mbe_jja_exp_ii_obs_q, 'title': '(g) {0}-{1} JJA'.format(exp_ii_up, obs)},
'Plot 8': {'data 0': mbe_son_exp_ii_obs_u, 'data 1': mbe_son_exp_ii_obs_v, 'data 2': uv_exp_ii_son, 'data 3': mbe_son_exp_ii_obs_q, 'title': '(h) {0}-{1} SON'.format(exp_ii_up, obs)}}

for ax, (key, value) in zip(axes, plot_data.items()):
	data_ = value['data 0']
	data_i = value['data 1']
	data_ii = value['data 2']
	data_iii = value['data 3']
	title = value['title']

	lons, lats = np.meshgrid(lon, lat)
	if var == 'uv':
		plt_map = ax.quiver(lons[::vector, ::vector], lats[::vector, ::vector], data_[::vector, ::vector], data_i[::vector, ::vector], data_ii[::vector, ::vector], cmap=dict_plot[var][2])     
	else:
		plt_map = ax.contourf(lon, lat, data_iii*1000, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	ax.set_title(title, loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax)   

# Set colobar
cbar = fig.colorbar(plt_map, ax=fig.axes, orientation='vertical', pad=0.025, aspect=50)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/figs/{1}'.format(path, exp_ii_tg)
name_out = 'pyplt_maps_bias_{0}_{1}_{2}_RegCM5_{3}.png'.format(var, level, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
exit()


