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
from import_climate_tools import compute_mbe
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

parser = argparse.ArgumentParser()
parser.add_argument('--var', choices=['uv', 'q'], required=True)
parser.add_argument('--level', choices=['850hPa', '500hPa', '200hPa'], required=True)
args = parser.parse_args()
var = args.var
level = args.level

dict_var = {var: ['u', 'v', 'q', 'ua', 'va', 'hus']}

domain = 'EUR-11'
dt = '2000-2009'
dataset = 'ERA5'
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11'


def import_obs(param, dataset, season):

	arq   = '{0}/postproc/obs/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, dataset, domain, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]

	if param == 'u' or param == 'v':
		if level == '200hPa':
			mean = var[:][0,22,:,:]
		elif level == '500hPa':
			mean = var[:][0,15,:,:]
		else:
			mean = var[:][0,6,:,:]
	else:
		if level == '200hPa':
			mean = var[:][0,14,:,:]
		elif level == '500hPa':
			mean = var[:][0,21,:,:]
		else:
			mean = var[:][0,30,:,:]

	return lat, lon, mean


def import_rcm(param, dataset, season):

	arq   = '{0}/postproc/rcm/{1}_{2}_{3}_{4}_lonlat.nc'.format(path, param, dataset, season, dt)	
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


# Import model and obs dataset 
lat, lon, u_obs_djf = import_obs(dict_var[var][0], dataset, 'DJF')
lat, lon, u_obs_mam = import_obs(dict_var[var][0], dataset, 'MAM')
lat, lon, u_obs_jja = import_obs(dict_var[var][0], dataset, 'JJA')
lat, lon, u_obs_son = import_obs(dict_var[var][0], dataset, 'SON')

lat, lon, v_obs_djf = import_obs(dict_var[var][1], dataset, 'DJF')
lat, lon, v_obs_mam = import_obs(dict_var[var][1], dataset, 'MAM')
lat, lon, v_obs_jja = import_obs(dict_var[var][1], dataset, 'JJA')
lat, lon, v_obs_son = import_obs(dict_var[var][1], dataset, 'SON')

lat, lon, q_obs_djf = import_obs(dict_var[var][2], dataset, 'DJF')
lat, lon, q_obs_mam = import_obs(dict_var[var][2], dataset, 'MAM')
lat, lon, q_obs_jja = import_obs(dict_var[var][2], dataset, 'JJA')
lat, lon, q_obs_son = import_obs(dict_var[var][2], dataset, 'SON')

lat, lon, u_noto_djf = import_rcm(dict_var[var][3], 'RegCM5_NoTo-EUR', 'DJF')
lat, lon, u_noto_mam = import_rcm(dict_var[var][3], 'RegCM5_NoTo-EUR', 'MAM')
lat, lon, u_noto_jja = import_rcm(dict_var[var][3], 'RegCM5_NoTo-EUR', 'JJA')
lat, lon, u_noto_son = import_rcm(dict_var[var][3], 'RegCM5_NoTo-EUR', 'SON')

lat, lon, u_wdm7_djf = import_rcm(dict_var[var][3], 'RegCM5_WDM7-EUR', 'DJF')
lat, lon, u_wdm7_mam = import_rcm(dict_var[var][3], 'RegCM5_WDM7-EUR', 'MAM')
lat, lon, u_wdm7_jja = import_rcm(dict_var[var][3], 'RegCM5_WDM7-EUR', 'JJA')
lat, lon, u_wdm7_son = import_rcm(dict_var[var][3], 'RegCM5_WDM7-EUR', 'SON')

lat, lon, u_wsm7_djf = import_rcm(dict_var[var][3], 'RegCM5_WSM7-EUR', 'DJF')
lat, lon, u_wsm7_mam = import_rcm(dict_var[var][3], 'RegCM5_WSM7-EUR', 'MAM')
lat, lon, u_wsm7_jja = import_rcm(dict_var[var][3], 'RegCM5_WSM7-EUR', 'JJA')
lat, lon, u_wsm7_son = import_rcm(dict_var[var][3], 'RegCM5_WSM7-EUR', 'SON')

lat, lon, u_wsm5_djf = import_rcm(dict_var[var][3], 'RegCM5_WSM5-EUR', 'DJF')
lat, lon, u_wsm5_mam = import_rcm(dict_var[var][3], 'RegCM5_WSM5-EUR', 'MAM')
lat, lon, u_wsm5_jja = import_rcm(dict_var[var][3], 'RegCM5_WSM5-EUR', 'JJA')
lat, lon, u_wsm5_son = import_rcm(dict_var[var][3], 'RegCM5_WSM5-EUR', 'SON')

lat, lon, v_noto_djf = import_rcm(dict_var[var][4], 'RegCM5_NoTo-EUR', 'DJF')
lat, lon, v_noto_mam = import_rcm(dict_var[var][4], 'RegCM5_NoTo-EUR', 'MAM')
lat, lon, v_noto_jja = import_rcm(dict_var[var][4], 'RegCM5_NoTo-EUR', 'JJA')
lat, lon, v_noto_son = import_rcm(dict_var[var][4], 'RegCM5_NoTo-EUR', 'SON')

lat, lon, v_wdm7_djf = import_rcm(dict_var[var][4], 'RegCM5_WDM7-EUR', 'DJF')
lat, lon, v_wdm7_mam = import_rcm(dict_var[var][4], 'RegCM5_WDM7-EUR', 'MAM')
lat, lon, v_wdm7_jja = import_rcm(dict_var[var][4], 'RegCM5_WDM7-EUR', 'JJA')
lat, lon, v_wdm7_son = import_rcm(dict_var[var][4], 'RegCM5_WDM7-EUR', 'SON')

lat, lon, v_wsm7_djf = import_rcm(dict_var[var][4], 'RegCM5_WSM7-EUR', 'DJF')
lat, lon, v_wsm7_mam = import_rcm(dict_var[var][4], 'RegCM5_WSM7-EUR', 'MAM')
lat, lon, v_wsm7_jja = import_rcm(dict_var[var][4], 'RegCM5_WSM7-EUR', 'JJA')
lat, lon, v_wsm7_son = import_rcm(dict_var[var][4], 'RegCM5_WSM7-EUR', 'SON')

lat, lon, v_wsm5_djf = import_rcm(dict_var[var][4], 'RegCM5_WSM5-EUR', 'DJF')
lat, lon, v_wsm5_mam = import_rcm(dict_var[var][4], 'RegCM5_WSM5-EUR', 'MAM')
lat, lon, v_wsm5_jja = import_rcm(dict_var[var][4], 'RegCM5_WSM5-EUR', 'JJA')
lat, lon, v_wsm5_son = import_rcm(dict_var[var][4], 'RegCM5_WSM5-EUR', 'SON')

lat, lon, q_noto_djf = import_rcm(dict_var[var][5], 'RegCM5_NoTo-EUR', 'DJF')
lat, lon, q_noto_mam = import_rcm(dict_var[var][5], 'RegCM5_NoTo-EUR', 'MAM')
lat, lon, q_noto_jja = import_rcm(dict_var[var][5], 'RegCM5_NoTo-EUR', 'JJA')
lat, lon, q_noto_son = import_rcm(dict_var[var][5], 'RegCM5_NoTo-EUR', 'SON')

lat, lon, q_wdm7_djf = import_rcm(dict_var[var][5], 'RegCM5_WDM7-EUR', 'DJF')
lat, lon, q_wdm7_mam = import_rcm(dict_var[var][5], 'RegCM5_WDM7-EUR', 'MAM')
lat, lon, q_wdm7_jja = import_rcm(dict_var[var][5], 'RegCM5_WDM7-EUR', 'JJA')
lat, lon, q_wdm7_son = import_rcm(dict_var[var][5], 'RegCM5_WDM7-EUR', 'SON')

lat, lon, q_wsm7_djf = import_rcm(dict_var[var][5], 'RegCM5_WSM7-EUR', 'DJF')
lat, lon, q_wsm7_mam = import_rcm(dict_var[var][5], 'RegCM5_WSM7-EUR', 'MAM')
lat, lon, q_wsm7_jja = import_rcm(dict_var[var][5], 'RegCM5_WSM7-EUR', 'JJA')
lat, lon, q_wsm7_son = import_rcm(dict_var[var][5], 'RegCM5_WSM7-EUR', 'SON')

lat, lon, q_wsm5_djf = import_rcm(dict_var[var][5], 'RegCM5_WSM5-EUR', 'DJF')
lat, lon, q_wsm5_mam = import_rcm(dict_var[var][5], 'RegCM5_WSM5-EUR', 'MAM')
lat, lon, q_wsm5_jja = import_rcm(dict_var[var][5], 'RegCM5_WSM5-EUR', 'JJA')
lat, lon, q_wsm5_son = import_rcm(dict_var[var][5], 'RegCM5_WSM5-EUR', 'SON')

u_mbe_djf_wdm7_obs = compute_mbe(u_wdm7_djf, u_obs_djf)
u_mbe_mam_wdm7_obs = compute_mbe(u_wdm7_mam, u_obs_mam)
u_mbe_jja_wdm7_obs = compute_mbe(u_wdm7_jja, u_obs_jja)
u_mbe_son_wdm7_obs = compute_mbe(u_wdm7_son, u_obs_son)

u_mbe_djf_wsm7_obs = compute_mbe(u_wsm7_djf, u_obs_djf)
u_mbe_mam_wsm7_obs = compute_mbe(u_wsm7_mam, u_obs_mam)
u_mbe_jja_wsm7_obs = compute_mbe(u_wsm7_jja, u_obs_jja)
u_mbe_son_wsm7_obs = compute_mbe(u_wsm7_son, u_obs_son)

u_mbe_djf_wsm5_obs = compute_mbe(u_wsm5_djf, u_obs_djf)
u_mbe_mam_wsm5_obs = compute_mbe(u_wsm5_mam, u_obs_mam)
u_mbe_jja_wsm5_obs = compute_mbe(u_wsm5_jja, u_obs_jja)
u_mbe_son_wsm5_obs = compute_mbe(u_wsm5_son, u_obs_son)

u_mbe_djf_noto_obs = compute_mbe(u_noto_djf, u_obs_djf)
u_mbe_mam_noto_obs = compute_mbe(u_noto_mam, u_obs_mam)
u_mbe_jja_noto_obs = compute_mbe(u_noto_jja, u_obs_jja)
u_mbe_son_noto_obs = compute_mbe(u_noto_son, u_obs_son)

v_mbe_djf_wdm7_obs = compute_mbe(v_wdm7_djf, v_obs_djf)
v_mbe_mam_wdm7_obs = compute_mbe(v_wdm7_mam, v_obs_mam)
v_mbe_jja_wdm7_obs = compute_mbe(v_wdm7_jja, v_obs_jja)
v_mbe_son_wdm7_obs = compute_mbe(v_wdm7_son, v_obs_son)

v_mbe_djf_wsm7_obs = compute_mbe(v_wsm7_djf, v_obs_djf)
v_mbe_mam_wsm7_obs = compute_mbe(v_wsm7_mam, v_obs_mam)
v_mbe_jja_wsm7_obs = compute_mbe(v_wsm7_jja, v_obs_jja)
v_mbe_son_wsm7_obs = compute_mbe(v_wsm7_son, v_obs_son)

v_mbe_djf_wsm5_obs = compute_mbe(v_wsm5_djf, v_obs_djf)
v_mbe_mam_wsm5_obs = compute_mbe(v_wsm5_mam, v_obs_mam)
v_mbe_jja_wsm5_obs = compute_mbe(v_wsm5_jja, v_obs_jja)
v_mbe_son_wsm5_obs = compute_mbe(v_wsm5_son, v_obs_son)

v_mbe_djf_noto_obs = compute_mbe(v_noto_djf, v_obs_djf)
v_mbe_mam_noto_obs = compute_mbe(v_noto_mam, v_obs_mam)
v_mbe_jja_noto_obs = compute_mbe(v_noto_jja, v_obs_jja)
v_mbe_son_noto_obs = compute_mbe(v_noto_son, v_obs_son)

uv_mbe_djf_wdm7_obs = compute_ws(u_mbe_djf_wdm7_obs, v_mbe_djf_wdm7_obs)
uv_mbe_mam_wdm7_obs = compute_ws(u_mbe_mam_wdm7_obs, v_mbe_mam_wdm7_obs)
uv_mbe_jja_wdm7_obs = compute_ws(u_mbe_jja_wdm7_obs, v_mbe_jja_wdm7_obs)
uv_mbe_son_wdm7_obs = compute_ws(u_mbe_son_wdm7_obs, v_mbe_son_wdm7_obs)

uv_mbe_djf_wsm7_obs = compute_ws(u_mbe_djf_wsm7_obs, v_mbe_djf_wsm7_obs)
uv_mbe_mam_wsm7_obs = compute_ws(u_mbe_mam_wsm7_obs, v_mbe_mam_wsm7_obs)
uv_mbe_jja_wsm7_obs = compute_ws(u_mbe_jja_wsm7_obs, v_mbe_jja_wsm7_obs)
uv_mbe_son_wsm7_obs = compute_ws(u_mbe_son_wsm7_obs, v_mbe_son_wsm7_obs)

uv_mbe_djf_wsm5_obs = compute_ws(u_mbe_djf_wsm5_obs, v_mbe_djf_wsm5_obs)
uv_mbe_mam_wsm5_obs = compute_ws(u_mbe_mam_wsm5_obs, v_mbe_mam_wsm5_obs)
uv_mbe_jja_wsm5_obs = compute_ws(u_mbe_jja_wsm5_obs, v_mbe_jja_wsm5_obs)
uv_mbe_son_wsm5_obs = compute_ws(u_mbe_son_wsm5_obs, v_mbe_son_wsm5_obs)

uv_mbe_djf_noto_obs = compute_ws(u_mbe_djf_noto_obs, v_mbe_djf_noto_obs)
uv_mbe_mam_noto_obs = compute_ws(u_mbe_mam_noto_obs, v_mbe_mam_noto_obs)
uv_mbe_jja_noto_obs = compute_ws(u_mbe_jja_noto_obs, v_mbe_jja_noto_obs)
uv_mbe_son_noto_obs = compute_ws(u_mbe_son_noto_obs, v_mbe_son_noto_obs)

q_mbe_djf_noto_obs = compute_mbe(q_noto_djf, q_obs_djf)
q_mbe_mam_noto_obs = compute_mbe(q_noto_mam, q_obs_mam)
q_mbe_jja_noto_obs = compute_mbe(q_noto_jja, q_obs_jja)
q_mbe_son_noto_obs = compute_mbe(q_noto_son, q_obs_son)

q_mbe_djf_wdm7_obs = compute_mbe(q_wdm7_djf, q_obs_djf)
q_mbe_mam_wdm7_obs = compute_mbe(q_wdm7_mam, q_obs_mam)
q_mbe_jja_wdm7_obs = compute_mbe(q_wdm7_jja, q_obs_jja)
q_mbe_son_wdm7_obs = compute_mbe(q_wdm7_son, q_obs_son)

q_mbe_djf_wsm7_obs = compute_mbe(q_wsm7_djf, q_obs_djf)
q_mbe_mam_wsm7_obs = compute_mbe(q_wsm7_mam, q_obs_mam)
q_mbe_jja_wsm7_obs = compute_mbe(q_wsm7_jja, q_obs_jja)
q_mbe_son_wsm7_obs = compute_mbe(q_wsm7_son, q_obs_son)

q_mbe_djf_wsm5_obs = compute_mbe(q_wsm5_djf, q_obs_djf)
q_mbe_mam_wsm5_obs = compute_mbe(q_wsm5_mam, q_obs_mam)
q_mbe_jja_wsm5_obs = compute_mbe(q_wsm5_jja, q_obs_jja)
q_mbe_son_wsm5_obs = compute_mbe(q_wsm5_son, q_obs_son)

# Plot figure   
def configure_subplot(ax):

	lon_min = np.round(np.min(lon), 1)
	lon_max = np.round(np.max(lon), 1)
	lat_min = np.round(np.min(lat), 1)
	lat_max = np.round(np.max(lat), 1)
	ax.set_extent([np.min(lon), np.max(lon), np.min(lat), np.max(lat)], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(lon_min,lon_max,20), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(lat_min,lat_max,15), crs=ccrs.PlateCarree())

	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(6)

	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.coastlines(linewidth=0.6)

fig, axes = plt.subplots(4, 4, figsize=(18, 8), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()
font_size = 6 
vector = 40

if level == '200hPa':
	dict_plot={'uv': ['Bias of wind speed {0} (m s$^-$$^1$)'.format(level), np.arange(-4, 4.25, 0.25), cm.Spectral],
	'q': ['Bias of specific humidity {0} (g kg$^-$$^1$)'.format(level), np.arange(-0.004, 0.0044, 0.0004), cm.twilight_shifted]}
elif level == '500hPa':
	dict_plot={'uv': ['Bias of wind speed {0} (m s$^-$$^1$)'.format(level), np.arange(-6, 6.5, 0.5), cm.Spectral],
	'q': ['Bias of specific humidity {0} (g kg$^-$$^1$)'.format(level), np.arange(-0.6, 0.66, 0.06), cm.twilight_shifted]}
else:
	dict_plot={'uv': ['Bias of wind speed {0} (m s$^-$$^1$)'.format(level), np.arange(-6, 6.5, 0.5), cm.Spectral],
	'q': ['Bias of specific humidity {0} (g kg$^-$$^1$)'.format(level), np.arange(-2, 2.2, 0.2), cm.twilight_shifted]}

plot_data = {'Plot 1': {'data 0': uv_mbe_djf_noto_obs, 'data 1': q_mbe_djf_noto_obs, 'title': '(a) NoTo-{0} DJF'.format(dataset)},
'Plot 2': {'data 0': uv_mbe_djf_wdm7_obs, 'data 1': q_mbe_djf_wdm7_obs, 'title': '(b) WDM7-{0} DJF'.format(dataset)},
'Plot 3': {'data 0': uv_mbe_djf_wsm7_obs, 'data 1': q_mbe_djf_wsm7_obs, 'title': '(c) WSM7-{0} DJF'.format(dataset)},
'Plot 4': {'data 0': uv_mbe_djf_wsm5_obs, 'data 1': q_mbe_djf_wsm5_obs, 'title': '(d) WSM5-{0} DJF'.format(dataset)},
'Plot 5': {'data 0': uv_mbe_mam_noto_obs, 'data 1': q_mbe_mam_noto_obs, 'title': '(e) NoTo-{0} MAM'.format(dataset)},
'Plot 6': {'data 0': uv_mbe_mam_wdm7_obs, 'data 1': q_mbe_mam_wdm7_obs, 'title': '(f) WDM7-{0} MAM'.format(dataset)},
'Plot 7': {'data 0': uv_mbe_mam_wsm7_obs, 'data 1': q_mbe_mam_wsm7_obs, 'title': '(g) WSM7-{0} MAM'.format(dataset)},
'Plot 8': {'data 0': uv_mbe_mam_wsm5_obs, 'data 1': q_mbe_mam_wsm5_obs, 'title': '(h) WSM5-{0} MAM'.format(dataset)},
'Plot 9': {'data 0': uv_mbe_jja_noto_obs, 'data 1': q_mbe_jja_noto_obs, 'title': '(i) NoTo-{0} JJA'.format(dataset)},
'Plot 10': {'data 0': uv_mbe_jja_wdm7_obs, 'data 1': q_mbe_jja_wdm7_obs, 'title': '(j) WDM7-{0} JJA'.format(dataset)},
'Plot 11': {'data 0': uv_mbe_jja_wsm7_obs, 'data 1': q_mbe_jja_wsm7_obs, 'title': '(k) WSM7-{0} JJA'.format(dataset)},
'Plot 12': {'data 0': uv_mbe_jja_wsm5_obs, 'data 1': q_mbe_jja_wsm5_obs, 'title': '(l) WSM5-{0} JJA'.format(dataset)},
'Plot 13': {'data 0': uv_mbe_son_noto_obs, 'data 1': q_mbe_son_noto_obs, 'title': '(m) NoTo-{0} SON'.format(dataset)},
'Plot 14': {'data 0': uv_mbe_son_wdm7_obs, 'data 1': q_mbe_son_wdm7_obs, 'title': '(n) WDM7-{0} SON'.format(dataset)},
'Plot 15': {'data 0': uv_mbe_son_wsm7_obs, 'data 1': q_mbe_son_wsm7_obs, 'title': '(o) WSM7-{0} SON'.format(dataset)},
'Plot 16': {'data 0': uv_mbe_son_wsm5_obs, 'data 1': q_mbe_son_wsm5_obs, 'title': '(p) WSM5-{0} SON'.format(dataset)}}

for ax, (key, value) in zip(axes, plot_data.items()):
	data_i = value['data 0']
	data_ii = value['data 1']
	title = value['title']

	lons, lats = np.meshgrid(lon, lat)
	if var == 'uv':
		contourf = ax.contourf(lons, lats, data_i, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	else:
		contourf = ax.contourf(lons, lats, data_ii*1000, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	ax.set_title(title, loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax)    

# Set colobar
cbar = fig.colorbar(contourf, ax=fig.axes, orientation='vertical', pad=0.025, aspect=50)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_{2}_RegCM5_{3}.png'.format(var, level, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
exit()







