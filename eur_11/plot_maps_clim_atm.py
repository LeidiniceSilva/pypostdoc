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

	if level == '200hPa':
		mean = var[:][0,22,:,:]
	elif level == '500hPa':
		mean = var[:][0,15,:,:]
	else:
		mean = var[:][0,7,:,:]

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

uv_obs_djf = compute_ws(u_obs_djf, v_obs_djf)
uv_obs_mam = compute_ws(u_obs_mam, v_obs_mam)
uv_obs_jja = compute_ws(u_obs_jja, v_obs_jja)
uv_obs_son = compute_ws(u_obs_son, v_obs_son)

uv_noto_djf = compute_ws(u_noto_djf, v_noto_djf)
uv_noto_mam = compute_ws(u_noto_mam, v_noto_mam)
uv_noto_jja = compute_ws(u_noto_jja, v_noto_jja)
uv_noto_son = compute_ws(u_noto_son, v_noto_son)

uv_wdm7_djf = compute_ws(u_wdm7_djf, v_wdm7_djf)
uv_wdm7_mam = compute_ws(u_wdm7_mam, v_wdm7_mam)
uv_wdm7_jja = compute_ws(u_wdm7_jja, v_wdm7_jja)
uv_wdm7_son = compute_ws(u_wdm7_son, v_wdm7_son)

uv_wsm7_djf = compute_ws(u_wsm7_djf, v_wsm7_djf)
uv_wsm7_mam = compute_ws(u_wsm7_mam, v_wsm7_mam)
uv_wsm7_jja = compute_ws(u_wsm7_jja, v_wsm7_jja)
uv_wsm7_son = compute_ws(u_wsm7_son, v_wsm7_son)

uv_wsm5_djf = compute_ws(u_wsm5_djf, v_wsm5_djf)
uv_wsm5_mam = compute_ws(u_wsm5_mam, v_wsm5_mam)
uv_wsm5_jja = compute_ws(u_wsm5_jja, v_wsm5_jja)
uv_wsm5_son = compute_ws(u_wsm5_son, v_wsm5_son)

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

fig, axes = plt.subplots(4, 5, figsize=(18, 8), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()
font_size = 6 
vector = 40

if level == '200hPa':
	dict_plot={'uv': ['Wind speed {0} (m s$^-$$^1$)'.format(level), np.arange(0, 60, 5), cm.viridis_r],
	'q': ['Specific humidity {0} (g kg$^-$$^1$)'.format(level), np.arange(0, 0.05, 0.0001), cm.hsv]}
elif level == '500hPa':
	dict_plot={'uv': ['Wind speed {0} (m s$^-$$^1$)'.format(level), np.arange(0, 30, 2.5), cm.viridis_r],
	'q': ['Specific humidity {0} (g kg$^-$$^1$)'.format(level), np.arange(0, 5.125, 0.125), cm.hsv]}
else:
	dict_plot={'uv': ['Wind speed {0} (m s$^-$$^1$)'.format(level), np.arange(0, 15, 1.25), cm.viridis_r],
	'q': ['Specific humidity {0} (g kg$^-$$^1$)'.format(level), np.arange(0, 20.25, 0.25), cm.hsv]}


plot_data = {'Plot 1': {'data 0': uv_obs_djf, 'data 1': u_obs_djf, 'data 2': v_obs_djf, 'data 3': q_obs_djf, 'title': '(a) {0} DJF'.format(dataset)},
	'Plot 2': {'data 0': uv_noto_djf, 'data 1': u_noto_djf, 'data 2': v_noto_djf, 'data 3': q_noto_djf, 'title': '(b) NoTo DJF'},
	'Plot 3': {'data 0': uv_wdm7_djf, 'data 1': u_wdm7_djf, 'data 2': v_wdm7_djf, 'data 3': q_wdm7_djf, 'title': '(c) WDM7 DJF'},
	'Plot 4': {'data 0': uv_wsm7_djf, 'data 1': u_wsm7_djf, 'data 2': v_wsm7_djf, 'data 3': q_wsm7_djf, 'title': '(d) WSM7 DJF'},
	'Plot 5': {'data 0': uv_wsm5_djf, 'data 1': u_wsm5_djf, 'data 2': v_wsm5_djf, 'data 3': q_wsm5_djf, 'title': '(e) WSM5 DJF'},
	'Plot 6': {'data 0': uv_obs_mam, 'data 1': u_obs_mam, 'data 2': v_obs_mam, 'data 3': q_obs_mam, 'title': '(f) {0} MAM'.format(dataset)},
	'Plot 7': {'data 0': uv_noto_mam, 'data 1': u_noto_mam, 'data 2': v_noto_mam, 'data 3': q_noto_mam, 'title': '(g) NoTo MAM'},
	'Plot 8': {'data 0': uv_wdm7_mam, 'data 1': u_wdm7_mam, 'data 2': v_wdm7_mam, 'data 3': q_wdm7_mam, 'title': '(h) WDM7 MAM'},
	'Plot 9': {'data 0': uv_wsm7_mam, 'data 1': u_wsm7_mam, 'data 2': v_wsm7_mam, 'data 3': q_wsm7_mam, 'title': '(i) WSM7 MAM'},
	'Plot 10': {'data 0': uv_wsm5_mam, 'data 1': u_wsm5_mam, 'data 2': v_wsm5_mam, 'data 3': q_wsm5_mam, 'title': '(j) WSM5 MAM'},
	'Plot 11': {'data 0': uv_obs_jja, 'data 1': u_obs_jja, 'data 2': v_obs_jja, 'data 3': q_obs_jja, 'title': '(k) {0} JJA'.format(dataset)},
	'Plot 12': {'data 0': uv_noto_jja, 'data 1': u_noto_jja, 'data 2': v_noto_jja, 'data 3': q_noto_jja, 'title': '(l) NoTo JJA'},
	'Plot 13': {'data 0': uv_wdm7_jja, 'data 1': u_wdm7_jja, 'data 2': v_wdm7_jja, 'data 3': q_wdm7_jja, 'title': '(m) WDM7 JJA'},
	'Plot 14': {'data 0': uv_wsm7_jja, 'data 1': u_wsm7_jja, 'data 2': v_wsm7_jja, 'data 3': q_wsm7_jja, 'title': '(n) WSM7 JJA'},
	'Plot 15': {'data 0': uv_wsm5_jja, 'data 1': u_wsm5_jja, 'data 2': v_wsm5_jja, 'data 3': q_wsm5_jja, 'title': '(o) WSM5 JJA'},
	'Plot 16': {'data 0': uv_obs_son, 'data 1': u_obs_son, 'data 2': v_obs_son, 'data 3': q_obs_son, 'title': '(p) {0} SON'.format(dataset)},
	'Plot 17': {'data 0': uv_noto_son, 'data 1': u_noto_son, 'data 2': v_noto_son, 'data 3': q_noto_son, 'title': '(q) NoTo SON'},
	'Plot 18': {'data 0': uv_wdm7_son, 'data 1': u_wdm7_son, 'data 2': v_wdm7_son, 'data 3': q_wdm7_son, 'title': '(r) WDM7 JJSONA'},
	'Plot 19': {'data 0': uv_wsm7_son, 'data 1': u_wsm7_son, 'data 2': v_wsm7_son, 'data 3': q_wsm7_son, 'title': '(s) WSM7 SON'},
	'Plot 20': {'data 0': uv_wsm5_son, 'data 1': u_wsm5_son, 'data 2': v_wsm5_son, 'data 3': q_wsm5_son, 'title': '(t) WSM5 SON'}}

for ax, (key, value) in zip(axes, plot_data.items()):
	data_ = value['data 0']
	data_i = value['data 1']
	data_ii = value['data 2']
	data_iii = value['data 3']
	title = value['title']

	lons, lats = np.meshgrid(lon, lat)
	if var == 'uv':
		contourf = ax.contourf(lons, lats, data_, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
		quiv = ax.quiver(lons[::vector, ::vector], lats[::vector, ::vector], data_i[::vector, ::vector], data_ii[::vector, ::vector], color='black')     
	else:
		contourf = ax.contourf(lons, lats, data_iii*1000, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	ax.set_title(title, loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax)    

# Set colobar
cbar = fig.colorbar(contourf, ax=fig.axes, orientation='vertical', pad=0.025, aspect=50)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_maps_clim_{0}_{1}_{2}_RegCM5_{3}.png'.format(var, level, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
exit()


