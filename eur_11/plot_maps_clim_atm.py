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

var = 'uv'
level = '200hPa'
domain = 'EUR-11'
dt = '2000-2001'
dataset = 'ERA5'
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11'


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
		mean = var[:][0,7,:,:]

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

	ax.set_xticks(np.arange(-40,65,25), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(30,85,15), crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.tick_params(axis='x', labelsize=6, labelcolor='black')
	ax.tick_params(axis='y', labelsize=6, labelcolor='black')
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.coastlines()


# Import model and obs dataset 
lat, lon, u_obs_djf = import_obs('u', dataset, 'DJF')
lat, lon, u_obs_mam = import_obs('u', dataset, 'MAM')
lat, lon, u_obs_jja = import_obs('u', dataset, 'JJA')
lat, lon, u_obs_son = import_obs('u', dataset, 'SON')

lat, lon, v_obs_djf = import_obs('v', dataset, 'DJF')
lat, lon, v_obs_mam = import_obs('v', dataset, 'MAM')
lat, lon, v_obs_jja = import_obs('v', dataset, 'JJA')
lat, lon, v_obs_son = import_obs('v', dataset, 'SON')

lat, lon, u_noto_djf = import_rcm('ua', 'NoTo-Europe_RegCM5', 'DJF')
lat, lon, u_noto_mam = import_rcm('ua', 'NoTo-Europe_RegCM5', 'MAM')
lat, lon, u_noto_jja = import_rcm('ua', 'NoTo-Europe_RegCM5', 'JJA')
lat, lon, u_noto_son = import_rcm('ua', 'NoTo-Europe_RegCM5', 'SON')

lat, lon, u_wdm7_djf = import_rcm('ua', 'WDM7-Europe_RegCM5', 'DJF')
lat, lon, u_wdm7_mam = import_rcm('ua', 'WDM7-Europe_RegCM5', 'MAM')
lat, lon, u_wdm7_jja = import_rcm('ua', 'WDM7-Europe_RegCM5', 'JJA')
lat, lon, u_wdm7_son = import_rcm('ua', 'WDM7-Europe_RegCM5', 'SON')

lat, lon, u_wsm7_djf = import_rcm('ua', 'WSM7-Europe_RegCM5', 'DJF')
lat, lon, u_wsm7_mam = import_rcm('ua', 'WSM7-Europe_RegCM5', 'MAM')
lat, lon, u_wsm7_jja = import_rcm('ua', 'WSM7-Europe_RegCM5', 'JJA')
lat, lon, u_wsm7_son = import_rcm('ua', 'WSM7-Europe_RegCM5', 'SON')

lat, lon, u_wsm5_djf = import_rcm('ua', 'WSM5-Europe_RegCM5', 'DJF')
lat, lon, u_wsm5_mam = import_rcm('ua', 'WSM5-Europe_RegCM5', 'MAM')
lat, lon, u_wsm5_jja = import_rcm('ua', 'WSM5-Europe_RegCM5', 'JJA')
lat, lon, u_wsm5_son = import_rcm('ua', 'WSM5-Europe_RegCM5', 'SON')

lat, lon, v_noto_djf = import_rcm('va', 'NoTo-Europe_RegCM5', 'DJF')
lat, lon, v_noto_mam = import_rcm('va', 'NoTo-Europe_RegCM5', 'MAM')
lat, lon, v_noto_jja = import_rcm('va', 'NoTo-Europe_RegCM5', 'JJA')
lat, lon, v_noto_son = import_rcm('va', 'NoTo-Europe_RegCM5', 'SON')

lat, lon, v_wdm7_djf = import_rcm('va', 'WDM7-Europe_RegCM5', 'DJF')
lat, lon, v_wdm7_mam = import_rcm('va', 'WDM7-Europe_RegCM5', 'MAM')
lat, lon, v_wdm7_jja = import_rcm('va', 'WDM7-Europe_RegCM5', 'JJA')
lat, lon, v_wdm7_son = import_rcm('va', 'WDM7-Europe_RegCM5', 'SON')

lat, lon, v_wsm7_djf = import_rcm('va', 'WSM7-Europe_RegCM5', 'DJF')
lat, lon, v_wsm7_mam = import_rcm('va', 'WSM7-Europe_RegCM5', 'MAM')
lat, lon, v_wsm7_jja = import_rcm('va', 'WSM7-Europe_RegCM5', 'JJA')
lat, lon, v_wsm7_son = import_rcm('va', 'WSM7-Europe_RegCM5', 'SON')

lat, lon, v_wsm5_djf = import_rcm('va', 'WSM5-Europe_RegCM5', 'DJF')
lat, lon, v_wsm5_mam = import_rcm('va', 'WSM5-Europe_RegCM5', 'MAM')
lat, lon, v_wsm5_jja = import_rcm('va', 'WSM5-Europe_RegCM5', 'JJA')
lat, lon, v_wsm5_son = import_rcm('va', 'WSM5-Europe_RegCM5', 'SON')

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
fig, axes = plt.subplots(4, 5, figsize=(15, 6), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()
font_size = 6 
vector = 40

if level == '200hPa':
	dict_plot={'uv': ['Wind speed {0} (m s$^-$$^1$)'.format(level), np.arange(0, 60, 5), cm.viridis_r]}
elif level == '500hPa':
	dict_plot={'uv': ['Wind speed {0} (m s$^-$$^1$)'.format(level), np.arange(0, 30, 2.5), cm.viridis_r]}
else:
	dict_plot={'uv': ['Wind speed {0} (m s$^-$$^1$)'.format(level), np.arange(0, 15, 1.25), cm.viridis_r]}

plot_data = {'Plot 1': {'data 0': uv_obs_djf, 'data 1': u_obs_djf, 'data 2': v_obs_djf, 'title': '(a) {0} DJF'.format(dataset)},
	'Plot 2': {'data 0': uv_noto_djf, 'data 1': u_noto_djf, 'data 2': v_noto_djf, 'title': '(b) NoTo DJF'},
	'Plot 3': {'data 0': uv_wdm7_djf, 'data 1': u_wdm7_djf, 'data 2': v_wdm7_djf, 'title': '(c) WDM7 DJF'},
	'Plot 4': {'data 0': uv_wsm7_djf, 'data 1': u_wsm7_djf, 'data 2': v_wsm7_djf, 'title': '(d) WSM7 DJF'},
	'Plot 5': {'data 0': uv_wsm5_djf, 'data 1': u_wsm5_djf, 'data 2': v_wsm5_djf, 'title': '(e) WSM5 DJF'},
	'Plot 6': {'data 0': uv_obs_mam, 'data 1': u_obs_mam, 'data 2': v_obs_mam, 'title': '(f) {0} MAM'.format(dataset)},
	'Plot 7': {'data 0': uv_noto_mam, 'data 1': u_noto_mam, 'data 2': v_noto_mam, 'title': '(g) NoTo MAM'},
	'Plot 8': {'data 0': uv_wdm7_mam, 'data 1': u_wdm7_mam, 'data 2': v_wdm7_mam, 'title': '(h) WDM7 MAM'},
	'Plot 9': {'data 0': uv_wsm7_mam, 'data 1': u_wsm7_mam, 'data 2': v_wsm7_mam, 'title': '(i) WSM7 MAM'},
	'Plot 10': {'data 0': uv_wsm5_mam, 'data 1': u_wsm5_mam, 'data 2': v_wsm5_mam, 'title': '(j) WSM5 MAM'},
	'Plot 11': {'data 0': uv_obs_jja, 'data 1': u_obs_jja, 'data 2': v_obs_jja, 'title': '(k) {0} JJA'.format(dataset)},
	'Plot 12': {'data 0': uv_noto_jja, 'data 1': u_noto_jja, 'data 2': v_noto_jja, 'title': '(l) NoTo JJA'},
	'Plot 13': {'data 0': uv_wdm7_jja, 'data 1': u_wdm7_jja, 'data 2': v_wdm7_jja, 'title': '(m) WDM7 JJA'},
	'Plot 14': {'data 0': uv_wsm7_jja, 'data 1': u_wsm7_jja, 'data 2': v_wsm7_jja, 'title': '(n) WSM7 JJA'},
	'Plot 15': {'data 0': uv_wsm5_jja, 'data 1': u_wsm5_jja, 'data 2': v_wsm5_jja, 'title': '(o) WSM5 JJA'},
	'Plot 16': {'data 0': uv_obs_son, 'data 1': u_obs_son, 'data 2': v_obs_son, 'title': '(p) {0} SON'.format(dataset)},
	'Plot 17': {'data 0': uv_noto_son, 'data 1': u_noto_son, 'data 2': v_noto_son, 'title': '(q) NoTo SON'},
	'Plot 18': {'data 0': uv_wdm7_son, 'data 1': u_wdm7_son, 'data 2': v_wdm7_son, 'title': '(r) WDM7 JJSONA'},
	'Plot 19': {'data 0': uv_wsm7_son, 'data 1': u_wsm7_son, 'data 2': v_wsm7_son, 'title': '(s) WSM7 SON'},
	'Plot 20': {'data 0': uv_wsm5_son, 'data 1': u_wsm5_son, 'data 2': v_wsm5_son, 'title': '(t) WSM5 SON'}}

for ax, (key, value) in zip(axes, plot_data.items()):
	data_ = value['data 0']
	data_i = value['data 1']
	data_ii = value['data 2']
	title = value['title']

	lons, lats = np.meshgrid(lon, lat)
	contour = ax.contourf(lons, lats, data_, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	quiv = ax.quiver(lons[::vector, ::vector], lats[::vector, ::vector], data_i[::vector, ::vector], data_ii[::vector, ::vector], color='white') 
	#stream = ax.streamplot(lons, lats, data_i, data_ii, color='white', density=0.25)
	ax.set_title(title, loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax)     

# Set colobar
cbar = fig.colorbar(contour, ax=fig.axes, orientation='vertical', pad=0.025, aspect=50)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/figs/totc'.format(path)
name_out = 'pyplt_maps_clim_{0}_{1}_{2}_RegCM5_{3}.png'.format(var, level, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


