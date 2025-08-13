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
parser.add_argument('--var', choices=['pr', 'tas', 'clt'], required=True)
args = parser.parse_args()
var = args.var

domain = 'EUR-11'
dt = '2000-2009'
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11'

if var == 'pr':
	dataset = 'EOBS'
else:
	dataset = 'ERA5'


def import_obs(param, dataset, season):

	arq   = '{0}/postproc/obs/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, dataset, domain, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

	return lat, lon, mean


def import_rcm(param, dataset, season):

	arq   = '{0}/postproc/rcm/{1}_{2}_{3}_{4}_lonlat.nc'.format(path, param, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]

	if param == 'tas':
		mean = var[:][0,0,:,:]
	else:
		mean = var[:][0,:,:]

	return lat, lon, mean
	

# Import model and obs dataset
dict_var = {'pr': ['rr', 'pre', 'precip', 'precipitation', 'tp'],
'tas': ['t2m'],
'clt': ['tcc']}

lat, lon, obs_djf = import_obs(dict_var[var][0], dataset, 'DJF')
lat, lon, obs_mam = import_obs(dict_var[var][0], dataset, 'MAM')
lat, lon, obs_jja = import_obs(dict_var[var][0], dataset, 'JJA')
lat, lon, obs_son = import_obs(dict_var[var][0], dataset, 'SON')

lat, lon, noto_djf = import_rcm(var, 'RegCM5_NoTo-EUR', 'DJF')
lat, lon, noto_mam = import_rcm(var, 'RegCM5_NoTo-EUR', 'MAM')
lat, lon, noto_jja = import_rcm(var, 'RegCM5_NoTo-EUR', 'JJA')
lat, lon, noto_son = import_rcm(var, 'RegCM5_NoTo-EUR', 'SON')

lat, lon, wdm7_djf = import_rcm(var, 'RegCM5_WDM7-EUR', 'DJF')
lat, lon, wdm7_mam = import_rcm(var, 'RegCM5_WDM7-EUR', 'MAM')
lat, lon, wdm7_jja = import_rcm(var, 'RegCM5_WDM7-EUR', 'JJA')
lat, lon, wdm7_son = import_rcm(var, 'RegCM5_WDM7-EUR', 'SON')

lat, lon, wsm7_djf = import_rcm(var, 'RegCM5_WSM7-EUR', 'DJF')
lat, lon, wsm7_mam = import_rcm(var, 'RegCM5_WSM7-EUR', 'MAM')
lat, lon, wsm7_jja = import_rcm(var, 'RegCM5_WSM7-EUR', 'JJA')
lat, lon, wsm7_son = import_rcm(var, 'RegCM5_WSM7-EUR', 'SON')

lat, lon, wsm5_djf = import_rcm(var, 'RegCM5_WSM5-EUR', 'DJF')
lat, lon, wsm5_mam = import_rcm(var, 'RegCM5_WSM5-EUR', 'MAM')
lat, lon, wsm5_jja = import_rcm(var, 'RegCM5_WSM5-EUR', 'JJA')
lat, lon, wsm5_son = import_rcm(var, 'RegCM5_WSM5-EUR', 'SON')

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

color=['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']

dict_plot={
'pr': ['Precipitation (mm d$^-$$^1$)', np.arange(0, 18.25, 0.25), matplotlib.colors.ListedColormap(color)],
'tas': ['Air temperature (°C)', np.arange(-20, 40, 2), cm.jet],
'clt': ['Total cloud cover (%)', np.arange(0, 1.05, 0.05), cm.Grays]}

plot_data = {'Plot 1': {'data': obs_djf, 'title': '(a) {0} DJF'.format(dataset)},
'Plot 2': {'data': noto_djf, 'title': '(b) NoTo DJF'},
'Plot 3': {'data': wdm7_djf, 'title': '(c) WDM7 DJF'},
'Plot 4': {'data': wsm7_djf, 'title': '(d) WSM7 DJF'},
'Plot 5': {'data': wsm5_djf, 'title': '(e) WSM5 DJF'},
'Plot 6': {'data': obs_mam, 'title': '(f) {0} MAM'.format(dataset)},
'Plot 7': {'data': noto_mam, 'title': '(g) NoTo MAM'},
'Plot 8': {'data': wdm7_mam, 'title': '(h) WDM7 MAM'},
'Plot 9': {'data': wsm7_mam, 'title': '(i) WSM7 MAM'},
'Plot 10': {'data': wsm5_mam, 'title': '(j) WSM5 MAM'},
'Plot 11': {'data': obs_jja, 'title': '(k) {0} JJA'.format(dataset)},
'Plot 12': {'data': noto_jja, 'title': '(l) NoTo JJA'},
'Plot 13': {'data': wdm7_jja, 'title': '(m) WDM7 JJA'},
'Plot 14': {'data': wsm7_jja, 'title': '(n) WSM7 JJA'},
'Plot 15': {'data': wsm5_jja, 'title': '(o) WSM5 JJA'},
'Plot 16': {'data': obs_son, 'title': '(p) {0} SON'.format(dataset)},
'Plot 17': {'data': noto_son, 'title': '(q) NoTo SON'},
'Plot 18': {'data': wdm7_son, 'title': '(r) WDM7 SON'},
'Plot 19': {'data': wsm7_son, 'title': '(s) WSM7 SON'},
'Plot 20': {'data': wsm5_son, 'title': '(t) WSM5 SON'}}

for ax, (key, value) in zip(axes, plot_data.items()):
	data = value['data']
	title = value['title']
    
	contour = ax.contourf(lon, lat, data, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	ax.set_title(title, loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax)

# Set colobar
cbar = fig.colorbar(contour, ax=fig.axes, orientation='vertical', pad=0.01, aspect=50)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_maps_clim_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
exit()


