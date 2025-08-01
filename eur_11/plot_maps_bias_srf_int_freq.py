# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot bias maps"

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

parser = argparse.ArgumentParser()
parser.add_argument('--stats', choices=['int', 'freq'], required=True)
args = parser.parse_args()
stats = args.stats

var = 'pr'
dt = '2000-2009'
domain = 'EUR-11'
dataset = 'EOBS'
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11'


def import_obs(param, dataset, season):

	arq   = '{0}/postproc/obs/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, param, stats, dataset, domain, season, dt)	

	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	
		
def import_rcm(param, dataset, season):

	arq   = '{0}/postproc/rcm/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, stats, dataset, season, dt)		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean


# Import model and obs dataset
dict_var = {'pr': ['precip', 'rr', 'tp']}

lat, lon, obs_djf = import_obs(dict_var[var][1], dataset, 'DJF')
lat, lon, obs_mam = import_obs(dict_var[var][1], dataset, 'MAM')
lat, lon, obs_jja = import_obs(dict_var[var][1], dataset, 'JJA')
lat, lon, obs_son = import_obs(dict_var[var][1], dataset, 'SON')

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

mbe_djf_noto_obs = compute_mbe(noto_djf, obs_djf)	
mbe_mam_noto_obs = compute_mbe(noto_mam, obs_mam)	
mbe_jja_noto_obs = compute_mbe(noto_jja, obs_jja)	
mbe_son_noto_obs = compute_mbe(noto_son, obs_son)  
	
mbe_djf_wdm7_obs = compute_mbe(wdm7_djf, obs_djf)	
mbe_mam_wdm7_obs = compute_mbe(wdm7_mam, obs_mam)	
mbe_jja_wdm7_obs = compute_mbe(wdm7_jja, obs_jja)	
mbe_son_wdm7_obs = compute_mbe(wdm7_son, obs_son) 

mbe_djf_wsm7_obs = compute_mbe(wsm7_djf, obs_djf)	
mbe_mam_wsm7_obs = compute_mbe(wsm7_mam, obs_mam)	
mbe_jja_wsm7_obs = compute_mbe(wsm7_jja, obs_jja)	
mbe_son_wsm7_obs = compute_mbe(wsm7_son, obs_son) 

mbe_djf_wsm5_obs = compute_mbe(wsm5_djf, obs_djf)	
mbe_mam_wsm5_obs = compute_mbe(wsm5_mam, obs_mam)	
mbe_jja_wsm5_obs = compute_mbe(wsm5_jja, obs_jja)	
mbe_son_wsm5_obs = compute_mbe(wsm5_son, obs_son)   

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

fig, axes = plt.subplots(4, 4, figsize=(15, 7), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()
font_size = 8

if stats == 'int':
	dict_plot = {'pr': ['Intensity of daily precipitation (mm d$^-$$^1$)', np.arange(-6, 6.5, 0.5), cm.BrBG]}

else:
	dict_plot = {'pr': ['Frequency of daily precipitation (%)', np.arange(-50, 55, 5), cm.BrBG]}

plot_data = {'Plot 1': {'data': mbe_djf_noto_obs[0][0], 'title': '(a) NoTo-{0} DJF'.format(dataset)},
'Plot 2': {'data': mbe_djf_wdm7_obs[0][0], 'title': '(b) WDM7-{0} DJF'.format(dataset)},
'Plot 3': {'data': mbe_djf_wsm7_obs[0][0], 'title': '(c) WSM7-{0} DJF'.format(dataset)},
'Plot 4': {'data': mbe_djf_wsm5_obs[0][0], 'title': '(d) WSM5-{0} DJF'.format(dataset)},
'Plot 5': {'data': mbe_mam_noto_obs[0][0], 'title': '(e) NoTo-{0} MAM'.format(dataset)},
'Plot 6': {'data': mbe_mam_wdm7_obs[0][0], 'title': '(f) WDM7-{0} MAM'.format(dataset)},
'Plot 7': {'data': mbe_mam_wsm7_obs[0][0], 'title': '(g) WSM7-{0} MAM'.format(dataset)},
'Plot 8': {'data': mbe_mam_wsm5_obs[0][0], 'title': '(h) WSM5-{0} MAM'.format(dataset)},
'Plot 9': {'data': mbe_jja_noto_obs[0][0], 'title': '(i) NoTo-{0} JJA'.format(dataset)},
'Plot 10': {'data': mbe_jja_wdm7_obs[0][0], 'title': '(j) WDM7-{0} JJA'.format(dataset)},
'Plot 11': {'data': mbe_jja_wsm7_obs[0][0], 'title': '(k) WSM7-{0} JJA'.format(dataset)},
'Plot 12': {'data': mbe_jja_wsm5_obs[0][0], 'title': '(l) WSM5-{0} JJA'.format(dataset)},
'Plot 13': {'data': mbe_son_noto_obs[0][0], 'title': '(m) NoTo-{0} SON'.format(dataset)},
'Plot 14': {'data': mbe_son_wdm7_obs[0][0], 'title': '(n) WDM7-{0} SON'.format(dataset)},
'Plot 15': {'data': mbe_son_wsm7_obs[0][0], 'title': '(o) WSM7-{0} SON'.format(dataset)},
'Plot 16': {'data': mbe_son_wsm5_obs[0][0], 'title': '(p) WSM5-{0} SON'.format(dataset)}}

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
name_out = 'pyplt_maps_bias_{0}_{1}_{2}_RegCM5_{3}.png'.format(var, stats, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
exit()

