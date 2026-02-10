# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Feb 02, 2026"
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

parser = argparse.ArgumentParser(description='Plot trend maps')
parser.add_argument('--var', required=True, help='Variable name')
parser.add_argument('--rcm_ii', required=True, help='GCM driven')
parser.add_argument('--gcm_iii', required=True, help='GCM name')
args = parser.parse_args()

# Configuration
var = args.var
rcm_i = 'RegCM5-ERA5_ICTP' 
rcm_ii = args.rcm_ii # RegCM5-ECEarth_ICTP RegCM5-MPI_ICTP RegCM5-Nor_USP
gcm_iii = args.gcm_iii # EC-Earth3-Veg MPI-ESM1-2-HR NorESM2-MM

obs = 'ERA5'
dt = '1971-1972'
domain = 'SAM-12'
latlon = [-105, -16, -57, 18]

dict_var = {'pr': ['tp'],
'tas': ['t2m']}

font_size = 6
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/{0}'.format(domain)


def import_obs(param, dataset, season):

	arq   = '{0}/postproc/obs/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

	return lat, lon, mean


def import_rcm_i(param, dataset, season):

	arq   = '{0}/postproc/rcm/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

	return lat, lon, mean


def import_rcm_ii(param, dataset, season):

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


def import_gcm(param, dataset, season):

	arq   = '{0}/postproc/gcm/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
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

lat, lon, rcm_i_djf = import_rcm_i(var, rcm_i, 'DJF')
lat, lon, rcm_i_mam = import_rcm_i(var, rcm_i, 'MAM')
lat, lon, rcm_i_jja = import_rcm_i(var, rcm_i, 'JJA')
lat, lon, rcm_i_son = import_rcm_i(var, rcm_i, 'SON')

lat, lon, rcm_ii_djf = import_rcm_ii(var, rcm_ii, 'DJF')
lat, lon, rcm_ii_mam = import_rcm_ii(var, rcm_ii, 'MAM')
lat, lon, rcm_ii_jja = import_rcm_ii(var, rcm_ii, 'JJA')
lat, lon, rcm_ii_son = import_rcm_ii(var, rcm_ii, 'SON')

lat, lon, gcm_iii_djf = import_gcm(var, gcm_iii, 'DJF')
lat, lon, gcm_iii_mam = import_gcm(var, gcm_iii, 'MAM')
lat, lon, gcm_iii_jja = import_gcm(var, gcm_iii, 'JJA')
lat, lon, gcm_iii_son = import_gcm(var, gcm_iii, 'SON')

mbe_djf_rcm_i_obs = compute_mbe(rcm_i_djf, obs_djf)
mbe_mam_rcm_i_obs = compute_mbe(rcm_i_mam, obs_mam)
mbe_jja_rcm_i_obs = compute_mbe(rcm_i_jja, obs_jja)
mbe_son_rcm_i_obs = compute_mbe(rcm_i_son, obs_son)

mbe_djf_rcm_ii_obs = compute_mbe(rcm_ii_djf, obs_djf)
mbe_mam_rcm_ii_obs = compute_mbe(rcm_ii_mam, obs_mam)
mbe_jja_rcm_ii_obs = compute_mbe(rcm_ii_jja, obs_jja)
mbe_son_rcm_ii_obs = compute_mbe(rcm_ii_son, obs_son)

mbe_djf_gcm_iii_obs = compute_mbe(gcm_iii_djf, obs_djf)
mbe_mam_gcm_iii_obs = compute_mbe(gcm_iii_mam, obs_mam)
mbe_jja_gcm_iii_obs = compute_mbe(gcm_iii_jja, obs_jja)
mbe_son_gcm_iii_obs = compute_mbe(gcm_iii_son, obs_son)

# Plot figure
fig, axes = plt.subplots(3, 4, figsize=(12, 8), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()

dict_plot = {'pr': ['Bias of precipitation (mm d$^-$$^1$)', np.arange(-8, 8.5, 0.5), cm.BrBG],
'tas': ['Bias of air temperature (°C)', np.arange(-8, 8.5, 0.5), cm.RdBu_r]}

plot_data = {'Plot 1': {'data': mbe_djf_rcm_i_obs, 'title': '(a) {0}-{1} DJF'.format(rcm_i, obs)},
'Plot 2': {'data': mbe_mam_rcm_i_obs, 'title': '(b) {0}-{1} MAM'.format(rcm_i, obs)},
'Plot 3': {'data': mbe_jja_rcm_i_obs, 'title': '(c) {0}-{1} JJA'.format(rcm_i, obs)},
'Plot 4': {'data': mbe_son_rcm_i_obs, 'title': '(d) {0}-{1} SON'.format(rcm_i, obs)},
'Plot 5': {'data': mbe_djf_rcm_ii_obs, 'title': '(e) {0}-{1} DJF'.format(rcm_ii, obs)},
'Plot 6': {'data': mbe_mam_rcm_ii_obs, 'title': '(f) {0}-{1} MAM'.format(rcm_ii, obs)},
'Plot 7': {'data': mbe_jja_rcm_ii_obs, 'title': '(g) {0}-{1} JJA'.format(rcm_ii, obs)},
'Plot 8': {'data': mbe_son_rcm_ii_obs, 'title': '(h) {0}-{1} SON'.format(rcm_ii, obs)},
'Plot 9': {'data': mbe_djf_gcm_iii_obs, 'title': '(e) {0}-{1} DJF'.format(gcm_iii, obs)},
'Plot 10': {'data': mbe_mam_gcm_iii_obs, 'title': '(f) {0}-{1} MAM'.format(gcm_iii, obs)},
'Plot 11': {'data': mbe_jja_gcm_iii_obs, 'title': '(g) {0}-{1} JJA'.format(gcm_iii, obs)},
'Plot 12': {'data': mbe_son_gcm_iii_obs, 'title': '(h) {0}-{1} SON'.format(gcm_iii, obs)}}

for ax, (key, value) in zip(axes, plot_data.items()):
	data = value['data']
	title = value['title']
    
	contour = ax.contourf(lon, lat, data, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both')
	ax.set_title(title, loc='left', fontsize=font_size, fontweight='bold')
	configure_subplot(ax)

# Set colobar
cbar = fig.colorbar(contour, ax=fig.axes, orientation='vertical', pad=0.025, aspect=50)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_{2}_{3}.png'.format(var, domain, rcm_ii, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
#plt.show()
exit()

