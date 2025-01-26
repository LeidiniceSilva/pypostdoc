# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
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

var = 'clh'
dataset = 'ERA5'
domain = 'CSAM-3'
idt, fdt = '2000', '2000'
dt = '{0}-{1}'.format(idt, fdt)

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5_MOIST'

		
def import_obs(param, dataset, season):

	if param == 'cll':
		param_ = 'lcc'
	elif param == 'clm':
		param_ = 'mcc'
	elif param == 'clh':
		param_ = 'hcc'
	else:
		param_ = param

	arq   = '{0}/postproc/obs/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param_][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean


def import_rcm(param, dataset, season):

	arq   = '{0}/postproc/rcm/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	    
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	
		
def configure_subplot(ax):

	ax.set_xticks(np.arange(-82,-34,12), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(-38,-10,6), crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.tick_params(axis='x', labelsize=6, labelcolor='black')
	ax.tick_params(axis='y', labelsize=6, labelcolor='black')
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.coastlines()
	

# Import model and obs dataset
dict_var = {'pr': ['pre', 'precip', 'cmorph', 'precipitation', 'pr'],
'tas': ['tmp', 'tas'],
'tasmax': ['tmx', 'tmax', 'tasmax'],
'tasmin': ['tmn', 'tmin', 'tasmin'],
'evspsblpot': ['pev'],
'rsnl': ['msnlwrf'],
'rsns': ['msnswrf'],
'clt': ['cld', 'clt'],
'cll': ['cll'],
'clm': ['clm'],
'clh': ['clh']}

if var == 'tas' or var == 'tasmax' or var == 'tasmin':
	lat, lon, obs_djf = import_obs(dict_var[var][0], dataset, 'DJF')
	lat, lon, obs_mam = import_obs(dict_var[var][0], dataset, 'MAM')
	lat, lon, obs_jja = import_obs(dict_var[var][0], dataset, 'JJA')
	lat, lon, obs_son = import_obs(dict_var[var][0], dataset, 'SON')

	lat, lon, rcm3_djf = import_rcm(var, 'RegCM5', 'DJF')
	lat, lon, rcm3_mam = import_rcm(var, 'RegCM5', 'MAM')
	lat, lon, rcm3_jja = import_rcm(var, 'RegCM5', 'JJA')
	lat, lon, rcm3_son = import_rcm(var, 'RegCM5', 'SON')

	mbe_djf_rcm3_obs = compute_mbe(rcm3_djf[0], obs_djf)
	mbe_mam_rcm3_obs = compute_mbe(rcm3_mam[0], obs_mam)
	mbe_jja_rcm3_obs = compute_mbe(rcm3_jja[0], obs_jja)
	mbe_son_rcm3_obs = compute_mbe(rcm3_son[0], obs_son)

elif var == 'pr' or var == 'clt':
	lat, lon, obs_djf = import_obs(dict_var[var][0], dataset, 'DJF')
	lat, lon, obs_mam = import_obs(dict_var[var][0], dataset, 'MAM')
	lat, lon, obs_jja = import_obs(dict_var[var][0], dataset, 'JJA')
	lat, lon, obs_son = import_obs(dict_var[var][0], dataset, 'SON')

	lat, lon, rcm3_djf = import_rcm(var, 'RegCM5', 'DJF')
	lat, lon, rcm3_mam = import_rcm(var, 'RegCM5', 'MAM')
	lat, lon, rcm3_jja = import_rcm(var, 'RegCM5', 'JJA')
	lat, lon, rcm3_son = import_rcm(var, 'RegCM5', 'SON')
		
	mbe_djf_rcm3_obs = compute_mbe(rcm3_djf, obs_djf)
	mbe_mam_rcm3_obs = compute_mbe(rcm3_mam, obs_mam)
	mbe_jja_rcm3_obs = compute_mbe(rcm3_jja, obs_jja)
	mbe_son_rcm3_obs = compute_mbe(rcm3_son, obs_son)

else:
	lat, lon, obs_djf = import_obs(dict_var[var][0], dataset, 'DJF')
	lat, lon, obs_mam = import_obs(dict_var[var][0], dataset, 'MAM')
	lat, lon, obs_jja = import_obs(dict_var[var][0], dataset, 'JJA')
	lat, lon, obs_son = import_obs(dict_var[var][0], dataset, 'SON')

	lat, lon, rcm3_djf = import_rcm(var, 'RegCM5', 'DJF')
	lat, lon, rcm3_mam = import_rcm(var, 'RegCM5', 'MAM')
	lat, lon, rcm3_jja = import_rcm(var, 'RegCM5', 'JJA')
	lat, lon, rcm3_son = import_rcm(var, 'RegCM5', 'SON')
		
	mbe_djf_rcm3_obs = compute_mbe(rcm3_djf/100, obs_djf)
	mbe_mam_rcm3_obs = compute_mbe(rcm3_mam/100, obs_mam)
	mbe_jja_rcm3_obs = compute_mbe(rcm3_jja/100, obs_jja)
	mbe_son_rcm3_obs = compute_mbe(rcm3_son/100, obs_son)	

print(mbe_djf_rcm3_obs.shape)

# Plot figure  
fig, axes = plt.subplots(1, 4, figsize=(10, 4), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()
font_size = 6

plot_data = {'Plot 1': {'data': mbe_djf_rcm3_obs, 'title': '(a) RCM3-OBS DJF'},
'Plot 2': {'data': mbe_mam_rcm3_obs, 'title': '(b) RCM3-OBS MAM'},
'Plot 3': {'data': mbe_jja_rcm3_obs, 'title': '(c) RCM3-OBS JJA'},
'Plot 4': {'data': mbe_son_rcm3_obs, 'title': '(d) RCM3-OBS SON'}}

dict_plot = {'pr': ['Precipitation (mm d$^-$$^1$)', np.arange(-10, 11, 1), cm.BrBG],
'tas': ['Air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'tasmax': ['Maximum air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'tasmin': ['Minimum air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'evspsblpot': ['Potential evaporation (mm d$^-$$^1$)', np.arange(-10, 11, 1), cm.seismic],
'rsnl': ['Surface net upward longwave flux (W mm$^-$$^2$)', np.arange(-60, 65, 5), cm.RdBu_r],
'rsns': ['Surface net downward shortwave flux (W mm$^-$$^2$)', np.arange(-60, 65, 5), cm.RdBu_r],
'clt': ['Total cloud cover (%)', np.arange(-90, 100, 10), cm.RdGy],
'cll': ['Low cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'clm': ['Medium cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'clh': ['High cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy]}

for ax, (key, value) in zip(axes, plot_data.items()):
        data = value['data']
        title = value['title']
    
        contour = ax.contourf(lon, lat, data[0], transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
        ax.set_title(title, loc='left', fontsize=font_size, fontweight='bold')
        configure_subplot(ax)

# Set colobar
cbar = fig.colorbar(contour, ax=fig.axes, orientation='horizontal', pad=0.1, aspect=50)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}_{3}.png'.format(var, domain, dataset, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


	
	

	
	
