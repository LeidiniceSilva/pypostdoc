# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
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

var = 'pr'
domain = 'CSAM-3'
idt, fdt = '2000', '2000'
dt = '{0}-{1}'.format(idt, fdt)

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5'

		
def import_obs(param, dataset, season):

	arq   = '{0}/postproc/obs/{1}_{2}_{3}_{4}_lonlat.nc'.format(path, param, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean


def import_rcm(param, dataset, season):

	arq   = '{0}/postproc/rcm/{1}_{2}_{3}_{4}_lonlat.nc'.format(path, param, dataset, season, dt)    
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
dict_var = {'pr': ['pr'],
'evspsblpot': ['mper'],
'CAPE': ['cape'],
'CIN': ['cin'],
'LI': ['lftx']}

lat, lon, era5_djf = import_obs(dict_var[var][0], 'CSAM-3_ERA5', 'DJF')
lat, lon, era5_mam = import_obs(dict_var[var][0], 'CSAM-3_ERA5', 'MAM')
lat, lon, era5_jja = import_obs(dict_var[var][0], 'CSAM-3_ERA5', 'JJA')
lat, lon, era5_son = import_obs(dict_var[var][0], 'CSAM-3_ERA5', 'SON')

lat, lon, rcm3_djf = import_rcm(var, 'CSAM-3_RegCM5', 'DJF')
lat, lon, rcm3_mam = import_rcm(var, 'CSAM-3_RegCM5', 'MAM')
lat, lon, rcm3_jja = import_rcm(var, 'CSAM-3_RegCM5', 'JJA')
lat, lon, rcm3_son = import_rcm(var, 'CSAM-3_RegCM5', 'SON')
	
# Plot figure   
fig, axes = plt.subplots(4, 2, figsize=(4, 6), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()
font_size = 6

# Step 1: Define a custom colormap similar to the image colors
colors = [
    (1, 1, 1),  # White for the lowest values
    (0.8, 0.9, 1),  # Light blue
    (0.6, 0.7, 1),  # Medium blue
    (0.4, 0.4, 1),  # Dark blue
    (1, 0.8, 0.6),  # Light orange
    (1, 0.5, 0.2),  # Orange
    (1, 0.2, 0),    # Red
    (0.6, 0, 0.6),  # Dark red/purple
    (0.4, 0, 0.4)   # Dark purple
]

# Create a colormap object
custom_cmap = LinearSegmentedColormap.from_list('cape_colormap', colors, N=256)

dict_plot = {'pr': ['Precipitation (mm d$^-$$^1$)', np.arange(0, 17, 1), cm.BrBG],
'evspsblpot': ['Potential evaporation (mm d$^-$$^1$)', np.arange(0, 8.5, 0.5), cm.jet],
'CAPE': ['Convective Available Potential Energy (J kg$^-$$^1$)', np.arange(0, 2550, 50), custom_cmap],
'CIN': ['Convective inhibition (J kg$^-$$^1$)', np.arange(0, 675, 25), custom_cmap],
'LI': ['Lifted Index (Kelvin)', np.arange(-100, 110, 10), cm.jet]}

plot_data = {'Plot 1': {'data': era5_djf, 'title': '(a) OBS DJF'},
'Plot 2': {'data': rcm3_djf, 'title': '(b) RCM3 DJF'},
'Plot 3': {'data': era5_mam, 'title': '(c) OBS JJA'},
'Plot 4': {'data': rcm3_mam, 'title': '(d) RCM3 SON'},
'Plot 5': {'data': era5_jja, 'title': '(e) OBS DJF'},
'Plot 6': {'data': rcm3_jja, 'title': '(f) RCM3 MAM'},
'Plot 7': {'data': era5_son, 'title': '(g) OBS JJA'},
'Plot 8': {'data': rcm3_son, 'title': '(h) RCM3 SON'}}

if var == 'evspsblpot':
	for ax, (key, value) in zip(axes, plot_data.items()):
        	data = value['data']
        	title = value['title']
    
        	contour = ax.contourf(lon, lat, data[0], transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
        	ax.set_title(title, loc='left', fontsize=font_size, fontweight='bold')
        	configure_subplot(ax)
elif var == 'LI':
	for ax, (key, value) in zip(axes, plot_data.items()):
        	data = value['data']
        	title = value['title']
    
        	contour = ax.contourf(lon, lat, data[0], transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
        	ax.set_title(title, loc='left', fontsize=font_size, fontweight='bold')
        	configure_subplot(ax)
else:
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
name_out = 'pyplt_maps_clim_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
