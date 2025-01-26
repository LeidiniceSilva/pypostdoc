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

    ax.set_xticks(np.arange(-40,65,25), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(30,85,15), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.tick_params(axis='x', labelsize=6, labelcolor='black')
    ax.tick_params(axis='y', labelsize=6, labelcolor='black')
    ax.grid(c='k', ls='--', alpha=0.4)
    ax.coastlines()
	
	
# Import model and obs dataset
dict_var = {'CAPE': ['cape'],
'CIN': ['cin'],
'evspsblpot': ['mper'],
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

dict_plot = {'CAPE': ['Convective Available Potential Energy (J kg$^-$$^1$)', np.arange(0, 2550, 50), custom_cmap],
'CIN': ['Convective inhibition (J kg$^-$$^1$)', np.arange(0, 675, 25), custom_cmap],
'LI': ['Lifted Index (Kelvin)', np.arange(-100, 110, 10), cm.jet],
'evspsblpot': ['Potential evaporation (mm d$^-$$^1$)', np.arange(0, 8.5, 0.5), cm.jet]}
	
if var == 'evspsblpot':
	fig = plt.figure(figsize=(4, 6))

	ax = fig.add_subplot(4, 2, 1)  
	map, xx, yy = basemap(lat, lon)
	era5_djf_mask = maskoceans(xx, yy, era5_djf[0])	
	plt_map = map.contourf(xx, yy, era5_djf_mask, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(a) ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'DJF', labelpad=20, fontsize=font_size, fontweight='bold')
	
	ax = fig.add_subplot(4, 2, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, rcm3_djf[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(b) CPM3', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 3)  
	map, xx, yy = basemap(lat, lon)
	era5_mam_mask = maskoceans(xx, yy, era5_mam[0])
	plt_map = map.contourf(xx, yy, era5_mam_mask, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(c) ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'MAM', labelpad=20, fontsize=font_size, fontweight='bold')
	
	ax = fig.add_subplot(4, 2, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, rcm3_mam[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(d) CPM3', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 5)  
	map, xx, yy = basemap(lat, lon)
	era5_jja_mask = maskoceans(xx, yy, era5_jja[0])
	plt_map = map.contourf(xx, yy, era5_jja_mask, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(e) ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'JJA', labelpad=20, fontsize=font_size, fontweight='bold')
	
	ax = fig.add_subplot(4, 2, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, rcm3_jja[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(f) CPM3', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 7)  
	map, xx, yy = basemap(lat, lon)
	era5_son_mask = maskoceans(xx, yy, era5_son[0])
	plt_map = map.contourf(xx, yy, era5_son_mask, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(g) ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'SON', labelpad=20, fontsize=font_size, fontweight='bold')
	
	ax = fig.add_subplot(4, 2, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, rcm3_son[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(h) CPM3', loc='left', fontsize=font_size, fontweight='bold')

	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.93, 0.3, 0.015, 0.4]))
	cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)
	
elif var == 'LI':
	fig = plt.figure(figsize=(4, 6))

	ax = fig.add_subplot(4, 2, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_djf[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(a) NCEP reanalysis 1', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'DJF', labelpad=20, fontsize=font_size, fontweight='bold')
	
	ax = fig.add_subplot(4, 2, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, rcm3_djf[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(b) CPM3', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_mam[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(c) NCEP reanalysis 1', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'MAM', labelpad=20, fontsize=font_size, fontweight='bold')
	
	ax = fig.add_subplot(4, 2, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, rcm3_mam[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(d) CPM3', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_jja[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(e) NCEP reanalysis 1', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'JJA', labelpad=20, fontsize=font_size, fontweight='bold')
	
	ax = fig.add_subplot(4, 2, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, rcm3_jja[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(f) CPM3', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_son[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(g) NCEP reanalysis 1', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'SON', labelpad=20, fontsize=font_size, fontweight='bold')
	
	ax = fig.add_subplot(4, 2, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, rcm3_son[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(h) CPM3', loc='left', fontsize=font_size, fontweight='bold')

	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.93, 0.3, 0.015, 0.4]))
	cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)

else:
	fig = plt.figure(figsize=(4, 6))

	ax = fig.add_subplot(4, 2, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_djf[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(a) ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'DJF', labelpad=20, fontsize=font_size, fontweight='bold')
	
	ax = fig.add_subplot(4, 2, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, rcm3_djf[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(b) CPM3', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_mam[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(c) ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'MAM', labelpad=20, fontsize=font_size, fontweight='bold')
	
	ax = fig.add_subplot(4, 2, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, rcm3_mam[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(d) CPM3', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_jja[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(e) ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'JJA', labelpad=20, fontsize=font_size, fontweight='bold')
	
	ax = fig.add_subplot(4, 2, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, rcm3_jja[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(f) CPM3', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_son[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(g) ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'SON', labelpad=20, fontsize=font_size, fontweight='bold')
	
	ax = fig.add_subplot(4, 2, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, rcm3_son[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(h) CPM3', loc='left', fontsize=font_size, fontweight='bold')

	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.93, 0.3, 0.015, 0.4]))
	cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)
			
# Path out to save figure
path_out = '{0}/user/mdasilva/CORDEX/figs'.format(path)
name_out = 'pyplt_maps_clim_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
