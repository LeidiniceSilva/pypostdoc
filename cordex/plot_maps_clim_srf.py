# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot bias maps"

import os
import netCDF4
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import maskoceans
from import_climate_tools import compute_mbe

var = 'evspsblpot'
domain = 'CSAM-3'
idt, fdt = '2000', '2009'
dt = '{0}-{1}'.format(idt, fdt)

path = '/marconi/home/userexternal/mdasilva'

		
def import_obs(param, dataset, season):

	arq   = '{0}/user/mdasilva/CORDEX/post_evaluate/obs/{1}_{2}_{3}_2000-2009_lonlat.nc'.format(path, param, dataset, season)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean


def import_rcm(param, dataset, season):

	arq   = '{0}/user/mdasilva/CORDEX/post_evaluate/rcm/{1}_{2}_{3}_2000-2005_lonlat.nc'.format(path, param, dataset, season)    
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	
		
def basemap(lat, lon):
	
	map = Basemap(projection='cyl', llcrnrlon=-80., llcrnrlat=-38., urcrnrlon=-34.,urcrnrlat=-10., resolution='c')
	map.drawmeridians(np.arange(-80., -34., 12.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-38., -8., 6.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black')
	
	lons, lats = np.meshgrid(lon, lat)
	xx, yy = map(lons,lats)
	
	return map, xx, yy
	
	
# Import model and obs dataset
dict_var = {'pr': ['pre', 'precip', 'hrf', 'precipitation', 'pr'],
'tas': ['tmp', 'tas'],
'tasmax': ['tmx', 'tmax', 'tasmax'],
'tasmin': ['tmn', 'tmin', 'tasmin'],
'evspsblpot': ['mper'],
'LI': ['lftx'],
'rsnl': ['msnlwrf'],
'rsns': ['msnswrf'],
'clt': ['cld', 'clt'],
'cll': ['lcc'],
'clm': ['mcc'],
'clh': ['hcc']}

if var == 'pr':	
	lat, lon, cru_djf = import_obs(dict_var[var][0], 'CSAM-3_CRU', 'DJF')
	lat, lon, cru_mam = import_obs(dict_var[var][0], 'CSAM-3_CRU', 'MAM')
	lat, lon, cru_jja = import_obs(dict_var[var][0], 'CSAM-3_CRU', 'JJA')
	lat, lon, cru_son = import_obs(dict_var[var][0], 'CSAM-3_CRU', 'SON')

	lat, lon, cpc_djf = import_obs(dict_var[var][1], 'CSAM-3_CPC', 'DJF')
	lat, lon, cpc_mam = import_obs(dict_var[var][1], 'CSAM-3_CPC', 'MAM')
	lat, lon, cpc_jja = import_obs(dict_var[var][1], 'CSAM-3_CPC', 'JJA')
	lat, lon, cpc_son = import_obs(dict_var[var][1], 'CSAM-3_CPC', 'SON')

	lat, lon, trmm_djf = import_obs(dict_var[var][2], 'CSAM-3_TRMM', 'DJF')
	lat, lon, trmm_mam = import_obs(dict_var[var][2], 'CSAM-3_TRMM', 'MAM')
	lat, lon, trmm_jja = import_obs(dict_var[var][2], 'CSAM-3_TRMM', 'JJA')
	lat, lon, trmm_son = import_obs(dict_var[var][2], 'CSAM-3_TRMM', 'SON')

	lat, lon, mswep_djf = import_obs(dict_var[var][3], 'CSAM-3_MSWEP', 'DJF')
	lat, lon, mswep_mam = import_obs(dict_var[var][3], 'CSAM-3_MSWEP', 'MAM')
	lat, lon, mswep_jja = import_obs(dict_var[var][3], 'CSAM-3_MSWEP', 'JJA')
	lat, lon, mswep_son = import_obs(dict_var[var][3], 'CSAM-3_MSWEP', 'SON')
	
	lat, lon, era5_djf = import_obs(dict_var[var][4], 'CSAM-3_ERA5', 'DJF')
	lat, lon, era5_mam = import_obs(dict_var[var][4], 'CSAM-3_ERA5', 'MAM')
	lat, lon, era5_jja = import_obs(dict_var[var][4], 'CSAM-3_ERA5', 'JJA')
	lat, lon, era5_son = import_obs(dict_var[var][4], 'CSAM-3_ERA5', 'SON')

	lat, lon, rcm3_djf = import_rcm(var, 'CSAM-3_RegCM5', 'DJF')
	lat, lon, rcm3_mam = import_rcm(var, 'CSAM-3_RegCM5', 'MAM')
	lat, lon, rcm3_jja = import_rcm(var, 'CSAM-3_RegCM5', 'JJA')
	lat, lon, rcm3_son = import_rcm(var, 'CSAM-3_RegCM5', 'SON')		
	
elif var == 'tas':
	lat, lon, cru_djf = import_obs(dict_var[var][0], 'CSAM-3_CRU', 'DJF')
	lat, lon, cru_mam = import_obs(dict_var[var][0], 'CSAM-3_CRU', 'MAM')
	lat, lon, cru_jja = import_obs(dict_var[var][0], 'CSAM-3_CRU', 'JJA')
	lat, lon, cru_son = import_obs(dict_var[var][0], 'CSAM-3_CRU', 'SON')

	lat, lon, era5_djf = import_obs(dict_var[var][1], 'CSAM-3_ERA5', 'DJF')
	lat, lon, era5_mam = import_obs(dict_var[var][1], 'CSAM-3_ERA5', 'MAM')
	lat, lon, era5_jja = import_obs(dict_var[var][1], 'CSAM-3_ERA5', 'JJA')
	lat, lon, era5_son = import_obs(dict_var[var][1], 'CSAM-3_ERA5', 'SON')

	lat, lon, rcm3_djf = import_rcm(var, 'CSAM-3_RegCM5', 'DJF')
	lat, lon, rcm3_mam = import_rcm(var, 'CSAM-3_RegCM5', 'MAM')
	lat, lon, rcm3_jja = import_rcm(var, 'CSAM-3_RegCM5', 'JJA')
	lat, lon, rcm3_son = import_rcm(var, 'CSAM-3_RegCM5', 'SON')	
	
elif var == 'tasmax' or var == 'tasmin':
	lat, lon, cru_djf = import_obs(dict_var[var][0], 'CSAM-3_CRU', 'DJF')
	lat, lon, cru_mam = import_obs(dict_var[var][0], 'CSAM-3_CRU', 'MAM')
	lat, lon, cru_jja = import_obs(dict_var[var][0], 'CSAM-3_CRU', 'JJA')
	lat, lon, cru_son = import_obs(dict_var[var][0], 'CSAM-3_CRU', 'SON')

	lat, lon, cpc_djf = import_obs(dict_var[var][1], 'CSAM-3_CPC', 'DJF')
	lat, lon, cpc_mam = import_obs(dict_var[var][1], 'CSAM-3_CPC', 'MAM')
	lat, lon, cpc_jja = import_obs(dict_var[var][1], 'CSAM-3_CPC', 'JJA')
	lat, lon, cpc_son = import_obs(dict_var[var][1], 'CSAM-3_CPC', 'SON')

	lat, lon, era5_djf = import_obs(dict_var[var][2], 'CSAM-3_ERA5', 'DJF')
	lat, lon, era5_mam = import_obs(dict_var[var][2], 'CSAM-3_ERA5', 'MAM')
	lat, lon, era5_jja = import_obs(dict_var[var][2], 'CSAM-3_ERA5', 'JJA')
	lat, lon, era5_son = import_obs(dict_var[var][2], 'CSAM-3_ERA5', 'SON')

	lat, lon, rcm3_djf = import_rcm(var, 'CSAM-3_RegCM5', 'DJF')
	lat, lon, rcm3_mam = import_rcm(var, 'CSAM-3_RegCM5', 'MAM')
	lat, lon, rcm3_jja = import_rcm(var, 'CSAM-3_RegCM5', 'JJA')
	lat, lon, rcm3_son = import_rcm(var, 'CSAM-3_RegCM5', 'SON')

elif var == 'LI' or var == 'evspsblpot' or var == 'rsnl' or var == 'rsns' or var == 'cll' or var == 'clm' or var == 'clh':
	lat, lon, era5_djf = import_obs(dict_var[var][0], 'CSAM-3_ERA5', 'DJF')
	lat, lon, era5_mam = import_obs(dict_var[var][0], 'CSAM-3_ERA5', 'MAM')
	lat, lon, era5_jja = import_obs(dict_var[var][0], 'CSAM-3_ERA5', 'JJA')
	lat, lon, era5_son = import_obs(dict_var[var][0], 'CSAM-3_ERA5', 'SON')

	lat, lon, rcm3_djf = import_rcm(var, 'CSAM-3_RegCM5', 'DJF')
	lat, lon, rcm3_mam = import_rcm(var, 'CSAM-3_RegCM5', 'MAM')
	lat, lon, rcm3_jja = import_rcm(var, 'CSAM-3_RegCM5', 'JJA')
	lat, lon, rcm3_son = import_rcm(var, 'CSAM-3_RegCM5', 'SON')
	
else:
	lat, lon, cru_djf = import_obs(dict_var[var][0], 'CSAM-3_CRU', 'DJF')
	lat, lon, cru_mam = import_obs(dict_var[var][0], 'CSAM-3_CRU', 'MAM')
	lat, lon, cru_jja = import_obs(dict_var[var][0], 'CSAM-3_CRU', 'JJA')
	lat, lon, cru_son = import_obs(dict_var[var][0], 'CSAM-3_CRU', 'SON')

	lat, lon, era5_djf = import_obs(dict_var[var][1], 'CSAM-3_ERA5', 'DJF')
	lat, lon, era5_mam = import_obs(dict_var[var][1], 'CSAM-3_ERA5', 'MAM')
	lat, lon, era5_jja = import_obs(dict_var[var][1], 'CSAM-3_ERA5', 'JJA')
	lat, lon, era5_son = import_obs(dict_var[var][1], 'CSAM-3_ERA5', 'SON')

	lat, lon, rcm3_djf = import_rcm(var, 'CSAM-3_RegCM5', 'DJF')
	lat, lon, rcm3_mam = import_rcm(var, 'CSAM-3_RegCM5', 'MAM')
	lat, lon, rcm3_jja = import_rcm(var, 'CSAM-3_RegCM5', 'JJA')
	lat, lon, rcm3_son = import_rcm(var, 'CSAM-3_RegCM5', 'SON')	

# Plot figure   
font_size = 6

dict_plot = {'pr': ['Precipitation (mm d$^-$$^1$)', np.arange(0, 17, 1), cm.Blues],
'tas': ['Air temperature (°C)', np.arange(10, 39, 2), cm.Reds],
'tasmax': ['Maximum air temperature (°C)', np.arange(10, 39, 2), cm.Reds],
'tasmin': ['Minimum air temperature (°C)', np.arange(10, 39, 2), cm.Reds],
'LI': ['Lifted Index (Kelvin)', np.arange(-100, 100, 10), cm.jet],
'evspsblpot': ['Potential evaporation (mm d$^-$$^1$)', np.arange(0, 8.5, 0.5), cm.jet],
'rsnl': ['Surface net upward longwave flux (W mm$^-$$^2$)', np.arange(0, 270, 10), cm.jet],
'rsns': ['Surface net downward shortwave flux (W mm$^-$$^2$)', np.arange(0, 270, 10), cm.jet],
'clt': ['Total cloud cover (0-1)', np.arange(0, 1, 0.1), cm.Greys],
'cll': ['Low cloud cover (0-1)', np.arange(0, 1, 0.1), cm.Greys],
'clm': ['Medium cloud cover (0-1)', np.arange(0, 1, 0.1), cm.Greys],
'clh': ['High cloud cover (0-1)', np.arange(0, 1, 0.1), cm.Greys]}

	
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
	
else:
	fig = plt.figure(figsize=(4, 6))

	ax = fig.add_subplot(4, 2, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_djf[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(a) NCEP reanalysis', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'DJF', labelpad=20, fontsize=font_size, fontweight='bold')
	
	ax = fig.add_subplot(4, 2, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, rcm3_djf[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(b) CPM3', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_mam[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(c) NCEP reanalysis', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'MAM', labelpad=20, fontsize=font_size, fontweight='bold')
	
	ax = fig.add_subplot(4, 2, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, rcm3_mam[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(d) CPM3', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_jja[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(e) NCEP reanalysis', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'JJA', labelpad=20, fontsize=font_size, fontweight='bold')
	
	ax = fig.add_subplot(4, 2, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, rcm3_jja[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(f) CPM3', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, era5_son[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(g) NCEP reanalysis', loc='left', fontsize=font_size, fontweight='bold')
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
