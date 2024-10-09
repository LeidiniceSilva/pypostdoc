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

var = 'pr'
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
dict_var = {'pr': ['pre', 'precip', 'cmorph', 'precipitation', 'pr'],
'tas': ['tmp', 'tas'],
'tasmax': ['tmx', 'tmax', 'tasmax'],
'tasmin': ['tmn', 'tmin', 'tasmin'],
'evspsblpot': ['pev'],
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

	lat, lon, cmorph_djf = import_obs(dict_var[var][2], 'CSAM-3_CMORPH', 'DJF')
	lat, lon, cmorph_mam = import_obs(dict_var[var][2], 'CSAM-3_CMORPH', 'MAM')
	lat, lon, cmorph_jja = import_obs(dict_var[var][2], 'CSAM-3_CMORPH', 'JJA')
	lat, lon, cmorph_son = import_obs(dict_var[var][2], 'CSAM-3_CMORPH', 'SON')

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
		
	mbe_djf_rcm_cru = compute_mbe(rcm3_djf, cru_djf)
	mbe_mam_rcm_cru = compute_mbe(rcm3_mam, cru_mam)
	mbe_jja_rcm_cru = compute_mbe(rcm3_jja, cru_jja)
	mbe_son_rcm_cru = compute_mbe(rcm3_son, cru_son)	

	mbe_djf_rcm_cpc = compute_mbe(rcm3_djf, cpc_djf)
	mbe_mam_rcm_cpc = compute_mbe(rcm3_mam, cpc_mam)
	mbe_jja_rcm_cpc = compute_mbe(rcm3_jja, cpc_jja)
	mbe_son_rcm_cpc = compute_mbe(rcm3_son, cpc_son)	
	
	mbe_djf_rcm_cmorph = compute_mbe(rcm3_djf, cmorph_djf)
	mbe_mam_rcm_cmorph = compute_mbe(rcm3_mam, cmorph_mam)
	mbe_jja_rcm_cmorph = compute_mbe(rcm3_jja, cmorph_jja)
	mbe_son_rcm_cmorph = compute_mbe(rcm3_son, cmorph_son)

	mbe_djf_rcm_mswep = compute_mbe(rcm3_djf, mswep_djf)
	mbe_mam_rcm_mswep = compute_mbe(rcm3_mam, mswep_mam)
	mbe_jja_rcm_mswep = compute_mbe(rcm3_jja, mswep_jja)
	mbe_son_rcm_mswep = compute_mbe(rcm3_son, mswep_son)
			
	mbe_djf_rcm_era5 = compute_mbe(rcm3_djf, era5_djf)
	mbe_mam_rcm_era5 = compute_mbe(rcm3_mam, era5_mam)
	mbe_jja_rcm_era5 = compute_mbe(rcm3_jja, era5_jja)
	mbe_son_rcm_era5 = compute_mbe(rcm3_son, era5_son)		
	
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

	mbe_djf_rcm_cru = compute_mbe(rcm3_djf[0], cru_djf)
	mbe_mam_rcm_cru = compute_mbe(rcm3_mam[0], cru_mam)
	mbe_jja_rcm_cru = compute_mbe(rcm3_jja[0], cru_jja)
	mbe_son_rcm_cru = compute_mbe(rcm3_son[0], cru_son)
		
	mbe_djf_rcm_era5 = compute_mbe(rcm3_djf[0], era5_djf)
	mbe_mam_rcm_era5 = compute_mbe(rcm3_mam[0], era5_mam)
	mbe_jja_rcm_era5 = compute_mbe(rcm3_jja[0], era5_jja)
	mbe_son_rcm_era5 = compute_mbe(rcm3_son[0], era5_son)	
	
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
	
	mbe_djf_rcm_cru = compute_mbe(rcm3_djf[0], cru_djf)
	mbe_mam_rcm_cru = compute_mbe(rcm3_mam[0], cru_mam)
	mbe_jja_rcm_cru = compute_mbe(rcm3_jja[0], cru_jja)
	mbe_son_rcm_cru = compute_mbe(rcm3_son[0], cru_son)

	mbe_djf_rcm_cpc = compute_mbe(rcm3_djf[0], cpc_djf)
	mbe_mam_rcm_cpc = compute_mbe(rcm3_mam[0], cpc_mam)
	mbe_jja_rcm_cpc = compute_mbe(rcm3_jja[0], cpc_jja)
	mbe_son_rcm_cpc = compute_mbe(rcm3_son[0], cpc_son)
			
	mbe_djf_rcm_era5 = compute_mbe(rcm3_djf[0], era5_djf)
	mbe_mam_rcm_era5 = compute_mbe(rcm3_mam[0], era5_mam)
	mbe_jja_rcm_era5 = compute_mbe(rcm3_jja[0], era5_jja)
	mbe_son_rcm_era5 = compute_mbe(rcm3_son[0], era5_son)	

elif var == 'evspsblpot' or var == 'rsnl' or var == 'rsns' or var == 'cll' or var == 'clm' or var == 'clh':
	lat, lon, era5_djf = import_obs(dict_var[var][0], 'CSAM-3_ERA5', 'DJF')
	lat, lon, era5_mam = import_obs(dict_var[var][0], 'CSAM-3_ERA5', 'MAM')
	lat, lon, era5_jja = import_obs(dict_var[var][0], 'CSAM-3_ERA5', 'JJA')
	lat, lon, era5_son = import_obs(dict_var[var][0], 'CSAM-3_ERA5', 'SON')

	lat, lon, rcm3_djf = import_rcm(var, 'CSAM-3_RegCM5', 'DJF')
	lat, lon, rcm3_mam = import_rcm(var, 'CSAM-3_RegCM5', 'MAM')
	lat, lon, rcm3_jja = import_rcm(var, 'CSAM-3_RegCM5', 'JJA')
	lat, lon, rcm3_son = import_rcm(var, 'CSAM-3_RegCM5', 'SON')
				
	mbe_djf_rcm_era5 = compute_mbe(rcm3_djf, era5_djf)
	mbe_mam_rcm_era5 = compute_mbe(rcm3_mam, era5_mam)
	mbe_jja_rcm_era5 = compute_mbe(rcm3_jja, era5_jja)
	mbe_son_rcm_era5 = compute_mbe(rcm3_son, era5_son)	
	
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
	
	mbe_djf_rcm_cru = compute_mbe(rcm3_djf, cru_djf)
	mbe_mam_rcm_cru = compute_mbe(rcm3_mam, cru_mam)
	mbe_jja_rcm_cru = compute_mbe(rcm3_jja, cru_jja)
	mbe_son_rcm_cru = compute_mbe(rcm3_son, cru_son)
		
	mbe_djf_rcm_era5 = compute_mbe(rcm3_djf, era5_djf)
	mbe_mam_rcm_era5 = compute_mbe(rcm3_mam, era5_mam)
	mbe_jja_rcm_era5 = compute_mbe(rcm3_jja, era5_jja)
	mbe_son_rcm_era5 = compute_mbe(rcm3_son, era5_son)	
	
# Plot figure   
font_size = 6

dict_plot = {'pr': ['Precipitation (mm d$^-$$^1$)', np.arange(-10, 11, 1), cm.BrBG],
'tas': ['Air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'tasmax': ['Maximum air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'tasmin': ['Minimum air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'evspsblpot': ['Potential evaporation (mm d$^-$$^1$)', np.arange(-10, 11, 1), cm.seismic],
'rsnl': ['Surface net upward longwave flux (W mm$^-$$^2$)', np.arange(-60, 65, 5), cm.RdBu_r],
'rsns': ['Surface net downward shortwave flux (W mm$^-$$^2$)', np.arange(-60, 65, 5), cm.RdBu_r],
'clt': ['Total cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'cll': ['Low cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'clm': ['Medium cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'clh': ['High cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy]}

if var == 'pr':
	fig = plt.figure(figsize=(10, 6))

	ax = fig.add_subplot(4, 5, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_rcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(a) CPM3 - CRU', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'DJF', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_rcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(b) CPM3 - CPC', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_rcm_cmorph[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(c) CPM3 - CMORPH', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_rcm_mswep[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(d) CPM3 - MSWEP', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_rcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(e) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_rcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(f) CPM3 - CRU', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'DJF', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_rcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(g) CPM3 - CPC', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_rcm_cmorph[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(h) CPM3 - CMORPH', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 9)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_rcm_mswep[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(i) CPM3 - MSWEP', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 10)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_rcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(j) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 11)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_rcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(k) CPM3 - CRU', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'DJF', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 12)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_rcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(l) CPM3 - CPC', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 13)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_rcm_cmorph[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(m) CPM3 - CMORPH', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 14)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_rcm_mswep[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(n) CPM3 - MSWEP', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 15)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_rcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(o) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 16)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_rcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(p) CPM3 - CRU', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'DJF', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 17)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_rcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(q) CPM3 - CPC', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 18)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_rcm_cmorph[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(r) CPM3 - CMORPH', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 19)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_rcm_mswep[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(s) CPM3 - MSWEP', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 20)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_rcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(t) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.3, 0.015, 0.4]))
	cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)
	
elif var == 'tas':
	fig = plt.figure(figsize=(4, 6))
	
	ax = fig.add_subplot(4, 2, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_rcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(a) CPM3 - CRU', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'DJF', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_rcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(b) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_rcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(c) CPM3 - CRU', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'MAM', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_rcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(d) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_rcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(e) CPM3 - CRU', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'JJA', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_rcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(f) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_rcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(g) CPM3 - CRU', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'SON', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_rcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(h) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.93, 0.3, 0.015, 0.4]))
	cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)
	
elif var == 'tasmax' or var == 'tasmin':
	fig = plt.figure(figsize=(8, 8))

	ax = fig.add_subplot(4, 3, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_rcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(a) CPM3 - CRU', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'DJF', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_rcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(b) CPM3 - CPC', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_rcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(c) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_rcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(d) CPM3 - CRU', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'MAM', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_rcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(e) CPM3 - CPC', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_rcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(f) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_rcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(g) CPM3 - CRU', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'JJA', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_rcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(h) CPM3 - CPC', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 9)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_rcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(i) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 10)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_rcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(j) CPM3 - CRU', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'SON', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 11)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_rcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(k) CPM3 - CPC', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 12)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_rcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(l) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.93, 0.3, 0.015, 0.4]))
	cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)
	
elif var == 'evspsblpot':
	fig = plt.figure(figsize=(4, 8))

	ax = fig.add_subplot(4, 1, 1)  
	map, xx, yy = basemap(lat, lon)
	mbe_djf_mask = maskoceans(xx, yy, mbe_djf_rcm_era5[0])
	plt_map = map.contourf(xx, yy, mbe_djf_mask, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(a) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'DJF', labelpad=20, fontsize=font_size, fontweight='bold')
	
	ax = fig.add_subplot(4, 1, 2)  
	map, xx, yy = basemap(lat, lon)
	mbe_mam_mask = maskoceans(xx, yy, mbe_mam_rcm_era5[0])
	plt_map = map.contourf(xx, yy, mbe_mam_mask, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(b) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'MAM', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 1, 3)  
	map, xx, yy = basemap(lat, lon)
	mbe_jja_mask = maskoceans(xx, yy, mbe_jja_rcm_era5[0])
	plt_map = map.contourf(xx, yy, mbe_jja_mask, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(c) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'JJA', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 1, 4)  
	map, xx, yy = basemap(lat, lon)
	mbe_son_mask = maskoceans(xx, yy, mbe_son_rcm_era5[0])
	plt_map = map.contourf(xx, yy, mbe_son_mask, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
	plt.title(u'(d) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'SON', labelpad=20, fontsize=font_size, fontweight='bold')

	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.85, 0.3, 0.015, 0.4]))
	cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)

elif var == 'rsnl' or var == 'rsns':
	fig = plt.figure(figsize=(4, 8))

	ax = fig.add_subplot(4, 1, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_rcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(a) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'DJF', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 1, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_rcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(b) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'MAM', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 1, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_rcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(c) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'JJA', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 1, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_rcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(d) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'SON', labelpad=20, fontsize=font_size, fontweight='bold')

	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.85, 0.3, 0.015, 0.4]))
	cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)

elif var == 'clt':
	fig = plt.figure(figsize=(4, 6))

	ax = fig.add_subplot(4, 2, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_rcm_cru[0]/100, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(a) CPM3 - CRU', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'DJF', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_rcm_era5[0]/100, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(b) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_rcm_cru[0]/100, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(c) CPM3 - CRU', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'MAM', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_rcm_era5[0]/100, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(d) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_rcm_cru[0]/100, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(e) CPM3 - CRU', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'JJA', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_rcm_era5[0]/100, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(f) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_rcm_cru[0]/100, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(g) CPM3 - CRU', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'SON', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_rcm_era5[0]/100, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(h) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.93, 0.3, 0.015, 0.4]))
	cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)

else:
	fig = plt.figure(figsize=(4, 8))

	ax = fig.add_subplot(4, 1, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_rcm_era5[0]/100, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(a) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'DJF', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 1, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_rcm_era5[0]/100, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(b) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'MAM', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 1, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_rcm_era5[0]/100, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(c) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'JJA', labelpad=20, fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 1, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_rcm_era5[0]/100, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(d) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
	plt.ylabel(u'SON', labelpad=20, fontsize=font_size, fontweight='bold')

	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.85, 0.3, 0.015, 0.4]))
	cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)
		
# Path out to save figure
path_out = '{0}/user/mdasilva/CORDEX/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
