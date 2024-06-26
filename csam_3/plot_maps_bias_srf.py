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
from import_climate_tools import compute_mbe

var = 'tasmin'
domain = 'CSAM-3'
idt, fdt = '2000', '2000'
dt = '{0}-{1}'.format(idt, fdt)

path = '/marconi/home/userexternal/mdasilva'

			
def import_obs(param, domain, dataset, season):

	arq   = '{0}/user/mdasilva/CORDEX/post_evaluate/obs/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean


def import_rcm(param, domain, dataset, season):

	arq   = '{0}/user/mdasilva/CORDEX/post_evaluate/rcm/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)    
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
dict_var = {
'pr': ['pre', 'pre', 'precip', 'sat_gauge_precip', 'tp'],
'tas': ['tmp', 'tmp', 't2m'],
'tasmax': ['tmx', 'tmax', 'mx2t'],
'tasmin': ['tmn', 'tmin', 'mn2t'],
'clt': ['cld', 'tcc'],
'rsnl': ['msnlwrf']
}

if var == 'pr':	
	lat, lon, cru_djf = import_obs(dict_var[var][1], domain, 'CRU', 'DJF')
	lat, lon, cru_mam = import_obs(dict_var[var][1], domain, 'CRU', 'MAM')
	lat, lon, cru_jja = import_obs(dict_var[var][1], domain, 'CRU', 'JJA')
	lat, lon, cru_son = import_obs(dict_var[var][1], domain, 'CRU', 'SON')

	lat, lon, cpc_djf = import_obs(dict_var[var][2], domain, 'CPC', 'DJF')
	lat, lon, cpc_mam = import_obs(dict_var[var][2], domain, 'CPC', 'MAM')
	lat, lon, cpc_jja = import_obs(dict_var[var][2], domain, 'CPC', 'JJA')
	lat, lon, cpc_son = import_obs(dict_var[var][2], domain, 'CPC', 'SON')

	lat, lon, gpcp_djf = import_obs(dict_var[var][3], domain, 'GPCP', 'DJF')
	lat, lon, gpcp_mam = import_obs(dict_var[var][3], domain, 'GPCP', 'MAM')
	lat, lon, gpcp_jja = import_obs(dict_var[var][3], domain, 'GPCP', 'JJA')
	lat, lon, gpcp_son = import_obs(dict_var[var][3], domain, 'GPCP', 'SON')

	lat, lon, era5_djf = import_obs(dict_var[var][4], var, domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_obs(dict_var[var][4], var, domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_obs(dict_var[var][4], var, domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_obs(dict_var[var][4], var, domain, 'ERA5', 'SON')

	lat, lon, regcm_djf = import_rcm(var, domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_rcm(var, domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_rcm(var, domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_rcm(var, domain, 'RegCM5', 'SON')
	
	mbe_djf_regcm_cru = compute_mbe(regcm_djf, cru_djf)
	mbe_mam_regcm_cru = compute_mbe(regcm_mam, cru_mam)
	mbe_jja_regcm_cru = compute_mbe(regcm_jja, cru_jja)
	mbe_son_regcm_cru = compute_mbe(regcm_son, cru_son)	

	mbe_djf_regcm_cpc = compute_mbe(regcm_djf, cpc_djf)
	mbe_mam_regcm_cpc = compute_mbe(regcm_mam, cpc_mam)
	mbe_jja_regcm_cpc = compute_mbe(regcm_jja, cpc_jja)
	mbe_son_regcm_cpc = compute_mbe(regcm_son, cpc_son)	
	
	mbe_djf_regcm_gpcp = compute_mbe(regcm_djf, gpcp_djf)
	mbe_mam_regcm_gpcp = compute_mbe(regcm_mam, gpcp_mam)
	mbe_jja_regcm_gpcp = compute_mbe(regcm_jja, gpcp_jja)
	mbe_son_regcm_gpcp = compute_mbe(regcm_son, gpcp_son)
		
	mbe_djf_regcm_era5 = compute_mbe(regcm_djf, era5_djf)
	mbe_mam_regcm_era5 = compute_mbe(regcm_mam, era5_mam)
	mbe_jja_regcm_era5 = compute_mbe(regcm_jja, era5_jja)
	mbe_son_regcm_era5 = compute_mbe(regcm_son, era5_son)		

elif var == 'tas':
	lat, lon, cru_djf = import_obs(dict_var[var][1], dict_var[var][1], domain, 'CRU', 'DJF')
	lat, lon, cru_mam = import_obs(dict_var[var][1], dict_var[var][1], domain, 'CRU', 'MAM')
	lat, lon, cru_jja = import_obs(dict_var[var][1], dict_var[var][1], domain, 'CRU', 'JJA')
	lat, lon, cru_son = import_obs(dict_var[var][1], dict_var[var][1], domain, 'CRU', 'SON')

	lat, lon, era5_djf = import_obs(dict_var[var][2], var, domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_obs(dict_var[var][2], var, domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_obs(dict_var[var][2], var, domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_obs(dict_var[var][2], var, domain, 'ERA5', 'SON')

	lat, lon, regcm_djf = import_rcm(var, domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_rcm(var, domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_rcm(var, domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_rcm(var, domain, 'RegCM5', 'SON')
	
	mbe_djf_regcm_cru = compute_mbe(regcm_djf[0], cru_djf)
	mbe_mam_regcm_cru = compute_mbe(regcm_mam[0], cru_mam)
	mbe_jja_regcm_cru = compute_mbe(regcm_jja[0], cru_jja)
	mbe_son_regcm_cru = compute_mbe(regcm_son[0], cru_son)	
		
	mbe_djf_regcm_era5 = compute_mbe(regcm_djf[0], era5_djf)
	mbe_mam_regcm_era5 = compute_mbe(regcm_mam[0], era5_mam)
	mbe_jja_regcm_era5 = compute_mbe(regcm_jja[0], era5_jja)
	mbe_son_regcm_era5 = compute_mbe(regcm_son[0], era5_son)
	
elif var == 'tasmax':
	lat, lon, cru_djf = import_obs(dict_var[var][0], domain, 'CRU', 'DJF')
	lat, lon, cru_mam = import_obs(dict_var[var][0], domain, 'CRU', 'MAM')
	lat, lon, cru_jja = import_obs(dict_var[var][0], domain, 'CRU', 'JJA')
	lat, lon, cru_son = import_obs(dict_var[var][0], domain, 'CRU', 'SON')

	lat, lon, cpc_djf = import_obs(dict_var[var][1], domain, 'CPC', 'DJF')
	lat, lon, cpc_mam = import_obs(dict_var[var][1], domain, 'CPC', 'MAM')
	lat, lon, cpc_jja = import_obs(dict_var[var][1], domain, 'CPC', 'JJA')
	lat, lon, cpc_son = import_obs(dict_var[var][1], domain, 'CPC', 'SON')

	lat, lon, regcm_djf = import_rcm(var, domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_rcm(var, domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_rcm(var, domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_rcm(var, domain, 'RegCM5', 'SON')
	
	mbe_djf_regcm_cru = compute_mbe(regcm_djf[0], cru_djf)
	mbe_mam_regcm_cru = compute_mbe(regcm_mam[0], cru_mam)
	mbe_jja_regcm_cru = compute_mbe(regcm_jja[0], cru_jja)
	mbe_son_regcm_cru = compute_mbe(regcm_son[0], cru_son)	

	mbe_djf_regcm_cpc = compute_mbe(regcm_djf[0], cpc_djf)
	mbe_mam_regcm_cpc = compute_mbe(regcm_mam[0], cpc_mam)
	mbe_jja_regcm_cpc = compute_mbe(regcm_jja[0], cpc_jja)
	mbe_son_regcm_cpc = compute_mbe(regcm_son[0], cpc_son)	

elif var == 'tasmin':
	lat, lon, cru_djf = import_obs(dict_var[var][0], domain, 'CRU', 'DJF')
	lat, lon, cru_mam = import_obs(dict_var[var][0], domain, 'CRU', 'MAM')
	lat, lon, cru_jja = import_obs(dict_var[var][0], domain, 'CRU', 'JJA')
	lat, lon, cru_son = import_obs(dict_var[var][0], domain, 'CRU', 'SON')

	lat, lon, cpc_djf = import_obs(dict_var[var][1], domain, 'CPC', 'DJF')
	lat, lon, cpc_mam = import_obs(dict_var[var][1], domain, 'CPC', 'MAM')
	lat, lon, cpc_jja = import_obs(dict_var[var][1], domain, 'CPC', 'JJA')
	lat, lon, cpc_son = import_obs(dict_var[var][1], domain, 'CPC', 'SON')

	lat, lon, regcm_djf = import_rcm(var, domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_rcm(var, domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_rcm(var, domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_rcm(var, domain, 'RegCM5', 'SON')
	
	mbe_djf_regcm_cru = compute_mbe(regcm_djf[0], cru_djf)
	mbe_mam_regcm_cru = compute_mbe(regcm_mam[0], cru_mam)
	mbe_jja_regcm_cru = compute_mbe(regcm_jja[0], cru_jja)
	mbe_son_regcm_cru = compute_mbe(regcm_son[0], cru_son)	

	mbe_djf_regcm_cpc = compute_mbe(regcm_djf[0], cpc_djf)
	mbe_mam_regcm_cpc = compute_mbe(regcm_mam[0], cpc_mam)
	mbe_jja_regcm_cpc = compute_mbe(regcm_jja[0], cpc_jja)
	mbe_son_regcm_cpc = compute_mbe(regcm_son[0], cpc_son)	
	
elif var == 'clt':
	lat, lon, cru_djf = import_obs(dict_var[var][0], domain, 'CRU', 'DJF')
	lat, lon, cru_mam = import_obs(dict_var[var][0], domain, 'CRU', 'MAM')
	lat, lon, cru_jja = import_obs(dict_var[var][0], domain, 'CRU', 'JJA')
	lat, lon, cru_son = import_obs(dict_var[var][0], domain, 'CRU', 'SON')

	lat, lon, era5_djf = import_obs(dict_var[var][1], domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_obs(dict_var[var][1], domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_obs(dict_var[var][1], domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_obs(dict_var[var][1], domain, 'ERA5', 'SON')

	lat, lon, regcm_djf = import_rcm(var, domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_rcm(var, domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_rcm(var, domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_rcm(var, domain, 'RegCM5', 'SON')
	
	mbe_djf_regcm_cru = compute_mbe(regcm_djf, cru_djf)
	mbe_mam_regcm_cru = compute_mbe(regcm_mam, cru_mam)
	mbe_jja_regcm_cru = compute_mbe(regcm_jja, cru_jja)
	mbe_son_regcm_cru = compute_mbe(regcm_son, cru_son)	
		
	mbe_djf_regcm_era5 = compute_mbe(regcm_djf, era5_djf)
	mbe_mam_regcm_era5 = compute_mbe(regcm_mam, era5_mam)
	mbe_jja_regcm_era5 = compute_mbe(regcm_jja, era5_jja)
	mbe_son_regcm_era5 = compute_mbe(regcm_son, era5_son)
	
else:
	lat, lon, era5_djf = import_obs(dict_var[var][0], domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_obs(dict_var[var][0], domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_obs(dict_var[var][0], domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_obs(dict_var[var][0], domain, 'ERA5', 'SON')

	lat, lon, regcm_djf = import_rcm(var, domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_rcm(var, domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_rcm(var, domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_rcm(var, domain, 'RegCM5', 'SON')
		
	mbe_djf_regcm_era5 = compute_mbe(regcm_djf, era5_djf)
	mbe_mam_regcm_era5 = compute_mbe(regcm_mam, era5_mam)
	mbe_jja_regcm_era5 = compute_mbe(regcm_jja, era5_jja)
	mbe_son_regcm_era5 = compute_mbe(regcm_son, era5_son)

# Plot figure   
dict_plot = {
'pr': ['Bias of  precipitation (mm d$^-$$^1$)', np.arange(-10, 11, 1), cm.BrBG],
'tas': ['Bias of air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'tasmax': ['Bias of maximum air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'tasmin': ['Bias of minimum air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'clt': ['Bias of total cloud cover (%)', np.arange(-70, 80, 10), cm.RdGy],
'rsnl': ['Bias of surface net upward longwave flux (W mm$^-$$^2$)', np.arange(-60, 65, 5), cm.RdBu_r]
}

font_size = 8
	
if var == 'pr':
	fig = plt.figure(figsize=(8, 8))

	ax = fig.add_subplot(4, 4, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(a) RegCM5-CRU DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(b) RegCM5-CPC DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_gpcp[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(c) RegCM5-GPCP DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(d) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(e) RegCM5-CRU MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(f) RegCM5-CPC MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_gpcp[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(g) RegCM5-GPCP MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(h) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 9)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(i) RegCM5-CRU JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 10)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(j) RegCM5-CPC JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 11)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_gpcp[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(k) RegCM5-GPCP JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 12)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(l) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 13)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(m) RegCM5-CRU SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 14)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(n) RegCM5-CPC SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 15)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_gpcp[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(o) RegCM5-GPCP SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 4, 16)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(p) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')

elif var == 'tas':
	fig = plt.figure(figsize=(6, 8))

	ax = fig.add_subplot(4, 2, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(a) RegCM5-CRU DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(b) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(c) RegCM5-CRU MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(d) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(e) RegCM5-CRU JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(f) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(g) RegCM5-CRU SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(h) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')

elif var == 'tasmax':
	fig = plt.figure(figsize=(6, 8))

	ax = fig.add_subplot(4, 2, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(a) RegCM5-CRU DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(b) RegCM5-CPC DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(c) RegCM5-CRU MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(d) RegCM5-CPC MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(e) RegCM5-CRU JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(f) RegCM5-CPC JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(g) RegCM5-CRU SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(h) RegCM5-CPC SON', loc='left', fontsize=font_size, fontweight='bold')

elif var == 'tasmin':
	fig = plt.figure(figsize=(6, 8))
	
	ax = fig.add_subplot(4, 2, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(a) RegCM5-CRU DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(b) RegCM5-CPC DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(c) RegCM5-CRU MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(d) RegCM5-CPC MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(e) RegCM5-CRU JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(f) RegCM5-CPC JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(g) RegCM5-CRU SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(h) RegCM5-CPC SON', loc='left', fontsize=font_size, fontweight='bold')

elif var == 'clt':
	fig = plt.figure(figsize=(6, 8))
	
	ax = fig.add_subplot(4, 2, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(a) RegCM5-CRU DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(b) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(c) RegCM5-CRU MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(d) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(e) RegCM5-CRU JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(f) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(g) RegCM5-CRU SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(h) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')

else:
	fig = plt.figure(figsize=(4, 8))

	ax = fig.add_subplot(4, 1, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(a) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 1, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(b) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 1, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(c) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 1, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
	plt.title(u'(d) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')

# Set colobar
if var == 'rsnl':
	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.84, 0.3, 0.03, 0.4]))
	cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)
else:
	cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.94, 0.3, 0.015, 0.4]))
	cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
	cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/user/mdasilva/CORDEX/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
