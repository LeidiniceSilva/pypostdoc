# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot p99 bias maps"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap
from import_climate_tools import compute_mbe

var = 'p99'
freq = 'hourly'
domain = 'SAM-3km'
path = '/marconi/home/userexternal/mdasilva'


if freq == 'hourly':
	dt = '1hr_2018-2021'
else:
	dt = '2018-2021'


def import_obs(param, domain, dataset):

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/obs/p99_{1}_{2}_{3}_lonlat.nc'.format(path, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	
		
def import_sam_3km(param, domain, dataset):

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/p99_{1}_{2}_{3}_lonlat.nc'.format(path, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean


def import_csam_4i(param, domain, dataset):

	arq   = '{0}/user/mdasilva/CSAM-4i/post_evaluate/p99_{1}_{2}_{3}_lonlat.nc'.format(path, domain, dataset, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]

	return lat, lon, mean
	

def basemap(lat, lon):

	map = Basemap(projection='cyl', llcrnrlon=-76., llcrnrlat=-35., urcrnrlon=-38.,urcrnrlat=-12., resolution='c')
	map.drawmeridians(np.arange(-76., -38., 10.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-35., -12., 5.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black')
	
	lons, lats = np.meshgrid(lon, lat)
	xx, yy = map(lons,lats)
	
	return map, xx, yy
	
	
# Import model and obs dataset
if freq == 'hourly':
	lat, lon, era5 = import_obs('tp', 'SAM-3km', 'ERA5')

	lat, lon, wrf_ucan = import_csam_4i('pr', 'CSAM-4i', 'ECMWF-ERA5_evaluation_UCAN-WRF433')
	lat, lon, wrf_ncar = import_csam_4i('pr', 'CSAM-4i', 'ECMWF-ERA5_evaluation_NCAR-WRF415')
	lat, lon, reg_usp = import_csam_4i('pr', 'CSAM-4i', 'ECMWF-ERA5_evaluation_USP-RegCM471')
	lat, lon, reg_ictp_4km = import_csam_4i('pr', 'CSAM-4i', 'ECMWF-ERA5_evaluation_ICTP-RegCM5pbl1')
	lat, lon, reg_ictp_3km = import_sam_3km('pr', 'SAM-3km', 'RegCM5')
	
	mbe_reg_ictp_3km_era5 = compute_mbe(reg_ictp_3km, era5)	
	mbe_reg_ictp_4km_era5 = compute_mbe(reg_ictp_4km, era5)	
	mbe_reg_usp_era5 = compute_mbe(reg_usp, era5)	
	mbe_wrf_ncar_era5 = compute_mbe(wrf_ncar, era5)  
	mbe_wrf_ucan_era5 = compute_mbe(wrf_ucan, era5)	
else:
	lat, lon, cpc = import_obs('precip', 'SAM-3km', 'CPC')
	lat, lon, era5 = import_obs('tp', 'SAM-3km', 'ERA5')

	lat, lon, wrf_ucan = import_csam_4i('pr', 'CSAM-4i', 'ECMWF-ERA5_evaluation_UCAN-WRF433')
	lat, lon, wrf_ncar = import_csam_4i('pr', 'CSAM-4i', 'ECMWF-ERA5_evaluation_NCAR-WRF415')
	lat, lon, reg_usp = import_csam_4i('pr', 'CSAM-4i', 'ECMWF-ERA5_evaluation_USP-RegCM471')
	lat, lon, reg_ictp_4km = import_csam_4i('pr', 'CSAM-4i', 'ECMWF-ERA5_evaluation_ICTP-RegCM5pbl1')
	lat, lon, reg_ictp_3km = import_sam_3km('pr', 'SAM-3km', 'RegCM5')
	
	mbe_reg_ictp_3km_cpc = compute_mbe(reg_ictp_3km, cpc)	
	mbe_reg_ictp_4km_cpc = compute_mbe(reg_ictp_4km, cpc)	
	mbe_reg_usp_cpc = compute_mbe(reg_usp, cpc)	
	mbe_wrf_ncar_cpc = compute_mbe(wrf_ncar, cpc)	
	mbe_wrf_ucan_cpc = compute_mbe(wrf_ucan, cpc)	

	mbe_reg_ictp_3km_era5 = compute_mbe(reg_ictp_3km, era5)	
	mbe_reg_ictp_4km_era5 = compute_mbe(reg_ictp_4km, era5)	
	mbe_reg_usp_era5 = compute_mbe(reg_usp, era5)	
	mbe_wrf_ncar_era5 = compute_mbe(wrf_ncar, era5)  
	mbe_wrf_ucan_era5 = compute_mbe(wrf_ucan, era5)	

# Plot figure
fig = plt.figure(figsize=(10, 3))   
font_size = 8

if freq == 'hourly':
	levs = np.arange(-5, 5.5, 0.5)
	legend = 'Hourly p99 (mm h$^-$$^1$)'

	ax = fig.add_subplot(1, 5, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_reg_ictp_3km_era5[0], levels=np.arange(-5, 5.5, 0.5), cmap=cm.BrBG, extend='both') 
	plt.title(u'(a) RegCM5-3km - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(1, 5, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_reg_ictp_4km_era5[0], levels=np.arange(-5, 5.5, 0.5), cmap=cm.BrBG, extend='both') 
	plt.title(u'(b) RegCM5-4km - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(1, 5, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_reg_usp_era5[0], levels=np.arange(-5, 5.5, 0.5), cmap=cm.BrBG, extend='both') 
	plt.title(u'(c) RegCM4 - CPC', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(1, 5, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_wrf_ncar_era5[0], levels=np.arange(-5, 5.5, 0.5), cmap=cm.BrBG, extend='both') 
	plt.title(u'(d) WRF-NCAR - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(1, 5, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_wrf_ucan_era5[0], levels=np.arange(-5, 5.5, 0.5), cmap=cm.BrBG, extend='both') 
	plt.title(u'(e) WRF-UCAN - ERA5', loc='left', fontsize=font_size, fontweight='bold')
else:
	levs = np.arange(-60, 70, 10)
	legend = 'Daily p99 (mm d$^-$$^1$)'

	ax = fig.add_subplot(2, 5, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_reg_ictp_3km_cpc[0], levels=levs, cmap=cm.BrBG, extend='both') 
	plt.title(u'(a) RegCM5-3km - CPC', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(2, 5, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_reg_ictp_4km_cpc[0], levels=levs, cmap=cm.BrBG, extend='both') 
	plt.title(u'(b) RegCM5-4km - CPC', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(2, 5, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_reg_usp_cpc[0], levels=levs, cmap=cm.BrBG, extend='both') 
	plt.title(u'(c) RegCM4 - CPC', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(2, 5, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_wrf_ncar_cpc[0], levels=levs, cmap=cm.BrBG, extend='both') 
	plt.title(u'(d) WRF-NCAR - CPC', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(2, 5, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_wrf_ucan_cpc[0], levels=levs, cmap=cm.BrBG, extend='both') 
	plt.title(u'(e) WRF-UCAN - CPC', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(2, 5, 6)
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_reg_ictp_3km_era5[0], levels=levs, cmap=cm.BrBG, extend='both')
	plt.title(u'(f) RegCM5-3km - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(2, 5, 7)
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_reg_ictp_4km_era5[0], levels=levs, cmap=cm.BrBG, extend='both')
	plt.title(u'(g) RegCM5-4km - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(2, 5, 8)
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_reg_usp_era5[0], levels=levs, cmap=cm.BrBG, extend='both')
	plt.title(u'(h) RegCM4 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(2, 5, 9)
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_wrf_ncar_era5[0], levels=levs, cmap=cm.BrBG, extend='both')
	plt.title(u'(i) WRF-NCAR - ERA5', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(2, 5, 10)
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_wrf_ucan_era5[0], levels=levs, cmap=cm.BrBG, extend='both')
	plt.title(u'(j) WRF-UCAN - ERA5', loc='left', fontsize=font_size, fontweight='bold')
		
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.3, 0.01, 0.4]))
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cp_3km-4km'.format(path)
name_out = 'pyplt_maps_bias_{0}_RegCM5_{1}_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
