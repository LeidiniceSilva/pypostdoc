# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot clim maps"

import os
import cmocean
import netCDF4
import numpy as np
import xarray as xr
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import maskoceans

var = 'pr'
season = 'JJA'
dt = '2018-2021'
path = '/marconi/home/userexternal/mdasilva'

skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,7793,100,105,106,107,112,117,124,135,137,139,
149,152,155,158,168,174,177,183,186,199,204,210,212,224,226,239,240,248,249,253,254,276,277,280,293,298,
303,305,306,308,319,334,335,341,343,359,362,364,384,393,396,398,399,400,402,413,416,417,422,423,426,427,
443,444,446,451,453,457,458,467,474,479,483,488,489,490,495,505,509,513,514,516,529,534,544,559,566]
	
	
def import_ws(param):
	
	yy, xx, mean = [], [], [] 
	for station in range(1, 567):
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		yy.append(inmet[station][2])
		xx.append(inmet[station][3])

		arq  = xr.open_dataset('{0}/OBS/BDMET/database/nc/hourly/{1}/'.format(path, param) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(param, inmet[station][0]))
		data = arq[param]
		time = data.sel(time=slice('2018-06-01','2021-05-31'))
		var  = time.groupby('time.season').mean(dim='time')
		mean.append(var.values*24)
		
	return yy, xx, mean
	
		
def import_obs(param, domain, dataset, season, dt):

	arq   = '{0}/user/mdasilva/{1}/post_evaluate/obs/{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, domain, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

	return lat, lon, mean


def import_sam_3km(param, domain, dataset, season, dt):

	arq   = '{0}/user/mdasilva/{1}/post_evaluate/rcm/{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, domain, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

	return lat, lon, mean
	

def import_csam_4i(param, domain, dataset, season, dt):

	arq   = '{0}/user/mdasilva/{1}/post_evaluate/{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, domain, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

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
dict_var = {'pr': ['pre', 'precip', 'sat_gauge_precip', 'tp']}

if season == 'DJF':
	lat_i, lon_i, inmet_ = import_ws(dict_var[var][0])
	inmet_i = []
	for i in range(0, 298):
		inmet_i.append(inmet_[i][0])
	lat, lon, cru = import_obs(dict_var[var][0], 'SAM-3km', 'CRU', 'DJF', dt)
	lat, lon, cpc = import_obs(dict_var[var][1], 'SAM-3km', 'CPC', 'DJF', dt)
	lat, lon, gpcp = import_obs(dict_var[var][2], 'SAM-3km', 'GPCP', 'DJF', dt)
	lat, lon, era5 = import_obs(dict_var[var][3], 'SAM-3km', 'ERA5', 'DJF', dt)
	lat, lon, wrf_ucan = import_csam_4i(var, 'CSAM-4i', 'ECMWF-ERA5_evaluation_UCAN-WRF433', 'DJF', dt)
	lat, lon, wrf_ncar = import_csam_4i(var, 'CSAM-4i', 'ECMWF-ERA5_evaluation_NCAR-WRF415', 'DJF', dt)
	lat, lon, reg_usp = import_csam_4i(var, 'CSAM-4i', 'ECMWF-ERA5_evaluation_USP-RegCM471', 'DJF', dt)
	lat, lon, reg_ictp_4km = import_csam_4i(var, 'CSAM-4i', 'ECMWF-ERA5_evaluation_ICTP-RegCM5pbl1', 'DJF', dt)
	lat, lon, reg_ictp_3km = import_sam_3km(var, 'SAM-3km', 'RegCM5', 'DJF', dt)
else:
	lat_i, lon_i, inmet_ = import_ws(dict_var[var][0])
	inmet_i = []
	for i in range(0, 298):
		inmet_i.append(inmet_[i][1])
	lat, lon, cru = import_obs(dict_var[var][0], 'SAM-3km', 'CRU', 'JJA', dt)
	lat, lon, cpc = import_obs(dict_var[var][1], 'SAM-3km', 'CPC', 'JJA', dt)
	lat, lon, gpcp = import_obs(dict_var[var][2], 'SAM-3km', 'GPCP', 'JJA', dt)
	lat, lon, era5 = import_obs(dict_var[var][3], 'SAM-3km', 'ERA5', 'JJA', dt)
	lat, lon, wrf_ucan = import_csam_4i(var, 'CSAM-4i', 'ECMWF-ERA5_evaluation_UCAN-WRF433', 'JJA', dt)
	lat, lon, wrf_ncar = import_csam_4i(var, 'CSAM-4i', 'ECMWF-ERA5_evaluation_NCAR-WRF415', 'JJA', dt)
	lat, lon, reg_usp = import_csam_4i(var, 'CSAM-4i', 'ECMWF-ERA5_evaluation_USP-RegCM471', 'JJA', dt)
	lat, lon, reg_ictp_4km = import_csam_4i(var, 'CSAM-4i', 'ECMWF-ERA5_evaluation_ICTP-RegCM5pbl1', 'JJA', dt)
	lat, lon, reg_ictp_3km = import_sam_3km(var, 'SAM-3km', 'RegCM5', 'JJA', dt)
	
# Plot figure
color = ['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']
dict_plot = {'pr': ['Precipitation (mm d$^-$$^1$)', np.arange(0, 18, 1), matplotlib.colors.ListedColormap(color)]}
font_size = 8
	
fig = plt.figure(figsize=(10, 3))

ax = fig.add_subplot(2, 5, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.scatter(lon_i, lat_i, 4, inmet_i, cmap=dict_plot[var][2], marker='o', vmin=0, vmax=18) 
plt.title(u'(a) INMET', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, cru, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(b) CRU', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, cpc, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(c) CPC', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 4)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, gpcp, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(d) GPCP', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 5)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, era5, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(e) ERA5', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 6)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, reg_ictp_3km, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(f) RegCM5 3km', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 7)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, reg_ictp_4km, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(g) RegCM5 4km', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 8)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, reg_usp, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(h) RegCM4', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 9)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, wrf_ncar, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(i) WRF-NCAR', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 10)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, wrf_ucan, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='max') 
plt.title(u'(j) WRF-UCAN', loc='left', fontsize=font_size, fontweight='bold')

# Set colobar
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.3, 0.01, 0.4]))
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cp_3km-4km'.format(path)
name_out = 'pyplt_maps_clim_{0}_RCM_SAM-3km_{1}_{2}.png'.format(var, season, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
