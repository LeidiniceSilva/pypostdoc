# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot bias maps"

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
from import_climate_tools import compute_mbe

var = 'pr'
season = 'JJA'
dt = '2018-2021'
path = '/marconi/home/userexternal/mdasilva'

skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,7793,100,105,106,107,112,117,124,135,137,139,
149,152,155,158,168,174,177,183,186,199,204,210,212,224,226,239,240,248,249,253,254,276,277,280,293,298,
303,305,306,308,319,334,335,341,343,359,362,364,384,393,396,398,399,400,402,413,416,417,422,423,426,427,
443,444,446,451,453,457,458,467,474,479,483,488,489,490,495,505,509,513,514,516,529,534,544,559,566]
	
	
def import_ws(param_i, param_ii):
	
	yy, xx = [], []
	mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi = [], [], [], [], [], []
	
	for station in range(1, 567):
		print(station, inmet[station][0])
		
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		yy.append(inmet[station][2])
		xx.append(inmet[station][3])

		arq_i  = xr.open_dataset('{0}/OBS/BDMET/database/nc/hourly/{1}/'.format(path, param_i) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(param_i, inmet[station][0]))
		data_i = arq_i[param_i]
		time_i = data_i.sel(time=slice('2018-06-01','2021-05-31'))
		var_i  = time_i.groupby('time.season').mean(dim='time')
		mean_i.append(var_i.values*24)

		arq_ii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/'.format(path) + '{0}_SAM-3km_RegCM5_{1}_{2}_lonlat.nc'.format(param_ii, season, dt))
		data_ii = arq_ii[param_ii]
		data_ii = data_ii.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		time_ii = data_ii.sel(time=slice('2018-06-01','2021-05-31'))
		var_ii  = time_ii.values
		mean_ii.append(var_ii)

		arq_iii  = xr.open_dataset('{0}/user/mdasilva/CSAM-4i/post_evaluate/'.format(path) + '{0}_CSAM-4i_ECMWF-ERA5_evaluation_ICTP-RegCM5pbl1_{1}_{2}_lonlat.nc'.format(param_ii, season, dt))
		data_iii = arq_iii[param_ii]
		data_iii = data_iii.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		time_iii = data_iii.sel(time=slice('2018-06-01','2021-05-31'))
		var_iii  = time_iii.values
		mean_iii.append(var_iii)

		arq_iv  = xr.open_dataset('{0}/user/mdasilva/CSAM-4i/post_evaluate/'.format(path) + '{0}_CSAM-4i_ECMWF-ERA5_evaluation_USP-RegCM471_{1}_{2}_lonlat.nc'.format(param_ii, season, dt))
		data_iv = arq_iv[param_ii]
		data_iv = data_iv.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		time_iv = data_iv.sel(time=slice('2018-06-01','2021-05-31'))
		var_iv  = time_iv.values
		mean_iv.append(var_iv)

		arq_v  = xr.open_dataset('{0}/user/mdasilva/CSAM-4i/post_evaluate/'.format(path) + '{0}_CSAM-4i_ECMWF-ERA5_evaluation_NCAR-WRF415_{1}_{2}_lonlat.nc'.format(param_ii, season, dt))
		data_v = arq_v[param_ii]
		data_v = data_v.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		time_v = data_v.sel(time=slice('2018-06-01','2021-05-31'))
		var_v  = time_v.values
		mean_v.append(var_v)

		arq_vi  = xr.open_dataset('{0}/user/mdasilva/CSAM-4i/post_evaluate/'.format(path) + '{0}_CSAM-4i_ECMWF-ERA5_evaluation_UCAN-WRF433_{1}_{2}_lonlat.nc'.format(param_ii, season, dt))
		data_vi = arq_vi[param_ii]
		data_vi = data_vi.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		time_vi = data_vi.sel(time=slice('2018-06-01','2021-05-31'))
		var_vi  = time_vi.values
		mean_vi.append(var_vi)
								
	return yy, xx, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi
	
		
def import_obs(param, domain, dataset, season):

	arq   = '{0}/user/mdasilva/{1}/post_evaluate/obs/{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, domain, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

	return lat, lon, mean


def import_sam_3km(param, domain, dataset, season):

	arq   = '{0}/user/mdasilva/{1}/post_evaluate/rcm/{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, domain, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

	return lat, lon, mean
	

def import_csam_4i(param, domain, dataset, season):

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

lat_i, lon_i, inmet_, reg_ictp_3km_, reg_ictp_4km_, reg_usp_, wrf_ncar_, wrf_ucan_ = import_ws(dict_var[var][0], var)
mbe_reg_ictp_3km_inmet, mbe_reg_ictp_4km_inmet, mbe_reg_usp_inmet, mbe_wrf_ncar_inmet, mbe_wrf_ucan_inmet = [], [], [], [], []

if season == 'DJF':
	idx = 0
else:
	idx = 1
	
for i in range(0, 298):
	mbe_reg_ictp_3km_inmet.append(compute_mbe(reg_ictp_3km_[i], inmet_[i][idx]))
	mbe_reg_ictp_4km_inmet.append(compute_mbe(reg_ictp_4km_[i], inmet_[i][idx]))
	mbe_reg_usp_inmet.append(compute_mbe(reg_usp_[i], inmet_[i][idx]))
	mbe_wrf_ncar_inmet.append(compute_mbe(wrf_ncar_[i], inmet_[i][idx]))
	mbe_wrf_ucan_inmet.append(compute_mbe(wrf_ucan_[i], inmet_[i][idx]))
		
lat, lon, cru = import_obs(dict_var[var][0], 'SAM-3km', 'CRU', season)
lat, lon, cpc = import_obs(dict_var[var][1], 'SAM-3km', 'CPC', season)
lat, lon, gpcp = import_obs(dict_var[var][2], 'SAM-3km', 'GPCP', season)
lat, lon, era5 = import_obs(dict_var[var][3], 'SAM-3km', 'ERA5', season)
lat, lon, wrf_ucan = import_csam_4i(var, 'CSAM-4i', 'ECMWF-ERA5_evaluation_UCAN-WRF433', season)
lat, lon, wrf_ncar = import_csam_4i(var, 'CSAM-4i', 'ECMWF-ERA5_evaluation_NCAR-WRF415', season)
lat, lon, reg_usp = import_csam_4i(var, 'CSAM-4i', 'ECMWF-ERA5_evaluation_USP-RegCM471', season)
lat, lon, reg_ictp_4km = import_csam_4i(var, 'CSAM-4i', 'ECMWF-ERA5_evaluation_ICTP-RegCM5pbl1', season)
lat, lon, reg_ictp_3km = import_sam_3km(var, 'SAM-3km', 'RegCM5', season)
	
mbe_reg_ictp_3km_cru = compute_mbe(reg_ictp_3km, cru)	
mbe_reg_ictp_4km_cru = compute_mbe(reg_ictp_4km, cru)	
mbe_reg_usp_cru = compute_mbe(reg_usp, cru)	
mbe_wrf_ncar_cru = compute_mbe(wrf_ncar, cru)	
mbe_wrf_ucan_cru = compute_mbe(wrf_ucan, cru)	

mbe_reg_ictp_3km_cpc = compute_mbe(reg_ictp_3km, cpc)	
mbe_reg_ictp_4km_cpc = compute_mbe(reg_ictp_4km, cpc)	
mbe_reg_usp_cpc = compute_mbe(reg_usp, cpc)	
mbe_wrf_ncar_cpc = compute_mbe(wrf_ncar, cpc)	
mbe_wrf_ucan_cpc = compute_mbe(wrf_ucan, cpc)	

mbe_reg_ictp_3km_gpcp = compute_mbe(reg_ictp_3km, gpcp)	
mbe_reg_ictp_4km_gpcp = compute_mbe(reg_ictp_4km, gpcp)	
mbe_reg_usp_gpcp = compute_mbe(reg_usp, gpcp)	
mbe_wrf_ncar_gpcp = compute_mbe(wrf_ncar, gpcp)	
mbe_wrf_ucan_gpcp = compute_mbe(wrf_ucan, gpcp)	

mbe_reg_ictp_3km_era5 = compute_mbe(reg_ictp_3km, era5)	
mbe_reg_ictp_4km_era5 = compute_mbe(reg_ictp_4km, era5)	
mbe_reg_usp_era5 = compute_mbe(reg_usp, era5)	
mbe_wrf_ncar_era5 = compute_mbe(wrf_ncar, era5)  
mbe_wrf_ucan_era5 = compute_mbe(wrf_ucan, era5)	

# Plot figure
dict_plot = {'pr': ['Bias of precipitation (mm d$^-$$^1$)', np.arange(-10, 11, 1), cm.BrBG]}
font_size = 8
	
fig = plt.figure(figsize=(10, 8))

ax = fig.add_subplot(5, 5, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.scatter(lon_i, lat_i, 4, mbe_reg_ictp_3km_inmet, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=11) 
plt.title(u'(a) RegCM5-3km - INMET', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.scatter(lon_i, lat_i, 4, mbe_reg_ictp_4km_inmet, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=11) 
plt.title(u'(b) RegCM5-4km - INMET', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.scatter(lon_i, lat_i, 4, mbe_reg_usp_inmet, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=11) 
plt.title(u'(c) RegCM4 - INMET', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 4)  
map, xx, yy = basemap(lat, lon)
plt_map = map.scatter(lon_i, lat_i, 4, mbe_wrf_ncar_inmet, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=11) 
plt.title(u'(d) WRF-NCAR - INMET', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 5)  
map, xx, yy = basemap(lat, lon)
plt_map = map.scatter(lon_i, lat_i, 4, mbe_wrf_ucan_inmet, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=11) 
plt.title(u'(e) WRF-UCAN - INMET', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 6)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_reg_ictp_3km_cru, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(f) RegCM5-3km - CRU', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 7)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_reg_ictp_4km_cru, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(g) RegCM5-4km - CRU', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 8)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_reg_usp_cru, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(h) RegCM4 - CRU', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 9)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wrf_ncar_cru, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(i) WRF-NCAR - CRU', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 10)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wrf_ucan_cru, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(j) WRF-UCAN- CRU', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 11)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_reg_ictp_3km_cpc, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(k) RegCM5-3km - CPC', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 12)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_reg_ictp_4km_cpc, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(l) RegCM5-4km - CPC', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 13)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_reg_usp_cpc, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(m) RegCM4 - CPC', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 14)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wrf_ncar_cpc, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(n) WRF-NCAR - CPC', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 15)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wrf_ucan_cpc, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(o) WRF-UCAN- CPC', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 16)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_reg_ictp_3km_gpcp, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(p) RegCM5-3km - GPCP', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 17)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_reg_ictp_4km_gpcp, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(q) RegCM5-4km - GPCP', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 18)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_reg_usp_gpcp, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(r) RegCM4 - GPCP', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 19)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wrf_ncar_gpcp, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(s) WRF-NCAR - GPCP', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 20)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wrf_ucan_gpcp, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(t) WRF-UCAN- GPCP', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 21)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_reg_ictp_3km_era5, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(u) RegCM5-3km - ERA5', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 22)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_reg_ictp_4km_era5, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(v) RegCM5-4km - ERA5', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 23)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_reg_usp_era5, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(w) RegCM4 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 24)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wrf_ncar_era5, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(x) WRF-NCAR - ERA5', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 5, 25)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_wrf_ucan_era5, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
plt.title(u'(y) WRF-UCAN- ERA5', loc='left', fontsize=font_size, fontweight='bold')

# Set colobar
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.3, 0.01, 0.4]))
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cp_3km-4km'.format(path)
name_out = 'pyplt_maps_bias_{0}_RCM_SAM-3km_{1}_{2}.png'.format(var, season, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
