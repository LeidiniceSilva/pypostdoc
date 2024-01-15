# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot bias maps"

import os
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from import_dataset_situ import import_situ
from import_dataset_grid import import_grid
from import_climate_tools import compute_mbe
from import_basemap import basemap

path='/marconi/home/userexternal/mdasilva'

domain = 'SAM-3km'

var = 'tasmax'

# Import model and obs dataset 
if var == 'pr':
	dict_var = {'pr': ['pre', 'pre', 'precip', 'sat_gauge_precip', 'tp']}

	lat_i, lon_i, inmet_i, regcm_i = import_situ(dict_var[var][0], var, domain, 'RegCM5')
	mbe_djf_regcm_inmet, mbe_mam_regcm_inmet, mbe_jja_regcm_inmet, mbe_son_regcm_inmet = [], [], [], []
	for i in range(0, 298):
		mbe_djf_regcm_inmet.append(compute_mbe(regcm_i[i][0], inmet_i[i][0]))
		mbe_mam_regcm_inmet.append(compute_mbe(regcm_i[i][1], inmet_i[i][1]))
		mbe_jja_regcm_inmet.append(compute_mbe(regcm_i[i][2], inmet_i[i][2]))
		mbe_son_regcm_inmet.append(compute_mbe(regcm_i[i][3], inmet_i[i][3]))
	
	lat, lon, cru_djf = import_grid(dict_var[var][1], domain, 'CRU', 'DJF')
	lat, lon, cru_mam = import_grid(dict_var[var][1], domain, 'CRU', 'MAM')
	lat, lon, cru_jja = import_grid(dict_var[var][1], domain, 'CRU', 'JJA')
	lat, lon, cru_son = import_grid(dict_var[var][1], domain, 'CRU', 'SON')

	lat, lon, cpc_djf = import_grid(dict_var[var][2], domain, 'CPC', 'DJF')
	lat, lon, cpc_mam = import_grid(dict_var[var][2], domain, 'CPC', 'MAM')
	lat, lon, cpc_jja = import_grid(dict_var[var][2], domain, 'CPC', 'JJA')
	lat, lon, cpc_son = import_grid(dict_var[var][2], domain, 'CPC', 'SON')

	lat, lon, gpcp_djf = import_grid(dict_var[var][3], domain, 'GPCP', 'DJF')
	lat, lon, gpcp_mam = import_grid(dict_var[var][3], domain, 'GPCP', 'MAM')
	lat, lon, gpcp_jja = import_grid(dict_var[var][3], domain, 'GPCP', 'JJA')
	lat, lon, gpcp_son = import_grid(dict_var[var][3], domain, 'GPCP', 'SON')

	lat, lon, era5_djf = import_grid(dict_var[var][4], domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_grid(dict_var[var][4], domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_grid(dict_var[var][4], domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_grid(dict_var[var][4], domain, 'ERA5', 'SON')

	lat, lon, regcm_djf = import_grid(var, domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_grid(var, domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_grid(var, domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_grid(var, domain, 'RegCM5', 'SON')
	
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
	dict_var = {'tas': ['tmp', 'tmp', 't2m']}

	lat_i, lon_i, inmet_i, regcm_i = import_situ(dict_var[var][0], var, domain, 'RegCM5')
	mbe_djf_regcm_inmet, mbe_mam_regcm_inmet, mbe_jja_regcm_inmet, mbe_son_regcm_inmet = [], [], [], []
	for i in range(0, 298):
		mbe_djf_regcm_inmet.append(compute_mbe(regcm_i[i][0], inmet_i[i][0]))
		mbe_mam_regcm_inmet.append(compute_mbe(regcm_i[i][1], inmet_i[i][1]))
		mbe_jja_regcm_inmet.append(compute_mbe(regcm_i[i][2], inmet_i[i][2]))
		mbe_son_regcm_inmet.append(compute_mbe(regcm_i[i][3], inmet_i[i][3]))
		
	lat, lon, cru_djf = import_grid(dict_var[var][1], domain, 'CRU', 'DJF')
	lat, lon, cru_mam = import_grid(dict_var[var][1], domain, 'CRU', 'MAM')
	lat, lon, cru_jja = import_grid(dict_var[var][1], domain, 'CRU', 'JJA')
	lat, lon, cru_son = import_grid(dict_var[var][1], domain, 'CRU', 'SON')

	lat, lon, cpc_djf = import_grid(dict_var[var][2], domain, 'CPC', 'DJF')
	lat, lon, cpc_mam = import_grid(dict_var[var][2], domain, 'CPC', 'MAM')
	lat, lon, cpc_jja = import_grid(dict_var[var][2], domain, 'CPC', 'JJA')
	lat, lon, cpc_son = import_grid(dict_var[var][2], domain, 'CPC', 'SON')

	lat, lon, era5_djf = import_grid(dict_var[var][4], domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_grid(dict_var[var][4], domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_grid(dict_var[var][4], domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_grid(dict_var[var][4], domain, 'ERA5', 'SON')

	lat, lon, regcm_djf = import_grid(var, domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_grid(var, domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_grid(var, domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_grid(var, domain, 'RegCM5', 'SON')
	
	mbe_djf_regcm_cru = compute_mbe(regcm_djf, cru_djf)
	mbe_mam_regcm_cru = compute_mbe(regcm_mam, cru_mam)
	mbe_jja_regcm_cru = compute_mbe(regcm_jja, cru_jja)
	mbe_son_regcm_cru = compute_mbe(regcm_son, cru_son)	

	mbe_djf_regcm_cpc = compute_mbe(regcm_djf, cpc_djf)
	mbe_mam_regcm_cpc = compute_mbe(regcm_mam, cpc_mam)
	mbe_jja_regcm_cpc = compute_mbe(regcm_jja, cpc_jja)
	mbe_son_regcm_cpc = compute_mbe(regcm_son, cpc_son)	
		
	mbe_djf_regcm_era5 = compute_mbe(regcm_djf, era5_djf)
	mbe_mam_regcm_era5 = compute_mbe(regcm_mam, era5_mam)
	mbe_jja_regcm_era5 = compute_mbe(regcm_jja, era5_jja)
	mbe_son_regcm_era5 = compute_mbe(regcm_son, era5_son)
	
elif var == 'tasmax':
	dict_var = {'tasmax': ['tmx', 'tmax', 'mx2t']}
	
	lat, lon, cru_djf = import_grid(dict_var[var][1], domain, 'CRU', 'DJF')
	lat, lon, cru_mam = import_grid(dict_var[var][1], domain, 'CRU', 'MAM')
	lat, lon, cru_jja = import_grid(dict_var[var][1], domain, 'CRU', 'JJA')
	lat, lon, cru_son = import_grid(dict_var[var][1], domain, 'CRU', 'SON')

	lat, lon, cpc_djf = import_grid(dict_var[var][2], domain, 'CPC', 'DJF')
	lat, lon, cpc_mam = import_grid(dict_var[var][2], domain, 'CPC', 'MAM')
	lat, lon, cpc_jja = import_grid(dict_var[var][2], domain, 'CPC', 'JJA')
	lat, lon, cpc_son = import_grid(dict_var[var][2], domain, 'CPC', 'SON')

	lat, lon, era5_djf = import_grid(dict_var[var][4], domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_grid(dict_var[var][4], domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_grid(dict_var[var][4], domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_grid(dict_var[var][4], domain, 'ERA5', 'SON')

	lat, lon, regcm_djf = import_grid(var, domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_grid(var, domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_grid(var, domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_grid(var, domain, 'RegCM5', 'SON')
	
	mbe_djf_regcm_cru = compute_mbe(regcm_djf, cru_djf)
	mbe_mam_regcm_cru = compute_mbe(regcm_mam, cru_mam)
	mbe_jja_regcm_cru = compute_mbe(regcm_jja, cru_jja)
	mbe_son_regcm_cru = compute_mbe(regcm_son, cru_son)	

	mbe_djf_regcm_cpc = compute_mbe(regcm_djf, cpc_djf)
	mbe_mam_regcm_cpc = compute_mbe(regcm_mam, cpc_mam)
	mbe_jja_regcm_cpc = compute_mbe(regcm_jja, cpc_jja)
	mbe_son_regcm_cpc = compute_mbe(regcm_son, cpc_son)	
		
	mbe_djf_regcm_era5 = compute_mbe(regcm_djf, era5_djf)
	mbe_mam_regcm_era5 = compute_mbe(regcm_mam, era5_mam)
	mbe_jja_regcm_era5 = compute_mbe(regcm_jja, era5_jja)
	mbe_son_regcm_era5 = compute_mbe(regcm_son, era5_son)	
	
elif var == 'tasmin':
	dict_var = {'tasmin': ['tmn', 'tmin', 'mn2t']}

	lat, lon, cru_djf = import_grid(dict_var[var][1], domain, 'CRU', 'DJF')
	lat, lon, cru_mam = import_grid(dict_var[var][1], domain, 'CRU', 'MAM')
	lat, lon, cru_jja = import_grid(dict_var[var][1], domain, 'CRU', 'JJA')
	lat, lon, cru_son = import_grid(dict_var[var][1], domain, 'CRU', 'SON')

	lat, lon, cpc_djf = import_grid(dict_var[var][2], domain, 'CPC', 'DJF')
	lat, lon, cpc_mam = import_grid(dict_var[var][2], domain, 'CPC', 'MAM')
	lat, lon, cpc_jja = import_grid(dict_var[var][2], domain, 'CPC', 'JJA')
	lat, lon, cpc_son = import_grid(dict_var[var][2], domain, 'CPC', 'SON')

	lat, lon, era5_djf = import_grid(dict_var[var][4], domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_grid(dict_var[var][4], domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_grid(dict_var[var][4], domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_grid(dict_var[var][4], domain, 'ERA5', 'SON')

	lat, lon, regcm_djf = import_grid(var, domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_grid(var, domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_grid(var, domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_grid(var, domain, 'RegCM5', 'SON')
	
	mbe_djf_regcm_cru = compute_mbe(regcm_djf, cru_djf)
	mbe_mam_regcm_cru = compute_mbe(regcm_mam, cru_mam)
	mbe_jja_regcm_cru = compute_mbe(regcm_jja, cru_jja)
	mbe_son_regcm_cru = compute_mbe(regcm_son, cru_son)	

	mbe_djf_regcm_cpc = compute_mbe(regcm_djf, cpc_djf)
	mbe_mam_regcm_cpc = compute_mbe(regcm_mam, cpc_mam)
	mbe_jja_regcm_cpc = compute_mbe(regcm_jja, cpc_jja)
	mbe_son_regcm_cpc = compute_mbe(regcm_son, cpc_son)	
		
	mbe_djf_regcm_era5 = compute_mbe(regcm_djf, era5_djf)
	mbe_mam_regcm_era5 = compute_mbe(regcm_mam, era5_mam)
	mbe_jja_regcm_era5 = compute_mbe(regcm_jja, era5_jja)
	mbe_son_regcm_era5 = compute_mbe(regcm_son, era5_son)

else:
	dict_var = {'clt': ['cld', 'tcc']}

	lat, lon, cru_djf = import_grid(dict_var[var][1], domain, 'CRU', 'DJF')
	lat, lon, cru_mam = import_grid(dict_var[var][1], domain, 'CRU', 'MAM')
	lat, lon, cru_jja = import_grid(dict_var[var][1], domain, 'CRU', 'JJA')
	lat, lon, cru_son = import_grid(dict_var[var][1], domain, 'CRU', 'SON')

	lat, lon, era5_djf = import_grid(dict_var[var][4], domain, 'ERA5', 'DJF')
	lat, lon, era5_mam = import_grid(dict_var[var][4], domain, 'ERA5', 'MAM')
	lat, lon, era5_jja = import_grid(dict_var[var][4], domain, 'ERA5', 'JJA')
	lat, lon, era5_son = import_grid(dict_var[var][4], domain, 'ERA5', 'SON')

	lat, lon, regcm_djf = import_grid(var, domain, 'RegCM5', 'DJF')
	lat, lon, regcm_mam = import_grid(var, domain, 'RegCM5', 'MAM')
	lat, lon, regcm_jja = import_grid(var, domain, 'RegCM5', 'JJA')
	lat, lon, regcm_son = import_grid(var, domain, 'RegCM5', 'SON')
	
	mbe_djf_regcm_cru = compute_mbe(regcm_djf, cru_djf)
	mbe_mam_regcm_cru = compute_mbe(regcm_mam, cru_mam)
	mbe_jja_regcm_cru = compute_mbe(regcm_jja, cru_jja)
	mbe_son_regcm_cru = compute_mbe(regcm_son, cru_son)	
		
	mbe_djf_regcm_era5 = compute_mbe(regcm_djf, era5_djf)
	mbe_mam_regcm_era5 = compute_mbe(regcm_mam, era5_mam)
	mbe_jja_regcm_era5 = compute_mbe(regcm_jja, era5_jja)
	mbe_son_regcm_era5 = compute_mbe(regcm_son, era5_son)

# Plot figure   
dict_plot = {
'pr': ['Bias of  precipitation (mm d⁻¹)', np.arange(-10, 11, 1), cm.BrBG],
'tas': ['Bias of air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'tasmax': ['Bias of maximum air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'tasmin': ['Bias of minimum air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'clt': ['Bias of total cloud cover (%)', np.arange(-60, 70, 10), cm.RdGy]
}

font_size = 8
	
if var == 'pr':
	fig = plt.figure(figsize=(10, 6))

	ax = fig.add_subplot(4, 5, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.scatter(lon_i, lat_i, 4, mbe_djf_regcm_inmet, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=10) 
	plt.title(u'(a) RegCM5-INMET DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(b) RegCM5-CRU DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(c) RegCM5-CPC DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_gpcp[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(d) RegCM5-GPCP DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(e) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.scatter(lon_i, lat_i, 4, mbe_mam_regcm_inmet, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=10) 
	plt.title(u'(f) RegCM5-INMET MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(g) RegCM5-CRU MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(h) RegCM5-CPC MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 9)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_gpcp[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(i) RegCM5-GPCP MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 10)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(j) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 11)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.scatter(lon_i, lat_i, 4, mbe_jja_regcm_inmet, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=10) 
	plt.title(u'(k) RegCM5-INMET JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 12)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(l) RegCM5-CRU JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 13)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(m) RegCM5-CPC JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 14)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_gpcp[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(n) RegCM5-GPCP JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 15)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(o) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 16)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.scatter(lon_i, lat_i, 4, mbe_son_regcm_inmet, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=10) 
	plt.title(u'(p) RegCM5-INMET SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 17)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(q) RegCM5-CRU SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 18)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(r) RegCM5-CPC SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 19)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_gpcp[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(s) RegCM5-GPCP SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 5, 20)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(t) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')

elif var == 'tas':
	fig = plt.figure(figsize=(8, 8))

	ax = fig.add_subplot(4, 3, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.scatter(lon_i, lat_i, 4, mbe_djf_regcm_inmet, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=10) 
	plt.title(u'(a) RegCM5-INMET DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(b) RegCM5-CRU DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(c) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.scatter(lon_i, lat_i, 4, mbe_mam_regcm_inmet, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=10) 
	plt.title(u'(d) RegCM5-INMET MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(e) RegCM5-CRU MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(f) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.scatter(lon_i, lat_i, 4, mbe_jja_regcm_inmet, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=10) 
	plt.title(u'(g) RegCM5-INMET JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(h) RegCM5-CRU JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 9)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(i) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 10)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.scatter(lon_i, lat_i, 4, mbe_son_regcm_inmet, cmap=dict_plot[var][2], marker='o', vmin=-10, vmax=10) 
	plt.title(u'(j) RegCM5-INMET SON', loc='left', fontsize=font_size, fontweight='bold')
	
	ax = fig.add_subplot(4, 3, 11)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(k) RegCM5-CRU SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 12)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(l) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')

elif var == 'tasmax':
	fig = plt.figure(figsize=(8, 8))

	ax = fig.add_subplot(4, 3, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(a) RegCM5-CRU DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(b) RegCM5-CPC DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(c) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(d) RegCM5-CRU MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(e) RegCM5-CPC MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(f) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(g) RegCM5-CRU JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(h) RegCM5-CPC JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 9)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(i) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 10)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(j) RegCM5-CRU SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 11)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(k) RegCM5-CPC SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 12)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(l) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')

elif var == 'tasmin':
	fig = plt.figure(figsize=(8, 8))

	ax = fig.add_subplot(4, 3, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(a) RegCM5-CRU DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(b) RegCM5-CPC DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(c) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(d) RegCM5-CRU MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(e) RegCM5-CPC MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(f) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(g) RegCM5-CRU JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(h) RegCM5-CPC JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 9)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(i) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 10)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_cru[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(j) RegCM5-CRU SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 11)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_cpc[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(k) RegCM5-CPC SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 3, 12)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_era5[0], levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(l) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')

else:
	fig = plt.figure(figsize=(6, 8))

	ax = fig.add_subplot(4, 2, 1)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_cru, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(a) RegCM5-CRU DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 2)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_djf_regcm_era5, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(b) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 3)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_cru, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(c) RegCM5-CRU MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 4)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_mam_regcm_era5, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(d) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 5)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_cru, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(e) RegCM5-CRU JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 6)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_jja_regcm_era5, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(f) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 7)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_cru, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(g) RegCM5-CRU SON', loc='left', fontsize=font_size, fontweight='bold')

	ax = fig.add_subplot(4, 2, 8)  
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, mbe_son_regcm_era5, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='both') 
	plt.title(u'(h) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')

# Set colobar
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.94, 0.3, 0.015, 0.4]))
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/user/mdasilva/sam_3km/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_SAM-3km_RegCM5_2018-2021.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
