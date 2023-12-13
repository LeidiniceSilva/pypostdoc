# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot bias maps"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from mpl_toolkits.basemap import Basemap
from import_dataset_situ import import_obs_situ, import_rcm_situ
from import_dataset_grid import import_obs_srf, import_rcm_srf

path='/marconi/home/userexternal/mdasilva'

domain = 'SAM-3km'

var = 'pr'

dict_var = {
'pr': ['pre', 'pre', 'precip', 'sat_gauge_precip', 'tp'],
'tas': ['tmp', 'tmp', 't2m'],
'tasmax': ['tmx', 'tmax', 'mx2t'],
'tasmin': ['tmn', 'tmin', 'mn2t'],
'clt': ['cld', 'tmp', 't2m']
}
	

def basemap(lat, lon):
	
	map = Basemap(projection='cyl', llcrnrlon=-80., llcrnrlat=-38., urcrnrlon=-34.,urcrnrlat=-8., resolution='c')
	map.drawmeridians(np.arange(-80., -34., 12.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-38., -8., 6.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
		
	lons, lats = np.meshgrid(lon, lat)
	xx, yy = map(lons,lats)

	# Import shapefile 	
	map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black')
	
	return map, xx, yy


# Import model and obs database 
lat_regcm, lon_regcm, regcm = import_rcm_situ(var)
lat_inmet, lon_inmet, inmet = import_obs_situ(dict_var[var][0])

lat, lon, reg = import_rcm_srf(var, domain, 'RegCM5')
lat_cru, lon_cru, cru = import_obs_srf(dict_var[var][1], domain, 'CRU')
lat_cpc, lon_cpc, cpc = import_obs_srf(dict_var[var][2], domain, 'CPC')
lat_gpcp, lon_gpcp, gpcp = import_obs_srf(dict_var[var][3], domain, 'GPCP')
lat_era5, lon_era5, era5 = import_obs_srf(dict_var[var][4], domain, 'ERA5')

bias_djf_reg_inmet, bias_mam_reg_inmet, bias_jja_reg_inmet, bias_son_reg_inmet = [], [], [], []

# Calculate seasonal regcm bias
for i in range(0, 298):
	
	bias_djf_reg_inmet.append(regcm[i][0] - inmet[i][0])
	bias_mam_reg_inmet.append(regcm[i][1] - inmet[i][1])
	bias_jja_reg_inmet.append(regcm[i][2] - inmet[i][2])
	bias_son_reg_inmet.append(regcm[i][3] - inmet[i][3])

bias_djf_reg_cru = reg[0] - cru[0]
bias_mam_reg_cru = reg[1] - cru[1]
bias_jja_reg_cru = reg[2] - cru[2]
bias_son_reg_cru = reg[3] - cru[3]

bias_djf_reg_cpc = reg[0] - cpc[0]
bias_mam_reg_cpc = reg[1] - cpc[1]
bias_jja_reg_cpc = reg[2] - cpc[2]
bias_son_reg_cpc = reg[3] - cpc[3]

bias_djf_reg_gpcp = reg[0] - gpcp[0]
bias_mam_reg_gpcp = reg[1] - gpcp[1]
bias_jja_reg_gpcp = reg[2] - gpcp[2]
bias_son_reg_gpcp = reg[3] - gpcp[3]

bias_djf_reg_era5 = reg[0] - era5[0]
bias_mam_reg_era5 = reg[1] - era5[1]
bias_jja_reg_era5 = reg[2] - era5[2]
bias_son_reg_era5 = reg[3] - era5[3]

# Plot figure   
fig = plt.figure(figsize=(10, 6))
font_size = 8

if var == 'pr':
	v_min = -10
	v_max = 10
	legend = 'Bias of  precipitation (mm d⁻¹)'
	levs = [-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
	color = cm.BrBG
elif var == 'tas':
	legend = 'Bias of air temperature (°C)'
	levs = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]
	color = cm.bwr
elif var == 'tasmax':
	legend = 'Bias of maximum air temperature (°C)'
	levs = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]
	color = cm.bwr
elif var == 'tasmin':
	legend = 'Bias of minimum air temperature (°C)'
	levs1 = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]
	color1 = cm.bwr	
else:
	legend = 'Bias of total cloud cover (%)'
	levs = [-50, -40, -30, -20, -10, -5, 0, 5, 10, 20, 30, 40, 50]
	color = cm.RdGy
	
ax = fig.add_subplot(4, 5, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.scatter(lon_regcm, lat_regcm, 4, bias_djf_reg_inmet, cmap=color, marker='o', vmin=v_min, vmax=v_max) 
plt.title(u'(a) RegCM5-INMET (DJF)', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_djf_reg_cru, levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(b) RegCM5-CRU (DJF)', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_djf_reg_cpc, levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(c) RegCM5-CPC (DJF)', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 4)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_djf_reg_gpcp, levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(d) RegCM5-GPCP (DJF)', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 5)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_djf_reg_era5, levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(e) RegCM5-ERA5 (DJF)', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 6)  
map, xx, yy = basemap(lat, lon)
plt_map = map.scatter(lon_regcm, lat_regcm, 4, bias_mam_reg_inmet, cmap=color, marker='o', vmin=v_min, vmax=v_max) 
plt.title(u'(f) RegCM5-INMET (MAM)', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 7)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_mam_reg_cru, levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(g) RegCM5-CRU (MAM)', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 8)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_mam_reg_cpc, levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(h) RegCM5-CPC (MAM)', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 9)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_mam_reg_gpcp, levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(i) RegCM5-GPCP (MAM)', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 10)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_mam_reg_era5, levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(j) RegCM5-ERA5 (MAM)', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 11)  
map, xx, yy = basemap(lat, lon)
plt_map = map.scatter(lon_regcm, lat_regcm, 4, bias_jja_reg_inmet, cmap=color, marker='o', vmin=v_min, vmax=v_max) 
plt.title(u'(k) RegCM5-INMET (JJA)', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 12)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_jja_reg_cru, levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(l) RegCM5-CRU (JJA)', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 13)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_jja_reg_cpc, levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(m) RegCM5-CPC (JJA)', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 14)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_jja_reg_gpcp, levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(n) RegCM5-GPCP (JJA)', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 15)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_jja_reg_era5, levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(o) RegCM5-ERA5 (JJA)', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 16)  
map, xx, yy = basemap(lat, lon)
plt_map = map.scatter(lon_regcm, lat_regcm, 4, bias_son_reg_inmet, cmap=color, marker='o', vmin=v_min, vmax=v_max) 
plt.title(u'(p) RegCM5-INMET (SON)', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=10, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 17)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_son_reg_cru, levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(q) RegCM5-CRU (SON)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=10, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 18)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_son_reg_cpc, levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(r) RegCM5-CPC (SON)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=10, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 19)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_son_reg_gpcp, levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(s) RegCM5-GPCP (SON)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=10, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 5, 20)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_son_reg_era5, levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(t) RegCM5-ERA5 (SON)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=10, fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.94, 0.3, 0.015, 0.4]))
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/user/mdasilva/sam_3km/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_SAM-3km_RegCM5_2018-2021.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
