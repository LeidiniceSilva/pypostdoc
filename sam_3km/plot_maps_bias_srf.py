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
from import_dataset_situ import import_obs_situ

path='/marconi/home/userexternal/mdasilva/github_projects'


def basemap(lat, lon):
	
	map = Basemap(projection='cyl', llcrnrlon=-80., llcrnrlat=-38., urcrnrlon=-34.,urcrnrlat=-8., resolution='c')
	map.drawmeridians(np.arange(-80., -34., 6.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-38., -8., 6.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')
		
	lons, lats = np.meshgrid(lon, lat)
	xx, yy = map(lons,lats)

	# Import shapefile 	
	map.readshapefile('{0}/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black')
	
	return map, xx, yy


# Import model and obs database 
var_rcm = 'pr'

if var_rcm == 'pr':
	var_obs = 'pre'
elif var_rcm == 'tas':
	var_obs = 'tmp'
elif var_rcm == 'tasmax':
	var_obs = 'tmx'
elif var_rcm == 'tasmin':
	var_obs = 'tmn'
elif var_rcm == 'clt':
	var_obs = 'cld'
else:
	var_ref = 'v'
	
lat_inmet, lon_inmet, situ_inmet = import_obs_situ(var_obs)

print(lat_inmet)
print()
print(lon_inmet)
print()
print(situ_inmet)

exit()

obs_situ_reg = import_rcm_situ(var_obs)





bias_reg_inmet = obs_situ_reg - obs_situ_inmet
print(bias_reg_inmet)
exit()


# Plot figure   
fig = plt.figure(figsize=(10, 8))
font_size = 8

if var_rcm == 'pr':
	legend = 'Bias of  precipitation (mm d⁻¹)'
	levs0 = [0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12]
	color0 = cm.Blues
	levs1 = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]
	color1 = cm.BrBG
elif var_rcm == 'clt':
	legend = 'Bias of total cloud cover (%)'
	levs0 = [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
	color0 = cm.Greys
	levs1 = [-50, -40, -30, -20, -10, -5, 0, 5, 10, 20, 30, 40, 50]
	color1 = cm.RdGy
elif var_rcm == 'tas':
	legend = 'Bias of air temperature (°C)'
	levs0 = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34]
	color0 = cm.Reds
	levs1 = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]
	color1 = cm.bwr
elif var_rcm == 'tasmax':
	legend = 'Bias of maximum air temperature (°C)'
	levs0 = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34]
	color0 = cm.Reds
	levs1 = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]
	color1 = cm.bwr
else:
	legend = 'Bias of minimum air temperature (°C)'
	levs0 = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34]
	color0 = cm.Reds
	levs1 = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]
	color1 = cm.bwr

ax = fig.add_subplot(1, 3, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, clim_ref, levels=levs0, latlon=True, cmap=color0, extend='max') 
plt.title(u'(a) CRU', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=25, fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(1, 3, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, clim_rcm, levels=levs0, latlon=True, cmap=color0, extend='max') 
plt.title(u'(b) Reg-3km', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.91, 0.3, 0.015, 0.4]))
cbar.ax.tick_params(labelsize=font_size)

ax = fig.add_subplot(1, 3, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_rcm_ref, levels=levs1, latlon=True, cmap=color1, extend='both') 
plt.title(u'(c) Reg-3km - CRU', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.96, 0.3, 0.015, 0.4]))
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/ICTP/figs/sam_3km'.format(path)
name_out = 'pyplt_maps_{0}_SAM-3km_RegCM5_cru_ts4.07_2018-2021.png'.format(var_rcm)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
