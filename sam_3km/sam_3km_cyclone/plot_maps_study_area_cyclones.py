# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot map of study area"


import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import metpy.calc as mpcalc
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap

font_size = 10
path='/marconi/home/userexternal/mdasilva'


def import_mean(param):

	arq   = '{0}/user/mdasilva/SAM-3km/post_cyclone/rcm/preproc/{1}.2018.00.nc'.format(path, param)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]	
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	

def basemap(lat, lon):
	
	map = Basemap(projection='cyl', llcrnrlon=-76., llcrnrlat=-34.5, urcrnrlon=-38.5,urcrnrlat=-15., resolution='c')
	map.drawmeridians(np.arange(-76., -38.5, 5.0), size=font_size, labels=[0,0,0,1], linewidth=0.01, color='gray')
	map.drawparallels(np.arange(-34.5, -15., 2.5), size=font_size, labels=[1,0,0,0], linewidth=0.01, color='gray')
	map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', linewidth=1., color='red')

	lons, lats = np.meshgrid(lon, lat)
	xx, yy = map(lons,lats)
	
	return map, xx, yy
	
	
# Import model and obs dataset 
lat, lon, psl_mean = import_mean('psl')
lat, lon, ua_mean = import_mean('ua')
lat, lon, va_mean = import_mean('va')

psl_t1 = psl_mean[0]
ua_t1 = ua_mean[0][0]
va_t1 = va_mean[0][0]

wind_speed = np.sqrt(ua_mean**2 + va_mean**2)
ws_t1 = wind_speed[0][0]

# Plot figure 
fig = plt.figure()

ax = fig.add_subplot(1, 1, 1) 
map, xx, yy = basemap(lat, lon)

plt_map1 = plt.contourf(xx, yy, ws_t1, np.arange(0, 18, 1), cmap='Blues', extend='max')
cbar = plt.colorbar(plt_map1, ax=ax, pad=0.02, aspect=16, shrink=0.8)

plt_map2 = map.contour(xx, yy, psl_t1/100, levels=np.arange(996, 1012, 2), colors='black', linewidths=0.75)
plt.clabel(plt_map2, inline=1, fontsize=font_size)

plt_map3 = map.quiver(xx[::1, ::1], yy[::1, ::1], ua_t1[::1, ::1], va_t1[::1, ::1], color='gray')
plt.quiverkey(plt_map3, X=0.9, Y=1.05, U=np.max(ua_t1), label='15 m/s', labelpos='E', fontproperties={'size': 8})

plt.title(u'RegCM5 2018-01-01 00UTC', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=font_size, labelpad=15, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=font_size, labelpad=40, fontweight='bold')

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cyclone'.format(path)
name_out = 'pyplt_study_area_cyclones_RegCM5_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
