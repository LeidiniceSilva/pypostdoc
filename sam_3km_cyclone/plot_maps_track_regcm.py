# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jan 02, 2024"
__description__ = "This script plot cyclone tracking"

import os
import cmocean
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import maskoceans

var = 'psl'
cyclone = 'cyclone_ii'
path = '/marconi/home/userexternal/mdasilva'

if cyclone == 'cyclone_i':
	t_i = 103
	t_f = 175
	t1 = pd.to_datetime('2023-06-14 00:00:00')
	t2 = pd.to_datetime('2023-06-22 21:00:00')
else:
	t_i = 327
	t_f = 367
	t1 = pd.to_datetime('2023-07-12 00:00:00')
	t2 = pd.to_datetime('2023-07-16 21:00:00')
	
	
def import_grid(param, domain, exp, dataset, freq, dt):

	arq   = '{0}/user/mdasilva/SAM-3km-cyclone/post/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, param, domain, exp, dataset, freq, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][t_i:t_f,:,:]

	return lat, lon, mean


def basemap(lat, lon):

	map = Basemap(projection='cyl', llcrnrlon=-82., llcrnrlat=-50., urcrnrlon=-32.,urcrnrlat=-10., resolution='c')
	map.drawmeridians(np.arange(-82., -32., 5.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-50., -10., 5.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black')
	
	lons, lats = np.meshgrid(lon, lat)
	xx, yy = map(lons,lats)
	
	return map, xx, yy


# Import model and obs dataset
lat, lon, regcm = import_grid(var, 'SAM-3km-cyclone', 'ECMWF-ERA5_evaluation_r1i1p1f1', 'ICTP-RegCM5', '3h', '2023060100-2023083100')

series = pd.date_range(t1,t2,freq='180min')
iso8601 = [t.strftime('%Y%m%dT%H:%M%SZ') for t in series]

for i in range(0, regcm.shape[0]):
	print(iso8601[i])

	# Plot figure
	fig = plt.figure()
	font_size = 8
	
	if var == 'pr':
		ax = fig.add_subplot(1, 1, 1)  
		map, xx, yy = basemap(lat, lon)
		plt_map = map.contourf(xx, yy, regcm[i], cmap=matplotlib.colors.ListedColormap(["#ffffffff","#d7f0fcff","#ade0f7ff","#86c4ebff","#60a5d6ff","#4794b3ff","#49a67cff","#55b848ff","#9ecf51ff","#ebe359ff","#f7be4aff","#f58433ff","#ed5a28ff","#de3728ff","#cc1f27ff","#b01a1fff","#911419ff"])) 
		plt.title(u'CP-RegCM5 Cyclone I {0}'.format(iso8601[i]), loc='left', fontsize=font_size, fontweight='bold')
	else:
		ax = fig.add_subplot(1, 1, 1)  
		map, xx, yy = basemap(lat, lon)
		plt_map = map.contour(xx, yy, regcm[i]/100, cmap=cm.jet, linewidths=0.5) 
		plt.clabel(plt_map, inline=1, fontsize=font_size)
		plt.title(u'CP-RegCM5 Cyclone II {0}'.format(iso8601[i]), loc='left', fontsize=font_size, fontweight='bold')
		
	# Path out to save figure
	path_out = '{0}/user/mdasilva/SAM-3km-cyclone/figs/{1}/{2}'.format(path, var, cyclone)
	name_out = 'pyplt_maps_tracking_cyclone_cp-regcm5_{0}_{1}.png'.format(var, iso8601[i])
	plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
	plt.close('all')
	plt.cla()
exit()
