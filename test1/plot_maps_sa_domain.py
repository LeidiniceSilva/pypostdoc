# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "September 16, 2023"
__description__ = "This script plot annual cycle"

import os
import netCDF4
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from mpl_toolkits.basemap import Basemap


def import_rcm(param):

	arq = xr.open_dataset('/home/nice/Downloads/' + '{0}_SAM-22_Reg5_mon_2000_lonlat.nc'.format(param))
	data = arq[param]
	time = data.sel(time=slice('2000-01-01','2000-12-31'))
	var = time.groupby('time.year').mean('time')
	lat = var.lat
	lon = var.lon
	mean = np.nanmean(var.values, axis=0)
	
	return lat, lon, mean
	
	
def basemap(lat, lon):
	
	aux_lon1 = []
	aux_lon2 = []
	for l in lon:
		if l <= 180:
			aux_lon1.append(l)
		else:
			aux_lon2.append(l-360)
		
	lon = np.array(aux_lon1[::-1] + aux_lon2[::-1])
	new_lat = lat
	new_lon = lon[::-1]

	map = Basemap(projection='cyl', llcrnrlon=-105., llcrnrlat=-57., urcrnrlon=-16.,urcrnrlat=18., resolution='c')
	map.drawmeridians(np.arange(-105., -16., 20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-57., 18., 15.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
		
	lons, lats = np.meshgrid(new_lon, new_lat)
	xx, yy = map(lons,lats)

	# Import shapefile 	
	path = '/home/nice/Documentos/github_projects/shp'
	map.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black')
	
	return map, xx, yy
	
	
# Import model and obs database 
var_rcm = 'tas'

lat, lon, clim_rcm = import_rcm(var_rcm)

# Plot figure   
fig = plt.figure()

levs = [-2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
color = cm.rainbow

ax = fig.add_subplot(1, 1, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, clim_rcm[0], levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'RegCM5 SA Domain (°C)', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=8, fontweight='bold')
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.90, 0.3, 0.02, 0.4]))
cbar.ax.tick_params(labelsize=8)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_domain.png'.format(var_rcm)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


