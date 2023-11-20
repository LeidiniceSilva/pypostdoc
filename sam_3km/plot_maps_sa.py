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

	arq = xr.open_dataset('/afs/ictp.it/home/m/mda_silv/Downloads/' + '{0}_SAM-3km_STS.2018010100_lonlat.nc'.format(param))
	data = arq[param]
	time = data.sel(time=slice('2018-01-01','2018-01-31'))
	var = time.groupby('time.year').mean('time')
	lat = var.lat
	lon = var.lon
	mean = var.values*86400
	
	return lat, lon, mean
	
	
def basemap(lat, lon):
	
	map = Basemap(projection='cyl', llcrnrlon=-90., llcrnrlat=-58., urcrnrlon=-30.,urcrnrlat=18., resolution='c')
	map.drawmeridians(np.arange(-90., -30., 10.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-58., 18., 10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')
		
	lons, lats = np.meshgrid(lon, lat)
	xx, yy = map(lons,lats)

	# Import shapefile 	
	path = '/afs/ictp.it/home/m/mda_silv/Documents/github_projects/shp'
	map.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black')
	
	return map, xx, yy


# Import model and obs database 
var_rcm = 'pr'

lat, lon, clim_rcm = import_rcm(var_rcm)

# Plot figure   
fig = plt.figure()
font_size = 8

if var_rcm == 'pr':
	legend = 'Precipitation (mm d⁻¹)'
	levs = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, 24, 26]
	color = cm.jet
else:
	legend = 'Temperature (C)'
	levs0 = [14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34]
	color0 = cm.Reds

ax = fig.add_subplot(1, 1, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, clim_rcm[0], levels=levs, latlon=True, cmap=color, extend='max') 
plt.title(u'(a) Jan-2018 RegCM5-CP-3km', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')
plt.text(-36, 11, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.78, 0.1, 0.02, 0.78]))
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '/afs/ictp.it/home/m/mda_silv/Documents/ICTP/figs/sam_3km'
name_out = 'pyplt_maps_{0}_sam_3km.png'.format(var_rcm)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


