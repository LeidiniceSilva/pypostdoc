# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot bias maps"

import os
import netCDF4
import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from mpl_toolkits.basemap import Basemap
from import_climate_tools import compute_mbe

var = 'pr'
domain = 'SAM-3km'
idt, fdt = '2018', '2018'
dt = '{0}-{1}'.format(idt, fdt)

path = '/marconi/home/userexternal/mdasilva'

skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,7793,100,105,106,107,112,117,124,135,137,139,
149,152,155,158,168,174,177,183,186,199,204,210,212,224,226,239,240,248,249,253,254,276,277,280,293,298,
303,305,306,308,319,334,335,341,343,359,362,364,384,393,396,398,399,400,402,413,416,417,422,423,426,427,
443,444,446,451,453,457,458,467,474,479,483,488,489,490,495,505,509,513,514,516,529,534,544,559,566]
	
	
def import_inmet(param_i, param_ii, domain, dataset):
	
	yy, xx = [], []
	mean_i, mean_ii = [], []
	
	for station in range(1, 567):
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue

		yy.append(inmet[station][2])
		xx.append(inmet[station][3])

		arq_i  = xr.open_dataset('{0}/OBS/BDMET/database/nc/hourly/{1}/'.format(path, param_i) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(param_i, inmet[station][0]))
		data_i = arq_i[param_i]
		time_i = data_i.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		time_i = time_i.groupby('time.month').mean('time')
		var_i  = time_i.values*24
		mean_i.append(var_i)

		arq_ii  = xr.open_dataset('{0}/user/mdasilva/SAM-3km_v5/post/rcm/'.format(path) + '{0}_{1}_{2}_mon_{3}_lonlat.nc'.format(param_ii, domain, dataset, dt))
		data_ii = arq_ii[param_ii]
		data_ii = data_ii.sel(lat=slice(inmet[station][2]-0.03,inmet[station][2]+0.03),lon=slice(inmet[station][3]-0.03,inmet[station][3]+0.03)).mean(('lat','lon'))
		time_ii = data_ii.sel(time=slice('{0}-01-01'.format(idt),'{0}-12-31'.format(fdt)))
		var_ii  = time_ii.values
		mean_ii.append(var_ii)
		
	return yy, xx, mean_i, mean_ii
	
	
def basemap():
	
	map = Basemap(projection='cyl', llcrnrlon=-80., llcrnrlat=-38., urcrnrlon=-34.,urcrnrlat=-8., resolution='c')
	map.drawmeridians(np.arange(-80., -34., 12.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-38., -8., 6.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black')
	
	return map
	
	
# Import model and obs dataset
dict_var = {'pr': ['pre', 'pre', 'precip', 'sat_gauge_precip', 'tp']}

lat_i, lon_i, inmet_i, regcm_i = import_inmet(dict_var[var][0], var, domain, 'RegCM5')

mbe_jan, mbe_feb, mbe_mar, mbe_apr, mbe_may, mbe_jun, mbe_jul, mbe_aug, mbe_sep = [], [], [], [], [], [], [], [], []
for i in range(0, 298):
	mbe_jan.append(compute_mbe(regcm_i[i][0], inmet_i[i][0]))
	mbe_feb.append(compute_mbe(regcm_i[i][1], inmet_i[i][1]))
	mbe_mar.append(compute_mbe(regcm_i[i][2], inmet_i[i][2]))
	mbe_apr.append(compute_mbe(regcm_i[i][3], inmet_i[i][3]))
	mbe_may.append(compute_mbe(regcm_i[i][4], inmet_i[i][4]))
	mbe_jun.append(compute_mbe(regcm_i[i][5], inmet_i[i][5]))
	mbe_jul.append(compute_mbe(regcm_i[i][6], inmet_i[i][6]))
	mbe_aug.append(compute_mbe(regcm_i[i][7], inmet_i[i][7]))
	mbe_sep.append(compute_mbe(regcm_i[i][8], inmet_i[i][8]))

# Plot figure  
fig = plt.figure(figsize=(10, 7))

font_size = 8
cmap = plt.cm.BrBG
norm = mpl.colors.BoundaryNorm(np.arange(-10, 11, 1), cmap.N)
dict_plot = {'pr': ['Bias of  precipitation (mm d$^-$$^1$)']}

ax = fig.add_subplot(3, 3, 1)  
map = basemap()
plt_map = map.scatter(lon_i, lat_i, 4, mbe_jan, marker='o', cmap=cmap, norm=norm) 
plt.title(u'(a) RegCM5-INMET Jan', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 2)  
map = basemap()
plt_map = map.scatter(lon_i, lat_i, 4, mbe_feb, marker='o', cmap=cmap, norm=norm) 
plt.title(u'(b) RegCM5-INMET Feb', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 3)  
map = basemap()
plt_map = map.scatter(lon_i, lat_i, 4, mbe_mar, marker='o', cmap=cmap, norm=norm) 
plt.title(u'(c) RegCM5-INMET Mar', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 4)  
map = basemap()
plt_map = map.scatter(lon_i, lat_i, 4, mbe_apr, marker='o', cmap=cmap, norm=norm) 
plt.title(u'(d) RegCM5-INMET Apr', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 5)  
map = basemap()
plt_map = map.scatter(lon_i, lat_i, 4, mbe_may, marker='o', cmap=cmap, norm=norm) 
plt.title(u'(e) RegCM5-INMET May', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 6)  
map = basemap()
plt_map = map.scatter(lon_i, lat_i, 4, mbe_jun, marker='o', cmap=cmap, norm=norm) 
plt.title(u'(f) RegCM5-INMET Jun', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 7)  
map = basemap()
plt_map = map.scatter(lon_i, lat_i, 4, mbe_jul, marker='o', cmap=cmap, norm=norm) 
plt.title(u'(g) RegCM5-INMET Jul', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 8)  
map = basemap()
plt_map = map.scatter(lon_i, lat_i, 4, mbe_aug, marker='o', cmap=cmap, norm=norm) 
plt.title(u'(h) RegCM5-INMET Aug', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 9)  
map = basemap()
plt_map = map.scatter(lon_i, lat_i, 4, mbe_sep, marker='o', cmap=cmap, norm=norm) 
plt.title(u'(i) RegCM5-INMET Sep', loc='left', fontsize=font_size, fontweight='bold')

cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.3, 0.018, 0.4]))
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km_v5/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}_ws.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
