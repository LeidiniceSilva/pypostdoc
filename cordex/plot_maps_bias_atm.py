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

from mpl_toolkits.basemap import Basemap

var = 'quv'
domain = 'CSAM-3'
level = '200'
path = '/marconi/home/userexternal/mdasilva'

idt, fdt = '2000', '2009'
dt = '{0}-{1}'.format(idt, fdt)


def import_obs(param, level, domain, dataset, season):

	if param == 'uwnd':
		param_ = 'u'
	elif param == 'vwnd':
                param_ = 'v'
	else:
		param_ = 'q'

	arq   = '{0}/user/mdasilva/CORDEX/post_evaluate/obs/{1}_{2}hPa_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, param, level, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param_][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,0,:,:]
	
	return lat, lon, mean


def import_rcm(param, level, domain, dataset, season):

	if level =='200':
		param_ = '{0}200'.format(param)
	else:
		param_ = '{0}850'.format(param)

	arq   = '{0}/user/mdasilva/CORDEX/post_evaluate/rcm/{1}{2}_{3}_{4}_{5}_2000-2005_lonlat.nc'.format(path, param, level, domain, dataset, season)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param_][:]
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

	return lat, lon, mean
	
	
def basemap(lat, lon):

        map = Basemap(projection='cyl', llcrnrlon=-80., llcrnrlat=-38., urcrnrlon=-34.,urcrnrlat=-10., resolution='c')
        map.drawmeridians(np.arange(-80., -34., 12.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
        map.drawparallels(np.arange(-38., -8., 6.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
        map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black')

        lons, lats = np.meshgrid(lon, lat)
        xx, yy = map(lons,lats)

        return map, xx, yy


# Import model and obs dataset 
lat, lon, u_era5_djf = import_obs('uwnd', level, domain, 'ERA5', 'DJF')
lat, lon, u_era5_mam = import_obs('uwnd', level, domain, 'ERA5', 'MAM')
lat, lon, u_era5_jja = import_obs('uwnd', level, domain, 'ERA5', 'JJA')
lat, lon, u_era5_son = import_obs('uwnd', level, domain, 'ERA5', 'SON')
lat, lon, u_regcm_djf = import_rcm('ua', level, domain, 'RegCM5', 'DJF')
lat, lon, u_regcm_mam = import_rcm('ua', level, domain, 'RegCM5', 'MAM')
lat, lon, u_regcm_jja = import_rcm('ua', level, domain, 'RegCM5', 'JJA')
lat, lon, u_regcm_son = import_rcm('ua', level, domain, 'RegCM5', 'SON')
mbe_djf_regcm_era5_u = u_regcm_djf - u_era5_djf
mbe_mam_regcm_era5_u = u_regcm_mam - u_era5_mam
mbe_jja_regcm_era5_u = u_regcm_jja - u_era5_jja
mbe_son_regcm_era5_u = u_regcm_son - u_era5_son	

lat, lon, v_era5_djf = import_obs('vwnd', level, domain, 'ERA5', 'DJF')
lat, lon, v_era5_mam = import_obs('vwnd', level, domain, 'ERA5', 'MAM')
lat, lon, v_era5_jja = import_obs('vwnd', level, domain, 'ERA5', 'JJA')
lat, lon, v_era5_son = import_obs('vwnd', level, domain, 'ERA5', 'SON')
lat, lon, v_regcm_djf = import_rcm('va', level, domain, 'RegCM5', 'DJF')
lat, lon, v_regcm_mam = import_rcm('va', level, domain, 'RegCM5', 'MAM')
lat, lon, v_regcm_jja = import_rcm('va', level, domain, 'RegCM5', 'JJA')
lat, lon, v_regcm_son = import_rcm('va', level, domain, 'RegCM5', 'SON')
mbe_djf_regcm_era5_v = v_regcm_djf - v_era5_djf
mbe_mam_regcm_era5_v = v_regcm_mam - v_era5_mam
mbe_jja_regcm_era5_v = v_regcm_jja - v_era5_jja
mbe_son_regcm_era5_v = v_regcm_son - v_era5_son	

lat, lon, q_era5_djf = import_obs('qhum', level, domain, 'ERA5', 'DJF')
lat, lon, q_era5_mam = import_obs('qhum', level, domain, 'ERA5', 'MAM')
lat, lon, q_era5_jja = import_obs('qhum', level, domain, 'ERA5', 'JJA')
lat, lon, q_era5_son = import_obs('qhum', level, domain, 'ERA5', 'SON')
lat, lon, q_regcm_djf = import_rcm('hus', level, domain, 'RegCM5', 'DJF')
lat, lon, q_regcm_mam = import_rcm('hus', level, domain, 'RegCM5', 'MAM')
lat, lon, q_regcm_jja = import_rcm('hus', level, domain, 'RegCM5', 'JJA')
lat, lon, q_regcm_son = import_rcm('hus', level, domain, 'RegCM5', 'SON')
mbe_djf_regcm_era5_q = q_regcm_djf - q_era5_djf
mbe_mam_regcm_era5_q = q_regcm_mam - q_era5_mam
mbe_jja_regcm_era5_q = q_regcm_jja - q_era5_jja
mbe_son_regcm_era5_q = q_regcm_son - q_era5_son

# Plot figure 
fig = plt.figure(figsize=(8, 8))
font_size = 8
mulc = 1000

if level == '200':
	dict_plot = {
	'u': [u'Bias of zonal wind {0} (m s$^-$$^1$)'.format(level), np.arange(-10, 11, 1), cm.RdBu_r],
	'v': [u'Bias of meridional wind {0} (m s$^-$$^1$)'.format(level), np.arange(-10, 11, 1), cm.RdBu_r],
	'q': [u'Bias of specific humidity {0} (10e$^-$$^3$kg kg$^-$$^1$)'.format(level), np.arange(-0.07, 0.09, 0.02), cm.BrBG]
	}
else:
	dict_plot = {
	'u': [u'Bias of zonal wind {0} (m s$^-$$^1$)'.format(level), np.arange(-10, 11, 1), cm.RdBu_r],
	'v': [u'Bias of meridional wind {0} (m s$^-$$^1$)'.format(level), np.arange(-10, 11, 1), cm.RdBu_r],
	'q': [u'Bias of specific humidity {0} (10e$^-$$^3$kg kg$^-$$^1$)'.format(level), np.arange(-8, 9, 1), cm.BrBG]
	}
	
ax = fig.add_subplot(4, 3, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_djf_regcm_era5_u, levels=dict_plot['u'][1], cmap=dict_plot['u'][2], extend='neither') 
plt.title(u'(a) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

# Set colobar
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.3, 0.015, 0.4]))
cbar.set_label('{0}'.format(dict_plot['u'][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

ax = fig.add_subplot(4, 3, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_djf_regcm_era5_v, levels=dict_plot['v'][1], cmap=dict_plot['v'][2], extend='neither') 
plt.title(u'(b) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

# Set colobar
cbar = plt.colorbar(plt_map, cax=fig.add_axes([1.03, 0.3, 0.015, 0.4]))
cbar.set_label('{0}'.format(dict_plot['v'][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

ax = fig.add_subplot(4, 3, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_djf_regcm_era5_q*mulc, levels=dict_plot['q'][1], cmap=dict_plot['q'][2], extend='neither') 
plt.title(u'(c) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

# Set colobar
cbar = plt.colorbar(plt_map, cax=fig.add_axes([1.13, 0.3, 0.015, 0.4]))
cbar.set_label('{0}'.format(dict_plot['q'][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

ax = fig.add_subplot(4, 3, 4)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mam_regcm_era5_u, levels=dict_plot['u'][1], cmap=dict_plot['u'][2], extend='neither') 
plt.title(u'(d) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 5)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mam_regcm_era5_v, levels=dict_plot['v'][1], cmap=dict_plot['v'][2], extend='neither') 
plt.title(u'(e) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 6)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mam_regcm_era5_q*mulc, levels=dict_plot['q'][1], cmap=dict_plot['q'][2], extend='neither') 
plt.title(u'(f) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 7)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_jja_regcm_era5_u, levels=dict_plot['u'][1], cmap=dict_plot['u'][2], extend='neither') 
plt.title(u'(g) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 8)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_jja_regcm_era5_v, levels=dict_plot['v'][1], cmap=dict_plot['v'][2], extend='neither') 
plt.title(u'(h) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 9)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_jja_regcm_era5_q*mulc, levels=dict_plot['q'][1], cmap=dict_plot['q'][2], extend='neither') 
plt.title(u'(i) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 10)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_son_regcm_era5_u, levels=dict_plot['u'][1], cmap=dict_plot['u'][2], extend='neither') 
plt.title(u'(j) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 11)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_son_regcm_era5_v, levels=dict_plot['v'][1], cmap=dict_plot['v'][2], extend='neither') 
plt.title(u'(k) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 12)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_son_regcm_era5_q*mulc, levels=dict_plot['q'][1], cmap=dict_plot['q'][2], extend='neither') 
plt.title(u'(l) CPM3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
	
# Path out to save figure
path_out = '{0}/user/mdasilva/CORDEX/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}hPa_CSAM-3_RegCM5_{2}.png'.format(var, level, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
