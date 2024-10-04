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

var = 'uv'
domain = 'CSAM-3'
level = '200'
path='/marconi/home/userexternal/mdasilva'

idt, fdt = '2000', '2009'
dt = '{0}-{1}'.format(idt, fdt)


def import_obs(param, level, domain, dataset, season):

        if param == 'uwnd':
                param_ = 'u'
        elif param == 'vwnd':
                param_ = 'v'
        else:
                param_ = 'q'

	arq   = '{0}/user/mdasilva/CORDEX/post_evaluate/obs/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, param, level, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,0,:,:]
	
	return lat, lon, mean
	

def import_rcm(param, level, domain, dataset, season):

        if level =='200':
                param_ = '{0}200'.format(param)
        else:
                param_ = '{0}850'.format(param)

	arq   = '{0}/user/mdasilva/SAM-3km/post_evaluate/rcm/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, param, level, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]
	
	return lat, lon, mean


def basemap(lat, lon):
	
	map = Basemap(projection='cyl', llcrnrlon=-80., llcrnrlat=-38., urcrnrlon=-34.,urcrnrlat=-8., resolution='c')
	map.drawmeridians(np.arange(-80., -34., 12.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-38., -8., 6.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black')
	
	lons, lats = np.meshgrid(lon, lat)
	xx, yy = map(lons,lats)
	
	return map, xx, yy
	
	
# Import model and obs dataset 
lat, lon, u_era5_djf = import_obs('u', level, domain, 'ERA5', 'DJF')
lat, lon, u_era5_mam = import_obs('u', level, domain, 'ERA5', 'MAM')
lat, lon, u_era5_jja = import_obs('u', level, domain, 'ERA5', 'JJA')
lat, lon, u_era5_son = import_obs('u', level, domain, 'ERA5', 'SON')

lat, lon, v_era5_djf = import_obs('v', level, domain, 'ERA5', 'DJF')
lat, lon, v_era5_mam = import_obs('v', level, domain, 'ERA5', 'MAM')
lat, lon, v_era5_jja = import_obs('v', level, domain, 'ERA5', 'JJA')
lat, lon, v_era5_son = import_obs('v', level, domain, 'ERA5', 'SON')

lat, lon, u_regcm_djf = import_rcm('ua', level, domain, 'RegCM5', 'DJF')
lat, lon, u_regcm_mam = import_rcm('ua', level, domain, 'RegCM5', 'MAM')
lat, lon, u_regcm_jja = import_rcm('ua', level, domain, 'RegCM5', 'JJA')
lat, lon, u_regcm_son = import_rcm('ua', level, domain, 'RegCM5', 'SON')

lat, lon, v_regcm_djf = import_rcm('va', level, domain, 'RegCM5', 'DJF')
lat, lon, v_regcm_mam = import_rcm('va', level, domain, 'RegCM5', 'MAM')
lat, lon, v_regcm_jja = import_rcm('va', level, domain, 'RegCM5', 'JJA')
lat, lon, v_regcm_son = import_rcm('va', level, domain, 'RegCM5', 'SON')

uv_era5_djf = np.sqrt(u_era5_djf**2+v_era5_djf**2)
uv_era5_mam = np.sqrt(u_era5_mam**2+v_era5_mam**2)
uv_era5_jja = np.sqrt(u_era5_jja**2+v_era5_jja**2)
uv_era5_son = np.sqrt(u_era5_son**2+v_era5_son**2)

uv_regcm_djf = np.sqrt(u_regcm_djf**2+v_regcm_djf**2)
uv_regcm_mam = np.sqrt(u_regcm_mam**2+v_regcm_mam**2)
uv_regcm_jja = np.sqrt(u_regcm_jja**2+v_regcm_jja**2)
uv_regcm_son = np.sqrt(u_regcm_son**2+v_regcm_son**2)

mbe_u_djf = u_regcm_djf - u_era5_djf
mbe_u_mam = u_regcm_mam - u_era5_mam
mbe_u_jja = u_regcm_jja - u_era5_jja
mbe_u_son = u_regcm_son - u_era5_son

mbe_v_djf = v_regcm_djf - v_era5_djf
mbe_v_mam = v_regcm_mam - v_era5_mam
mbe_v_jja = v_regcm_jja - v_era5_jja
mbe_v_son = v_regcm_son - v_era5_son

mbe_uv_djf = np.sqrt(mbe_u_djf**2+mbe_v_djf**2)
mbe_uv_mam = np.sqrt(mbe_u_mam**2+mbe_v_mam**2)
mbe_uv_jja = np.sqrt(mbe_u_jja**2+mbe_v_jja**2)
mbe_uv_son = np.sqrt(mbe_u_son**2+mbe_v_son**2)

# Plot figure 
fig = plt.figure(figsize=(8, 8))
color = plt.cm.jet
font_size = 8
mag = 50

if level == '200hPa':
	vmin, vmax = 0, 85
	bmin, bmax = 0, 45
elif level == '500hPa':
	vmin, vmax = 0, 45
	bmin, bmax = 0, 25
else:
	vmin, vmax = 0, 16
	bmin, bmax = 0, 8	

ax = fig.add_subplot(4, 3, 1)
map, xx, yy = basemap(lat, lon)
plt_map = plt.quiver(xx[::mag, ::mag], yy[::mag, ::mag], u_era5_djf[::mag, ::mag], v_era5_djf[::mag, ::mag], uv_era5_djf[::mag, ::mag], cmap=color)
plt.title(u'(a) ERA5 {0} Wind DJF'.format(level), loc='left', fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(plt_map, cmap=plt.cm.jet, cax=fig.add_axes([0.92, 0.3, 0.015, 0.4]))
cbar.set_label('Wind speed (m s$^-$$^1$)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
plt.clim(vmin, vmax)

ax = fig.add_subplot(4, 3, 2)
map, xx, yy = basemap(lat, lon)
plt_map = plt.quiver(xx[::mag, ::mag], yy[::mag, ::mag], u_regcm_djf[::mag, ::mag], v_regcm_djf[::mag, ::mag], uv_regcm_djf[::mag, ::mag], cmap=color) 
plt.title(u'(b) RegCM5 {0} Wind DJF'.format(level), loc='left', fontsize=font_size, fontweight='bold')
plt.clim(vmin, vmax)

ax = fig.add_subplot(4, 3, 3)
map, xx, yy = basemap(lat, lon)
plt_map = plt.quiver(xx[::mag, ::mag], yy[::mag, ::mag], mbe_u_djf[::mag, ::mag], mbe_v_djf[::mag, ::mag], mbe_uv_djf[::mag, ::mag], cmap=color) 
plt.title(u'(c) RegCM5-ERA5', loc='left', fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(plt_map, cmap=plt.cm.jet, cax=fig.add_axes([1.03, 0.3, 0.015, 0.4]))
cbar.set_label('Bias of wind speed (m s$^-$$^1$)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
plt.clim(bmin, bmax)

ax = fig.add_subplot(4, 3, 4)
map, xx, yy = basemap(lat, lon)
plt_map = plt.quiver(xx[::mag, ::mag], yy[::mag, ::mag], u_era5_mam[::mag, ::mag], v_era5_mam[::mag, ::mag], uv_era5_mam[::mag, ::mag], cmap=color) 
plt.title(u'(d) ERA5 {0} Wind MAM'.format(level), loc='left', fontsize=font_size, fontweight='bold')
plt.clim(vmin, vmax)

ax = fig.add_subplot(4, 3, 5)
map, xx, yy = basemap(lat, lon)
plt_map = plt.quiver(xx[::mag, ::mag], yy[::mag, ::mag], u_regcm_mam[::mag, ::mag], v_regcm_mam[::mag, ::mag], uv_regcm_mam[::mag, ::mag], cmap=color) 
plt.title(u'(e) RegCM5 {0} Wind MAM'.format(level), loc='left', fontsize=font_size, fontweight='bold')
plt.clim(vmin, vmax)

ax = fig.add_subplot(4, 3, 6)
map, xx, yy = basemap(lat, lon)
plt_map = plt.quiver(xx[::mag, ::mag], yy[::mag, ::mag], mbe_u_mam[::mag, ::mag], mbe_v_mam[::mag, ::mag], mbe_uv_mam[::mag, ::mag], cmap=color) 
plt.title(u'(f) RegCM5-ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.clim(bmin, bmax)

ax = fig.add_subplot(4, 3, 7)
map, xx, yy = basemap(lat, lon)
plt_map = plt.quiver(xx[::mag, ::mag], yy[::mag, ::mag], u_era5_jja[::mag, ::mag], v_era5_jja[::mag, ::mag], uv_era5_jja[::mag, ::mag], cmap=color) 
plt.title(u'(g) ERA5 {0} Wind JJA'.format(level), loc='left', fontsize=font_size, fontweight='bold')
plt.clim(vmin, vmax)

ax = fig.add_subplot(4, 3, 8)
map, xx, yy = basemap(lat, lon)
plt_map = plt.quiver(xx[::mag, ::mag], yy[::mag, ::mag], u_regcm_jja[::mag, ::mag], v_regcm_jja[::mag, ::mag], uv_regcm_jja[::mag, ::mag], cmap=color) 
plt.title(u'(h) RegCM5 {0} Wind JJA'.format(level), loc='left', fontsize=font_size, fontweight='bold')
plt.clim(vmin, vmax)

ax = fig.add_subplot(4, 3, 9)
map, xx, yy = basemap(lat, lon)
plt_map = plt.quiver(xx[::mag, ::mag], yy[::mag, ::mag], mbe_u_jja[::mag, ::mag], mbe_v_jja[::mag, ::mag], mbe_uv_jja[::mag, ::mag], cmap=color) 
plt.title(u'(i) RegCM5-ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.clim(bmin, bmax)

ax = fig.add_subplot(4, 3, 10)
map, xx, yy = basemap(lat, lon)
plt_map = plt.quiver(xx[::mag, ::mag], yy[::mag, ::mag], u_era5_son[::mag, ::mag], v_era5_son[::mag, ::mag], uv_era5_son[::mag, ::mag], cmap=color) 
plt.title(u'(j) ERA5 {0} Wind SON'.format(level), loc='left', fontsize=font_size, fontweight='bold')
plt.clim(vmin, vmax)

ax = fig.add_subplot(4, 3, 11)
map, xx, yy = basemap(lat, lon)
plt_map = plt.quiver(xx[::mag, ::mag], yy[::mag, ::mag], u_regcm_son[::mag, ::mag], v_regcm_son[::mag, ::mag], uv_regcm_son[::mag, ::mag], cmap=color) 
plt.title(u'(k) RegCM5 {0} Wind SON'.format(level), loc='left', fontsize=font_size, fontweight='bold')
plt.clim(vmin, vmax)

ax = fig.add_subplot(4, 3, 12)
map, xx, yy = basemap(lat, lon)
plt_map = plt.quiver(xx[::mag, ::mag], yy[::mag, ::mag], mbe_u_son[::mag, ::mag], mbe_v_son[::mag, ::mag], mbe_uv_son[::mag, ::mag], cmap=color) 
plt.title(u'(l) RegCM5-ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.clim(bmin, bmax)

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/evaluate'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_{2}_RegCM5_{2}.png'.format(var, level, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
