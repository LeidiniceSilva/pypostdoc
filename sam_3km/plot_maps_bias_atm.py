# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot bias maps"

import os
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from import_dataset_grid import import_grid_atm
from import_climate_tools import compute_mbe
from import_basemap import basemap

path='/marconi/home/userexternal/mdasilva'

domain = 'SAM-3km'

var = 'uvq'

level = '200hPa'

# Import model and obs dataset 
lat, lon, u_era5_djf = import_grid_atm('u', level, domain, 'ERA5', 'DJF')
lat, lon, u_era5_mam = import_grid_atm('u', level, domain, 'ERA5', 'MAM')
lat, lon, u_era5_jja = import_grid_atm('u', level, domain, 'ERA5', 'JJA')
lat, lon, u_era5_son = import_grid_atm('u', level, domain, 'ERA5', 'SON')
lat, lon, u_regcm_djf = import_grid_atm('ua', level, domain, 'RegCM5', 'DJF')
lat, lon, u_regcm_mam = import_grid_atm('ua', level, domain, 'RegCM5', 'MAM')
lat, lon, u_regcm_jja = import_grid_atm('ua', level, domain, 'RegCM5', 'JJA')
lat, lon, u_regcm_son = import_grid_atm('ua', level, domain, 'RegCM5', 'SON')
mbe_djf_regcm_era5_u = compute_mbe(u_regcm_djf[0], u_era5_djf[0])
mbe_mam_regcm_era5_u = compute_mbe(u_regcm_mam[0], u_era5_mam[0])
mbe_jja_regcm_era5_u = compute_mbe(u_regcm_jja[0], u_era5_jja[0])
mbe_son_regcm_era5_u = compute_mbe(u_regcm_son[0], u_era5_son[0])	
	
lat, lon, v_era5_djf = import_grid_atm('v', level, domain, 'ERA5', 'DJF')
lat, lon, v_era5_mam = import_grid_atm('v', level, domain, 'ERA5', 'MAM')
lat, lon, v_era5_jja = import_grid_atm('v', level, domain, 'ERA5', 'JJA')
lat, lon, v_era5_son = import_grid_atm('v', level, domain, 'ERA5', 'SON')
lat, lon, v_regcm_djf = import_grid_atm('va', level, domain, 'RegCM5', 'DJF')
lat, lon, v_regcm_mam = import_grid_atm('va', level, domain, 'RegCM5', 'MAM')
lat, lon, v_regcm_jja = import_grid_atm('va', level, domain, 'RegCM5', 'JJA')
lat, lon, v_regcm_son = import_grid_atm('va', level, domain, 'RegCM5', 'SON')
mbe_djf_regcm_era5_v = compute_mbe(v_regcm_djf[0], v_era5_djf[0])
mbe_mam_regcm_era5_v = compute_mbe(v_regcm_mam[0], v_era5_mam[0])
mbe_jja_regcm_era5_v = compute_mbe(v_regcm_jja[0], v_era5_jja[0])
mbe_son_regcm_era5_v = compute_mbe(v_regcm_son[0], v_era5_son[0])	

lat, lon, q_era5_djf = import_grid_atm('q', level, domain, 'ERA5', 'DJF')
lat, lon, q_era5_mam = import_grid_atm('q', level, domain, 'ERA5', 'MAM')
lat, lon, q_era5_jja = import_grid_atm('q', level, domain, 'ERA5', 'JJA')
lat, lon, q_era5_son = import_grid_atm('q', level, domain, 'ERA5', 'SON')
lat, lon, q_regcm_djf = import_grid_atm('hus', level, domain, 'RegCM5', 'DJF')
lat, lon, q_regcm_mam = import_grid_atm('hus', level, domain, 'RegCM5', 'MAM')
lat, lon, q_regcm_jja = import_grid_atm('hus', level, domain, 'RegCM5', 'JJA')
lat, lon, q_regcm_son = import_grid_atm('hus', level, domain, 'RegCM5', 'SON')
mbe_djf_regcm_era5_q = compute_mbe(q_regcm_djf[0], q_era5_djf[0])
mbe_mam_regcm_era5_q = compute_mbe(q_regcm_mam[0], q_era5_mam[0])
mbe_jja_regcm_era5_q = compute_mbe(q_regcm_jja[0], q_era5_jja[0])
mbe_son_regcm_era5_q = compute_mbe(q_regcm_son[0], q_era5_son[0])


# Plot figure   
dict_plot = {
'u': [u'Bias of zonal wind {0} (m s⁻¹)'.format(level), np.arange(-20, 22, 2), cm.RdBu_r],
'v': [u'Bias of meridional wind {0} (m s⁻¹)'.format(level), np.arange(-20, 22, 2), cm.RdBu_r],
'q': [u'Bias of specific humidity {0} (kg kg⁻¹)'.format(level), np.arange(-0.02, 0.03, 0.005), cm.BrBG]
}

font_size = 8
fig = plt.figure(figsize=(8, 8))

ax = fig.add_subplot(4, 3, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_djf_regcm_era5_u[0], levels=dict_plot['u'][1], cmap=dict_plot['u'][2], extend='both') 
plt.title(u'(a) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')

# Set colobar
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.94, 0.3, 0.015, 0.4]))
cbar.set_label('{0}'.format(dict_plot['u'][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

ax = fig.add_subplot(4, 3, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_djf_regcm_era5_v[0], levels=dict_plot['v'][1], cmap=dict_plot['v'][2], extend='both') 
plt.title(u'(b) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')

# Set colobar
cbar = plt.colorbar(plt_map, cax=fig.add_axes([1.05, 0.3, 0.015, 0.4]))
cbar.set_label('{0}'.format(dict_plot['v'][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

ax = fig.add_subplot(4, 3, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_djf_regcm_era5_q[0], levels=dict_plot['q'][1], cmap=dict_plot['q'][2], extend='both') 
plt.title(u'(c) RegCM5-ERA5 DJF', loc='left', fontsize=font_size, fontweight='bold')

# Set colobar
cbar = plt.colorbar(plt_map, cax=fig.add_axes([1.15, 0.3, 0.015, 0.4]))
cbar.set_label('{0}'.format(dict_plot['q'][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

ax = fig.add_subplot(4, 3, 4)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mam_regcm_era5_u[0], levels=dict_plot['u'][1], cmap=dict_plot['u'][2], extend='both') 
plt.title(u'(d) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 5)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mam_regcm_era5_v[0], levels=dict_plot['v'][1], cmap=dict_plot['v'][2], extend='both') 
plt.title(u'(e) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 6)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_mam_regcm_era5_q[0], levels=dict_plot['q'][1], cmap=dict_plot['q'][2], extend='both') 
plt.title(u'(f) RegCM5-ERA5 MAM', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 7)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_jja_regcm_era5_u[0], levels=dict_plot['u'][1], cmap=dict_plot['u'][2], extend='both') 
plt.title(u'(g) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 8)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_jja_regcm_era5_v[0], levels=dict_plot['v'][1], cmap=dict_plot['v'][2], extend='both') 
plt.title(u'(h) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 9)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_jja_regcm_era5_q[0], levels=dict_plot['q'][1], cmap=dict_plot['q'][2], extend='both') 
plt.title(u'(i) RegCM5-ERA5 JJA', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 10)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_son_regcm_era5_u[0], levels=dict_plot['u'][1], cmap=dict_plot['u'][2], extend='both') 
plt.title(u'(j) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 11)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_son_regcm_era5_v[0], levels=dict_plot['v'][1], cmap=dict_plot['v'][2], extend='both') 
plt.title(u'(k) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(4, 3, 12)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mbe_son_regcm_era5_q[0], levels=dict_plot['q'][1], cmap=dict_plot['q'][2], extend='both') 
plt.title(u'(l) RegCM5-ERA5 SON', loc='left', fontsize=font_size, fontweight='bold')
	
# Path out to save figure
path_out = '{0}/user/mdasilva/sam_3km/figs'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_SAM-3km_RegCM5_2018-2021.png'.format(var, level)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
