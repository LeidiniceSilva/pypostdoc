# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot clim maps"

import os
import netCDF4
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat

from cartopy import config
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

var = 'pr'
domain = 'EUR-11'
dt = '2000-2001'
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11'
	
			
def import_obs(param, dataset, season):

	arq   = '{0}/postproc/obs/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

	return lat, lon, mean


def import_rcm(param, dataset, season):

	arq   = '{0}/postproc/rcm/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

	return lat, lon, mean
	

def configure_subplot(ax):

    ax.set_xticks(np.arange(-40,65,25), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(30,85,15), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.tick_params(axis='x', labelsize=6, labelcolor='black') 
    ax.tick_params(axis='y', labelsize=6, labelcolor='black') 
    ax.grid(c='k', ls='--', alpha=0.4)
    ax.coastlines()

# Import model and obs dataset
dict_var = {'pr': ['rr', 'precipitation', 'precip']}

lat, lon, eobs_djf = import_obs(dict_var[var][0], 'EOBS', 'DJF')
lat, lon, eobs_mam = import_obs(dict_var[var][0], 'EOBS', 'MAM')
lat, lon, eobs_jja = import_obs(dict_var[var][0], 'EOBS', 'JJA')
lat, lon, eobs_son = import_obs(dict_var[var][0], 'EOBS', 'SON')

lat, lon, mswep_djf = import_obs(dict_var[var][1], 'MSWEP', 'DJF')
lat, lon, mswep_mam = import_obs(dict_var[var][1], 'MSWEP', 'MAM')
lat, lon, mswep_jja = import_obs(dict_var[var][1], 'MSWEP', 'JJA')
lat, lon, mswep_son = import_obs(dict_var[var][1], 'MSWEP', 'SON')

lat, lon, cpc_djf = import_obs(dict_var[var][2], 'CPC', 'DJF')
lat, lon, cpc_mam = import_obs(dict_var[var][2], 'CPC', 'MAM')
lat, lon, cpc_jja = import_obs(dict_var[var][2], 'CPC', 'JJA')
lat, lon, cpc_son = import_obs(dict_var[var][2], 'CPC', 'SON')

lat, lon, noto_djf = import_rcm(var, 'NoTo-Europe_RegCM5', 'DJF')
lat, lon, noto_mam = import_rcm(var, 'NoTo-Europe_RegCM5', 'MAM')
lat, lon, noto_jja = import_rcm(var, 'NoTo-Europe_RegCM5', 'JJA')
lat, lon, noto_son = import_rcm(var, 'NoTo-Europe_RegCM5', 'SON')

lat, lon, wdm7_djf = import_rcm(var, 'WDM7-Europe_RegCM5', 'DJF')
lat, lon, wdm7_mam = import_rcm(var, 'WDM7-Europe_RegCM5', 'MAM')
lat, lon, wdm7_jja = import_rcm(var, 'WDM7-Europe_RegCM5', 'JJA')
lat, lon, wdm7_son = import_rcm(var, 'WDM7-Europe_RegCM5', 'SON')

lat, lon, wsm7_djf = import_rcm(var, 'WSM7-Europe_RegCM5', 'DJF')
lat, lon, wsm7_mam = import_rcm(var, 'WSM7-Europe_RegCM5', 'MAM')
lat, lon, wsm7_jja = import_rcm(var, 'WSM7-Europe_RegCM5', 'JJA')
lat, lon, wsm7_son = import_rcm(var, 'WSM7-Europe_RegCM5', 'SON')

lat, lon, wsm5_djf = import_rcm(var, 'WSM5-Europe_RegCM5', 'DJF')
lat, lon, wsm5_mam = import_rcm(var, 'WSM5-Europe_RegCM5', 'MAM')
lat, lon, wsm5_jja = import_rcm(var, 'WSM5-Europe_RegCM5', 'JJA')
lat, lon, wsm5_son = import_rcm(var, 'WSM5-Europe_RegCM5', 'SON')

# Plot figure
fig = plt.figure(figsize=(15, 6))

color = ['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']
dict_plot = {'pr': ['Precipitation (mm d$^-$$^1$)', np.arange(0, 18, 1), matplotlib.colors.ListedColormap(color)]}
font_size = 6 

ax1 = fig.add_subplot(4, 5, 1, projection=ccrs.PlateCarree())
plt1 = ax1.contourf(lon, lat, eobs_djf, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
configure_subplot(ax1)
ax1.set_title(u'(a) EOBS DJF', loc='left', fontsize=font_size, fontweight='bold')

ax2 = fig.add_subplot(4, 5, 2, projection=ccrs.PlateCarree())
ax2.contourf(lon, lat, noto_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
configure_subplot(ax2)
ax2.set_title(u'(b) NoTo DJF', loc='left', fontsize=font_size, fontweight='bold')

ax3 = fig.add_subplot(4, 5, 3, projection=ccrs.PlateCarree())
ax3.contourf(lon, lat, wdm7_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
configure_subplot(ax3)
ax3.set_title(u'(c) WDM7 DJF', loc='left', fontsize=font_size, fontweight='bold')

ax4 = fig.add_subplot(4, 5, 4, projection=ccrs.PlateCarree())
ax4.contourf(lon, lat, wsm7_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
configure_subplot(ax4)
ax4.set_title(u'(d) WSM7 DJF', loc='left', fontsize=font_size, fontweight='bold')

ax5 = fig.add_subplot(4, 5, 5, projection=ccrs.PlateCarree())
ax5.contourf(lon, lat, wsm5_djf, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
configure_subplot(ax5)
ax5.set_title(u'(e) WSM5 DJF', loc='left', fontsize=font_size, fontweight='bold')

ax6 = fig.add_subplot(4, 5, 6, projection=ccrs.PlateCarree())
ax6.contourf(lon, lat, eobs_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
configure_subplot(ax6)
ax6.set_title(u'(f) EOBS MAM', loc='left', fontsize=font_size, fontweight='bold')

ax7 = fig.add_subplot(4, 5, 7, projection=ccrs.PlateCarree())
ax7.contourf(lon, lat, noto_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
configure_subplot(ax7)
ax7.set_title(u'(g) NoTo MAM', loc='left', fontsize=font_size, fontweight='bold')

ax8 = fig.add_subplot(4, 5, 8, projection=ccrs.PlateCarree())
ax8.contourf(lon, lat, wdm7_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
configure_subplot(ax8)
ax8.set_title(u'(h) WDM7 MAM', loc='left', fontsize=font_size, fontweight='bold')

ax9 = fig.add_subplot(4, 5, 9, projection=ccrs.PlateCarree())
ax9.contourf(lon, lat, wsm7_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
configure_subplot(ax9)
ax9.set_title(u'(i) WSM7 MAM', loc='left', fontsize=font_size, fontweight='bold')

ax10 = fig.add_subplot(4, 5, 10, projection=ccrs.PlateCarree())
ax10.contourf(lon, lat, wsm5_mam, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
configure_subplot(ax10)
ax10.set_title(u'(j) WSM5 MAM', loc='left', fontsize=font_size, fontweight='bold')

ax11 = fig.add_subplot(4, 5, 11, projection=ccrs.PlateCarree())
ax11.contourf(lon, lat, eobs_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
configure_subplot(ax11)
ax11.set_title(u'(k) EOBS JJA', loc='left', fontsize=font_size, fontweight='bold')

ax12 = fig.add_subplot(4, 5, 12, projection=ccrs.PlateCarree())
ax12.contourf(lon, lat, noto_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
configure_subplot(ax12)
ax12.set_title(u'(l) NoTo JJA', loc='left', fontsize=font_size, fontweight='bold')

ax13 = fig.add_subplot(4, 5, 13, projection=ccrs.PlateCarree())
ax13.contourf(lon, lat, wdm7_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
configure_subplot(ax13)
ax13.set_title(u'(m) WDM7 JJA', loc='left', fontsize=font_size, fontweight='bold')

ax14 = fig.add_subplot(4, 5, 14, projection=ccrs.PlateCarree())
ax14.contourf(lon, lat, wsm7_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
configure_subplot(ax14)
ax14.set_title(u'(n) WSM7 JJA', loc='left', fontsize=font_size, fontweight='bold')

ax15 = fig.add_subplot(4, 5, 15, projection=ccrs.PlateCarree())
ax15.contourf(lon, lat, wsm5_jja, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
configure_subplot(ax15)
ax15.set_title(u'(o) WSM5 JJA', loc='left', fontsize=font_size, fontweight='bold')

ax16 = fig.add_subplot(4, 5, 16, projection=ccrs.PlateCarree())
ax16.contourf(lon, lat, eobs_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
configure_subplot(ax16)
ax16.set_title(u'(p) EOBS SON', loc='left', fontsize=font_size, fontweight='bold')

ax17 = fig.add_subplot(4, 5, 17, projection=ccrs.PlateCarree())
ax17.contourf(lon, lat, noto_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
configure_subplot(ax17)
ax17.set_title(u'(q) NoTo SON', loc='left', fontsize=font_size, fontweight='bold')

ax18 = fig.add_subplot(4, 5, 18, projection=ccrs.PlateCarree())
ax18.contourf(lon, lat, wdm7_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
configure_subplot(ax18)
ax18.set_title(u'(r) WDM7 SON', loc='left', fontsize=font_size, fontweight='bold')

ax19 = fig.add_subplot(4, 5, 19, projection=ccrs.PlateCarree())
ax19.contourf(lon, lat, wsm7_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
configure_subplot(ax19)
ax19.set_title(u'(s) WSM7 SON', loc='left', fontsize=font_size, fontweight='bold')

ax20 = fig.add_subplot(4, 5, 20, projection=ccrs.PlateCarree())
ax20.contourf(lon, lat, wsm5_son, levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
configure_subplot(ax20)
ax20.set_title(u'(t) WSM5 SON', loc='left', fontsize=font_size, fontweight='bold')

# Set colobar
cbar = fig.colorbar(plt1, ax=fig.axes, pad=0.025, aspect=50)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/figs/ctrl'.format(path)
name_out = 'pyplt_maps_clim_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


