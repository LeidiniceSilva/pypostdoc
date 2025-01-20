# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot bias maps"

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
from import_climate_tools import compute_mbe

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

mbe_djf_noto_eobs = compute_mbe(noto_djf, eobs_djf)
mbe_mam_noto_eobs = compute_mbe(noto_mam, eobs_mam)
mbe_jja_noto_eobs = compute_mbe(noto_jja, eobs_jja)
mbe_son_noto_eobs = compute_mbe(noto_son, eobs_son)

mbe_djf_wdm7_eobs = compute_mbe(wdm7_djf, eobs_djf)
mbe_mam_wdm7_eobs = compute_mbe(wdm7_mam, eobs_mam)
mbe_jja_wdm7_eobs = compute_mbe(wdm7_jja, eobs_jja)
mbe_son_wdm7_eobs = compute_mbe(wdm7_son, eobs_son)

mbe_djf_wsm7_eobs = compute_mbe(wsm7_djf, eobs_djf)
mbe_mam_wsm7_eobs = compute_mbe(wsm7_mam, eobs_mam)
mbe_jja_wsm7_eobs = compute_mbe(wsm7_jja, eobs_jja)
mbe_son_wsm7_eobs = compute_mbe(wsm7_son, eobs_son)

mbe_djf_wsm5_eobs = compute_mbe(wsm5_djf, eobs_djf)
mbe_mam_wsm5_eobs = compute_mbe(wsm5_mam, eobs_mam)
mbe_jja_wsm5_eobs = compute_mbe(wsm5_jja, eobs_jja)
mbe_son_wsm5_eobs = compute_mbe(wsm5_son, eobs_son)

# Plot figure
fig, axes = plt.subplots(4, 4, figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()

dict_plot = {'pr': ['Bias of precipitation (mm d$^-$$^1$)', np.arange(-10, 11, 1), cm.BrBG]}
font_size = 8

plot_data = {
     'Plot 1': {'data': mbe_djf_noto_eobs, 'title': '(a) NoTo-EOBS DJF'},
     'Plot 2': {'data': mbe_djf_wdm7_eobs, 'title': '(b) WDM7-EOBS DJF'},
     'Plot 3': {'data': mbe_djf_wsm7_eobs, 'title': '(c) WSM7-EOBS DJF'},
     'Plot 4': {'data': mbe_djf_wsm5_eobs, 'title': '(d) WSM5-EOBS DJF'},
     'Plot 5': {'data': mbe_mam_noto_eobs, 'title': '(e) NoTo-EOBS MAM'},
     'Plot 6': {'data': mbe_mam_wdm7_eobs, 'title': '(f) WDM7-EOBS MAM'},
     'Plot 7': {'data': mbe_mam_wsm7_eobs, 'title': '(g) WSM7-EOBS MAM'},
     'Plot 8': {'data': mbe_mam_wsm5_eobs, 'title': '(h) WSM5-EOBS MAM'},
     'Plot 9': {'data': mbe_jja_noto_eobs, 'title': '(i) NoTo-EOBS JJA'},
     'Plot 10': {'data': mbe_jja_wdm7_eobs, 'title': '(j) WDM7-EOBS JJA'},
     'Plot 11': {'data': mbe_jja_wsm7_eobs, 'title': '(k) WSM7-EOBS JJA'},
     'Plot 12': {'data': mbe_jja_wsm5_eobs, 'title': '(l) WSM5-EOBS JJA'},
     'Plot 13': {'data': mbe_son_noto_eobs, 'title': '(m) NoTo-EOBS SON'},
     'Plot 14': {'data': mbe_son_wdm7_eobs, 'title': '(n) WDM7-EOBS SON'},
     'Plot 15': {'data': mbe_son_wsm7_eobs, 'title': '(o) WSM7-EOBS SON'},
     'Plot 16': {'data': mbe_son_wsm5_eobs, 'title': '(p) WSM5-EOBS SON'}
}

for ax, (key, value) in zip(axes, plot_data.items()):
    data = value['data']
    title = value['title']
    
    contour = ax.contourf(lon, lat, data, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
    ax.set_title(title, loc='left', fontsize=font_size, fontweight='bold')
    configure_subplot(ax)

# Set colobar
cbar = fig.colorbar(contour, ax=fig.axes, orientation='vertical', pad=0.025, aspect=50)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/figs/ctrl'.format(path)
name_out = 'pyplt_maps_bias_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
