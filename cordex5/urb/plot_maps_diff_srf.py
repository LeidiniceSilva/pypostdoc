# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
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

var = 'clt'
domain = 'CSAM-3'
idt, fdt = '2000', '2000'
dt = '{0}-{1}'.format(idt, fdt)
font_size = 6

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5'


def import_rcm_v1(param, domain, dataset, season):

        arq   = '{0}/postproc/evaluate/rcm/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)
        data  = netCDF4.Dataset(arq)
        var   = data.variables[param][:]
        lat   = data.variables['lat'][:]
        lon   = data.variables['lon'][:]
        mean = var[:][0,:,:]

        return lat, lon, mean


def import_rcm_v2(param, domain, dataset, season):

	arq   = '{0}/postproc/evaluate/rcm_urb/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]

	return lat, lon, mean


# Import model and obs dataset
dict_var = {'pr': ['pre', 'precip', 'cmorph', 'precipitation', 'tp'],
'tas': ['tmp', 't2m'],
'tasmax': ['tmx', 'tmax', 'tasmax'],
'tasmin': ['tmn', 'tmin', 'tasmin'],
'clt': ['cld', 'tcc'],
'cll': ['lcc'],
'clm': ['mcc'],
'clh': ['hcc'],
'evspsblpot': ['pev'],
'rlds': ['msdwlwrf']}

lat, lon, regcm_djf_v1 = import_rcm_v1(var, domain, 'RegCM5', 'DJF')
lat, lon, regcm_mam_v1 = import_rcm_v1(var, domain, 'RegCM5', 'MAM')
lat, lon, regcm_jja_v1 = import_rcm_v1(var, domain, 'RegCM5', 'JJA')
lat, lon, regcm_son_v1 = import_rcm_v1(var, domain, 'RegCM5', 'SON')
        
lat, lon, regcm_djf_v2 = import_rcm_v2(var, domain, 'RegCM5', 'DJF')
lat, lon, regcm_mam_v2 = import_rcm_v2(var, domain, 'RegCM5', 'MAM')
lat, lon, regcm_jja_v2 = import_rcm_v2(var, domain, 'RegCM5', 'JJA')
lat, lon, regcm_son_v2 = import_rcm_v2(var, domain, 'RegCM5', 'SON')

mbe_djf_regcm = compute_mbe(regcm_djf_v2, regcm_djf_v1)
mbe_mam_regcm = compute_mbe(regcm_mam_v2, regcm_mam_v1)
mbe_jja_regcm = compute_mbe(regcm_jja_v2, regcm_jja_v1)
mbe_son_regcm = compute_mbe(regcm_son_v2, regcm_son_v1)	
	
# Plot figure   
def configure_subplot(ax):

	lon_min = np.round(np.min(lon), 1)
	lon_max = np.round(np.max(lon), 1)
	lat_min = np.round(np.min(lat), 1)
	lat_max = np.round(np.max(lat), 1)
	ax.set_extent([np.min(lon), np.max(lon), np.min(lat), np.max(lat)], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(lon_min,lon_max,10), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(lat_min,lat_max,5), crs=ccrs.PlateCarree())

	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(6)

	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.add_feature(cfeat.BORDERS, linewidth=0.5)
	ax.coastlines(linewidth=0.5)

	if var == 'evspsblpot':
		ax.add_feature(cfeat.OCEAN, facecolor='white', zorder=1) 

fig, axes = plt.subplots(4, 1, figsize=(10, 10), subplot_kw={'projection': ccrs.PlateCarree()})

dict_plot = {'pr': ['Diff of  precipitation (mm d$^-$$^1$)', np.arange(-5, 5.5, 0.5), cm.BrBG],
'tas': ['Diff of air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'tasmax': ['Diff of maximum air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'tasmin': ['Diff of minimum air temperature (°C)', np.arange(-10, 11, 1), cm.bwr],
'clt': ['Diff of total cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'cll': ['Diff of low cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'clm': ['Diff of medium cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'clh': ['Diff of high cloud cover (0-1)', np.arange(-0.7, 0.8, 0.1), cm.RdGy],
'evspsblpot': ['Diff of potential evapotranspiration (mm d$^-$$^1$)', np.arange(-5, 5.5, 0.5), cm.bwr],
'rlds': ['Diff of surface downwelling longwave radiation (W mm$^-$$^2$)', np.arange(-60, 55, 5), cm.RdBu_r]}

ax1 = axes[0]
plt_map = ax1.contourf(lon, lat, mbe_djf_regcm, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
ax1.set_title(u'(a) RegCM5_URB - RegCM5 DJF', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

ax2 = axes[1]
plt_map = ax2.contourf(lon, lat, mbe_mam_regcm, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
ax2.set_title(u'(b) RegCM5_URB - RegCM5 MAM', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

ax3 = axes[2] 
plt_map = ax3.contourf(lon, lat, mbe_jja_regcm, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
ax3.set_title(u'(c) RegCM5_URB - RegCM5 JJA', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

ax4 = axes[3]
plt_map = ax4.contourf(lon, lat, mbe_son_regcm, transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither') 
ax4.set_title(u'(d) RegCM5_URB - RegCM5 SON', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax4)

# Set colobar
cbar = fig.colorbar(plt_map, ax=fig.axes, pad=0.02, aspect=50)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '{0}/figs/evaluate/rcm_urb'.format(path)
name_out = 'pyplt_maps_diff_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()



	
	

	
	
