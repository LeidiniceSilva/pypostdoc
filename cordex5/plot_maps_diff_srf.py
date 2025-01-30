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

var = 'clh'
domain = 'CSAM-3'
idt, fdt = '1999', '1999'
dt = '{0}-{1}'.format(idt, fdt)

path = '/leonardo/home/userexternal/mdasilva/leonardo_work'

	
def import_rcm(exp, param, dataset, season):

	arq   = '{0}/{1}/postproc/rcm/{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, exp, param, domain, dataset, season, dt)	    
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	
		
def configure_subplot(ax):

	ax.set_xticks(np.arange(-82,-34,12), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(-38,-10,6), crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.tick_params(axis='x', labelsize=6, labelcolor='black')
	ax.tick_params(axis='y', labelsize=6, labelcolor='black')
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.coastlines()
	

# Import model and obs dataset
dict_var = {'pr': ['pre', 'precip', 'cmorph', 'precipitation', 'pr'],
'tas': ['tmp', 'tas'],
'tasmax': ['tmx', 'tmax', 'tasmax'],
'tasmin': ['tmn', 'tmin', 'tasmin'],
'evspsblpot': ['pev'],
'rsnl': ['msnlwrf'],
'rsns': ['msnswrf'],
'clt': ['cld', 'clt'],
'cll': ['cll'],
'clm': ['clm'],
'clh': ['clh']}

if var == 'tas' or var == 'tasmax' or var == 'tasmin':
	lat, lon, rcm3_djf = import_rcm('CORDEX5_v1', var, 'RegCM5', 'DJF')
	lat, lon, rcm3_mam = import_rcm('CORDEX5_v1', var, 'RegCM5', 'MAM')
	lat, lon, rcm3_jja = import_rcm('CORDEX5_v1', var, 'RegCM5', 'JJA')
	lat, lon, rcm3_son = import_rcm('CORDEX5_v1', var, 'RegCM5', 'SON')

	lat, lon, rcm3_moist_djf = import_rcm('CORDEX5_MOIST', var, 'RegCM5', 'DJF')
	lat, lon, rcm3_moist_mam = import_rcm('CORDEX5_MOIST', var, 'RegCM5', 'MAM')
	lat, lon, rcm3_moist_jja = import_rcm('CORDEX5_MOIST', var, 'RegCM5', 'JJA')
	lat, lon, rcm3_moist_son = import_rcm('CORDEX5_MOIST', var, 'RegCM5', 'SON')

	mbe_djf_rcm3 = compute_mbe(rcm3_moist_djf[0], rcm3_djf[0])
	mbe_mam_rcm3 = compute_mbe(rcm3_moist_mam[0], rcm3_mam[0])
	mbe_jja_rcm3 = compute_mbe(rcm3_moist_jja[0], rcm3_jja[0])
	mbe_son_rcm3 = compute_mbe(rcm3_moist_son[0], rcm3_son[0])

elif var == 'pr' or var == 'clt':
	lat, lon, rcm3_djf = import_rcm('CORDEX5_v1', var, 'RegCM5', 'DJF')
	lat, lon, rcm3_mam = import_rcm('CORDEX5_v1', var, 'RegCM5', 'MAM')
	lat, lon, rcm3_jja = import_rcm('CORDEX5_v1', var, 'RegCM5', 'JJA')
	lat, lon, rcm3_son = import_rcm('CORDEX5_v1', var, 'RegCM5', 'SON')

	lat, lon, rcm3_moist_djf = import_rcm('CORDEX5_MOIST', var, 'RegCM5', 'DJF')
	lat, lon, rcm3_moist_mam = import_rcm('CORDEX5_MOIST', var, 'RegCM5', 'MAM')
	lat, lon, rcm3_moist_jja = import_rcm('CORDEX5_MOIST', var, 'RegCM5', 'JJA')
	lat, lon, rcm3_moist_son = import_rcm('CORDEX5_MOIST', var, 'RegCM5', 'SON')
		
	mbe_djf_rcm3 = compute_mbe(rcm3_moist_djf, rcm3_djf)
	mbe_mam_rcm3 = compute_mbe(rcm3_moist_mam, rcm3_mam)
	mbe_jja_rcm3 = compute_mbe(rcm3_moist_jja, rcm3_jja)
	mbe_son_rcm3 = compute_mbe(rcm3_moist_son, rcm3_son)

else:
	lat, lon, rcm3_djf = import_rcm('CORDEX5_v1', var, 'RegCM5', 'DJF')
	lat, lon, rcm3_mam = import_rcm('CORDEX5_v1', var, 'RegCM5', 'MAM')
	lat, lon, rcm3_jja = import_rcm('CORDEX5_v1', var, 'RegCM5', 'JJA')
	lat, lon, rcm3_son = import_rcm('CORDEX5_v1', var, 'RegCM5', 'SON')

	lat, lon, rcm3_moist_djf = import_rcm('CORDEX5_MOIST', var, 'RegCM5', 'DJF')
	lat, lon, rcm3_moist_mam = import_rcm('CORDEX5_MOIST', var, 'RegCM5', 'MAM')
	lat, lon, rcm3_moist_jja = import_rcm('CORDEX5_MOIST', var, 'RegCM5', 'JJA')
	lat, lon, rcm3_moist_son = import_rcm('CORDEX5_MOIST', var, 'RegCM5', 'SON')
		
	mbe_djf_rcm3 = compute_mbe(rcm3_moist_djf/100, rcm3_djf/100)
	mbe_mam_rcm3 = compute_mbe(rcm3_moist_mam/100, rcm3_mam/100)
	mbe_jja_rcm3 = compute_mbe(rcm3_moist_jja/100, rcm3_jja/100)
	mbe_son_rcm3 = compute_mbe(rcm3_moist_son/100, rcm3_son/100)

# Plot figure  
fig, axes = plt.subplots(1, 4, figsize=(10, 4), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()
font_size = 6

plot_data = {'Plot 1': {'data': mbe_djf_rcm3, 'title': '(a) RCM3(moist)-RCM3(crtl) DJF'},
'Plot 2': {'data': mbe_mam_rcm3, 'title': '(b) RCM3(moist)-RCM3(crtl) MAM'},
'Plot 3': {'data': mbe_jja_rcm3, 'title': '(c) RCM3(moist)-RCM3(crtl) JJA'},
'Plot 4': {'data': mbe_son_rcm3, 'title': '(d) RCM3(moist)-RCM3(crtl) SON'}}

dict_plot = {'pr': ['Precipitation (mm d$^-$$^1$)', np.arange(-4, 4.2, 0.2), cm.PiYG],
'tas': ['Air temperature (°C)', np.arange(-1, 1.1, 0.1), cm.PiYG],
'tasmax': ['Maximum air temperature (°C)', np.arange(-1, 1.1, 0.1), cm.PiYG],
'tasmin': ['Minimum air temperature (°C)', np.arange(-1, 1.1, 0.1), cm.PiYG],
'evspsblpot': ['Potential evaporation (mm d$^-$$^1$)', np.arange(-2, 2.2, 0.2), cm.PiYG],
'rsnl': ['Surface net upward longwave flux (W mm$^-$$^2$)', np.arange(-2, 2.2, 0.2), cm.PiYG],
'rsns': ['Surface net downward shortwave flux (W mm$^-$$^2$)', np.arange(-2, 2.2, 0.2), cm.PiYG],
'clt': ['Total cloud cover (%)', np.arange(-10, 11, 1), cm.PiYG],
'cll': ['Low cloud cover (0-1)', np.arange(-0.5, 0.55, 0.05), cm.PiYG],
'clm': ['Medium cloud cover (0-1)', np.arange(-0.5, 0.55, 0.05), cm.PiYG],
'clh': ['High cloud cover (0-1)', np.arange(-0.5, 0.55, 0.05), cm.PiYG]}

for ax, (key, value) in zip(axes, plot_data.items()):
        data = value['data']
        title = value['title']
    
        contour = ax.contourf(lon, lat, data[0], transform=ccrs.PlateCarree(), levels=dict_plot[var][1], cmap=dict_plot[var][2], extend='neither')
        ax.set_title(title, loc='left', fontsize=font_size, fontweight='bold')
        configure_subplot(ax)

# Set colobar
cbar = fig.colorbar(contour, ax=fig.axes, orientation='horizontal', pad=0.1, aspect=50)
cbar.set_label('{0}'.format(dict_plot[var][0]), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/CORDEX5_MOIST/figs'.format(path)
name_out = 'pyplt_maps_diff_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


	
	

	
	
