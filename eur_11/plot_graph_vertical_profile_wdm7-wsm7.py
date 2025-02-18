# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 12, 2024"
__description__ = "This script plot vertical profile"

import os
import netCDF4
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt

var = 'hail'
domain = 'EUR-11'
dt = '2000-2001'
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11'


def import_obs(param, dataset, season):
		
	arq   = '{0}/postproc/obs/{1}_{2}_FPS_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)	     
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	value = var[:][:,:,:,:]
	mean = np.nanmean(np.nanmean(value, axis=2), axis=2)
	
	return mean


def import_rcm(param, dataset, season):

	arq   = '{0}/postproc/rcm/{1}_{2}_FPS_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	value = var[:][:,:,:,:]
	mean = np.nanmean(np.nanmean(value, axis=2), axis=2)
	
	return mean
	
	
# Import model and obs dataset
wdm7_djf = import_rcm(var, 'WDM7-Europe', 'DJF')
wdm7_mam = import_rcm(var, 'WDM7-Europe', 'MAM')
wdm7_jja = import_rcm(var, 'WDM7-Europe', 'JJA')
wdm7_son = import_rcm(var, 'WDM7-Europe', 'SON')

wsm7_djf = import_rcm(var, 'WSM7-Europe', 'DJF')
wsm7_mam = import_rcm(var, 'WSM7-Europe', 'MAM')
wsm7_jja = import_rcm(var, 'WSM7-Europe', 'JJA')
wsm7_son = import_rcm(var, 'WSM7-Europe', 'SON')

factor = 1000000
wdm7_djf_ = wdm7_djf[0]*factor
wdm7_mam_ = wdm7_mam[0]*factor
wdm7_jja_ = wdm7_jja[0]*factor
wdm7_son_ = wdm7_son[0]*factor

wsm7_djf_ = wsm7_djf[0]*factor
wsm7_mam_ = wsm7_mam[0]*factor
wsm7_jja_ = wsm7_jja[0]*factor
wsm7_son_ = wsm7_son[0]*factor

# Plot figure  
fig = plt.figure(figsize=(12, 4))
font_size = 8

dict_plot = {
'gra': ['Mass fraction of graupel (mg kg$^-$$^1$)', 0, 1, np.arange(0, 1.1, 0.1)],
'hail': ['Mass fraction of hail (mg kg$^-$$^1$)', 0, 1, np.arange(0, 1.1, 0.1)]
}

levels_ii = (1000,925,850,700,600,500,400,300,250,200,150,100)

ax = fig.add_subplot(1, 4, 1)
plt.plot(wdm7_djf_, levels_ii, color='green', label='WDM7', linewidth=1)
plt.plot(wsm7_djf_, levels_ii, color='magenta', label='WSM7', linewidth=1)
plt.title(u'(a) DJF', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')
plt.ylabel('Level pressure (hPa)', fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.yticks(fontsize=font_size)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()
plt.legend(loc=1, ncol=1, fontsize=font_size)

ax = fig.add_subplot(1, 4, 2)
plt.plot(wdm7_mam_, levels_ii, color='green', label='WDM7', linewidth=1)
plt.plot(wsm7_mam_, levels_ii, color='magenta', label='WSM7', linewidth=1)
plt.title(u'(b) MAM', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.yticks(fontsize=font_size)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()

ax = fig.add_subplot(1, 4, 3)
plt.plot(wdm7_jja_, levels_ii, color='green', label='WDM7', linewidth=1)
plt.plot(wsm7_jja_, levels_ii, color='magenta', label='WSM7', linewidth=1)
plt.title(u'(c) JJA', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.yticks(fontsize=font_size)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()

ax = fig.add_subplot(1, 4, 4)
plt.plot(wdm7_son_, levels_ii, color='green', label='WDM7', linewidth=1)
plt.plot(wsm7_son_, levels_ii, color='magenta', label='WSM7', linewidth=1)
plt.title(u'(d) SON', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.yticks(fontsize=font_size)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()

# Path out to save figure
path_out = '{0}/figs/totc'.format(path)
name_out = 'pyplt_graph_vertical_profile_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
