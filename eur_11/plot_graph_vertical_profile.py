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

var = 'cl'
domain = 'EUR-11'
dt = '2000-2001'
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11'


def import_obs(param, dataset):

	if param == 'clfrac':
		param_ = 'cc'
	elif param == 'clliq':
		param_ = 'clwc'
	else:
		param_ = 'ciwc'
		
	arq   = '{0}/obs/{1}_{2}_FPS_{3}_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)	     
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param_][:] 
	value = var[:][:,:,:,:]
	mean = np.nanmean(np.nanmean(value, axis=2), axis=2)
	
	return mean


def import_rcm(param, dataset):

	arq   = '{0}/rcm/{1}_{2}_FPS_{3}_{4}_lonlat.nc'.format(path, param, domain, dataset, dt)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	value = var[:][:,:,:,:]
	mean = np.nanmean(np.nanmean(value, axis=2), axis=2)
	
	return mean
	
		
# Import model and obs dataset
dict_var = {'cl': ['clfrac'], 
'clw': ['clliq'],
'cli': ['clice']}
 
if var == 'cl':
	obs = import_obs(dict_var[var][0], 'ERA5')
	obs_lev = obs[0]

	noto = import_rcm(var, 'NoTo-Europe_RegCM5')
	noto_lev = noto[0]*100

	wdm7 = import_rcm(var, 'WDM7-Europe_RegCM5')
	wdm7_lev = wdm7[0]*100

	wsm7 = import_rcm(var, 'WSM7-Europe_RegCM5')
	wsm7_lev = wsm7[0]*100

	wsm5 = import_rcm(var, 'WSM5-Europe_RegCM5')
	wsm5_lev = wsm5[0]*100
else:
	obs = import_obs(dict_var[var][0], 'ERA5')
	obs_lev = obs[0]*1000000

	noto = import_rcm(var, 'NoTo-Europe_RegCM5')
	noto_lev = noto[0]*1000000

	wdm7 = import_rcm(var, 'WDM7-Europe_RegCM5')
	wdm7_lev = wdm7[0]*1000000

	wsm7 = import_rcm(var, 'WSM7-Europe_RegCM5')
	wsm7_lev = wsm7[0]*1000000

	wsm5 = import_rcm(var, 'WSM5-Europe_RegCM5')
	wsm5_lev = wsm5[0]*1000000
	
# Plot figure   
fig = plt.figure(figsize=(12, 4))
font_size = 8

dict_plot = {
'cl': ['Cloud fraction (%)', 0, 40, np.arange(0, 44, 4)],
'clw': ['Cloud liquid water (mg kg$^-$$^1$)', 0, 100, np.arange(0, 110, 10)],
'cli': ['Cloud ice (mg kg$^-$$^1$)', 0, 5, np.arange(0, 5.5, 0.5)]
}

levels_i = (1000,975,950,925,900,875,850,825,800,775,750,700,650,600,550,500,450,400,350,300,250,225,200,175,150,125,100,70,50,30,20,10,7,5,3,2,1)
levels_ii = (1000,925,850,700,600,500,400,300,250,200,150,100)

ax = fig.add_subplot(1, 4, 1)
plt.plot(era5_jan_lev[::-1], levels_i, color='black', label='ERA5', linewidth=1, linestyle='--')
plt.plot(noto_jan_lev, levels_ii, color='red', label='RegCM5', linewidth=1)
plt.title(u'(a) NoTo', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')
plt.ylabel('Level pressure (hPa)', fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()
plt.legend(loc=1, ncol=1, fontsize=font_size)

ax = fig.add_subplot(1, 4, 2)
ax.plot(era5_jan_lev[::-1], levels_i, color='black', label='ERA5', linewidth=1, linestyle='--')
plt.plot(wsm5_jan_lev, levels_ii, color='red', label='RegCM5', linewidth=1)
plt.title(u'(b) WSM5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.grid(linestyle='--')
plt.setp(ax.get_yticklabels(), visible=False)
plt.gca().invert_yaxis()

ax = fig.add_subplot(1, 4, 3)
ax.plot(era5_jan_lev[::-1], levels_i, color='black', label='ERA5', linewidth=1, linestyle='--')
plt.plot(wsm7_jan_lev, levels_ii, color='red', label='RegCM5', linewidth=1)
plt.title(u'(c) WSM7', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.grid(linestyle='--')
plt.setp(ax.get_yticklabels(), visible=False)
plt.gca().invert_yaxis()

ax = fig.add_subplot(1, 4, 4)
ax.plot(era5_jan_lev[::-1], levels_i, color='black', label='ERA5', linewidth=1, linestyle='--')
plt.plot(wdm7_jan_lev, levels_ii, color='red', label='RegCM5', linewidth=1)
plt.title(u'(d) WDM7', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.grid(linestyle='--')
plt.setp(ax.get_yticklabels(), visible=False)
plt.gca().invert_yaxis()

# Path out to save figure
path_out = '{0}/figs/ctrl'.format(path)
name_out = 'pyplt_graph_vertical_profile_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
