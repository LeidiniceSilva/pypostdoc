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

var = 'cls'
domain = 'EUR-11'
dt = '2000-2004'
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11'


def import_obs(param, dataset, season):
    
        arq   = '{0}/postproc/obs/{1}_{2}_FPS_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)
        data  = netCDF4.Dataset(arq)
        var   = data.variables[param][:]
        
        if param == 'crwc' or param == 'cswc':
                value = var[:][:,:,:,:]
        else:
            value = var[:][:,::-1,:,:]
        
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
dict_var = {'cl': ['cc'], 
'clw': ['clwc'],
'cli': ['ciwc'],
'clr': ['crwc'],
'cls': ['cswc'],
'rh': ['r'],
'hus': ['q']}

obs_djf = import_obs(dict_var[var][0], 'ERA5', 'DJF')
obs_mam = import_obs(dict_var[var][0], 'ERA5', 'MAM')
obs_jja = import_obs(dict_var[var][0], 'ERA5', 'JJA')
obs_son = import_obs(dict_var[var][0], 'ERA5', 'SON')

noto_djf = import_rcm(var, 'NoTo-Europe', 'DJF')
noto_mam = import_rcm(var, 'NoTo-Europe', 'MAM')
noto_jja = import_rcm(var, 'NoTo-Europe', 'JJA')
noto_son = import_rcm(var, 'NoTo-Europe', 'SON')

wdm7_djf = import_rcm(var, 'WDM7-Europe', 'DJF')
wdm7_mam = import_rcm(var, 'WDM7-Europe', 'MAM')
wdm7_jja = import_rcm(var, 'WDM7-Europe', 'JJA')
wdm7_son = import_rcm(var, 'WDM7-Europe', 'SON')

wsm7_djf = import_rcm(var, 'WSM7-Europe', 'DJF')
wsm7_mam = import_rcm(var, 'WSM7-Europe', 'MAM')
wsm7_jja = import_rcm(var, 'WSM7-Europe', 'JJA')
wsm7_son = import_rcm(var, 'WSM7-Europe', 'SON')

wsm5_djf = import_rcm(var, 'WSM5-Europe', 'DJF')
wsm5_mam = import_rcm(var, 'WSM5-Europe', 'MAM')
wsm5_jja = import_rcm(var, 'WSM5-Europe', 'JJA')
wsm5_son = import_rcm(var, 'WSM5-Europe', 'SON')

if var == 'cl':
	factor = 100
elif var == 'rh':
	factor = 1
elif var == 'hus':
	factor = 1000
else:
	factor = 1000000

obs_djf_ = obs_djf[0]*factor
obs_mam_ = obs_mam[0]*factor
obs_jja_ = obs_jja[0]*factor
obs_son_ = obs_son[0]*factor

noto_djf_ = noto_djf[0]*factor
noto_mam_ = noto_mam[0]*factor
noto_jja_ = noto_jja[0]*factor
noto_son_ = noto_son[0]*factor

wdm7_djf_ = wdm7_djf[0]*factor
wdm7_mam_ = wdm7_mam[0]*factor
wdm7_jja_ = wdm7_jja[0]*factor
wdm7_son_ = wdm7_son[0]*factor

wsm7_djf_ = wsm7_djf[0]*factor
wsm7_mam_ = wsm7_mam[0]*factor
wsm7_jja_ = wsm7_jja[0]*factor
wsm7_son_ = wsm7_son[0]*factor

wsm5_djf_ = wsm5_djf[0]*factor
wsm5_mam_ = wsm5_mam[0]*factor
wsm5_jja_ = wsm5_jja[0]*factor
wsm5_son_ = wsm5_son[0]*factor

# Plot figure  
fig = plt.figure(figsize=(12, 4))
font_size = 8

dict_plot = {
'cl': ['Cloud fraction (%)', 0, 40, np.arange(0, 44, 4)],
'clw': ['Cloud liquid water (mg kg$^-$$^1$)', 0, 50, np.arange(0, 55, 5)],
'cli': ['Cloud liquid ice (mg kg$^-$$^1$)', 0, 20, np.arange(0, 22, 2)],
'clr': ['Cloud liquid rain (mg kg$^-$$^1$)', 0, 10, np.arange(0, 11, 1)],
'cls': ['Cloud liquid snow (mg kg$^-$$^1$)', 0, 20, np.arange(0, 22, 2)],
'hus': ['Specific humidity (g kg$^-$$^1$)', 0, 10, np.arange(0, 11, 1)],
'rh': ['Relative humidity (%)', 0, 100, np.arange(0, 110, 10)]
}

levels_i = (1000,975,950,925,900,875,850,825,800,775,750,700,650,600,550,500,450,400,350,300,250,225,200,175,150,125,100,70,50,30,20,10,7,5,3,2,1)
levels_ii = (1000,925,850,700,600,500,400,300,250,200,150,100)

ax = fig.add_subplot(1, 4, 1)
plt.plot(obs_djf_, levels_i, color='black', label='ERA5', linewidth=1)
plt.plot(noto_djf_, levels_ii, color='blue', label='NoTo', linewidth=1)
plt.plot(wdm7_djf_, levels_ii, color='green', label='WDM7', linewidth=1)
plt.plot(wsm7_djf_, levels_ii, color='magenta', label='WSM7', linewidth=1)
plt.plot(wsm5_djf_, levels_ii, color='red',  label='WSM5', linewidth=1)
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
plt.plot(obs_mam_, levels_i, color='black', label='ERA5', linewidth=1)
plt.plot(noto_mam_, levels_ii, color='blue', label='NoTo', linewidth=1)
plt.plot(wdm7_mam_, levels_ii, color='green', label='WDM7', linewidth=1)
plt.plot(wsm7_mam_, levels_ii, color='magenta', label='WSM7', linewidth=1)
plt.plot(wsm5_mam_, levels_ii, color='red',  label='WSM5', linewidth=1)
plt.title(u'(b) MAM', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.yticks(fontsize=font_size)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()

ax = fig.add_subplot(1, 4, 3)
plt.plot(obs_jja_, levels_i, color='black', label='ERA5', linewidth=1)
plt.plot(noto_jja_, levels_ii, color='blue', label='NoTo', linewidth=1)
plt.plot(wdm7_jja_, levels_ii, color='green', label='WDM7', linewidth=1)
plt.plot(wsm7_jja_, levels_ii, color='magenta', label='WSM7', linewidth=1)
plt.plot(wsm5_jja_, levels_ii, color='red',  label='WSM5', linewidth=1)
plt.title(u'(c) JJA', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.yticks(fontsize=font_size)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()

ax = fig.add_subplot(1, 4, 4)
plt.plot(obs_son_, levels_i, color='black', label='ERA5', linewidth=1)
plt.plot(noto_son_, levels_ii, color='blue', label='NoTo', linewidth=1)
plt.plot(wdm7_son_, levels_ii, color='green', label='WDM7', linewidth=1)
plt.plot(wsm7_son_, levels_ii, color='magenta', label='WSM7', linewidth=1)
plt.plot(wsm5_son_, levels_ii, color='red',  label='WSM5', linewidth=1)
plt.title(u'(d) SON', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.yticks(fontsize=font_size)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()

# Path out to save figure
path_out = '{0}/figs/ctrl'.format(path)
name_out = 'pyplt_graph_vertical_profile_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
