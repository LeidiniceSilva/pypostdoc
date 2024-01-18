# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 04, 2023"
__description__ = "This script plot vertical profile"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

var = 'clw'
domain = 'SESA'
path='/marconi/home/userexternal/mdasilva'


def import_grid(param, domain, dataset, season):

	arq   = '{0}/user/mdasilva/sam_3km/post/{1}_{2}_{3}_{4}_2018-2021_lonlat.nc'.format(path, param, domain, dataset, season)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	value = var[:][:,:,:,:]
	mean = np.nanmean(np.nanmean(value, axis=2), axis=2)
	
	return mean
	
	
# Import model and obs dataset 
if var == 'cl':
	dict_var = {'cl': ['cc']}

	regcm_djf = import_grid(var, domain, 'RegCM5', 'DJF')
	regcm_mam = import_grid(var, domain, 'RegCM5', 'MAM')
	regcm_jja = import_grid(var, domain, 'RegCM5', 'JJA')
	regcm_son = import_grid(var, domain, 'RegCM5', 'SON')
		
	era5_djf = import_grid(dict_var[var][0], domain, 'ERA5', 'DJF')
	era5_mam = import_grid(dict_var[var][0], domain, 'ERA5', 'MAM')
	era5_jja = import_grid(dict_var[var][0], domain, 'ERA5', 'JJA')
	era5_son = import_grid(dict_var[var][0], domain, 'ERA5', 'SON')
	
	era5_djf_lev = era5_djf[0]*100
	era5_mam_lev = era5_mam[0]*100
	era5_jja_lev = era5_jja[0]*100
	era5_son_lev = era5_son[0]*100

	regcm_djf_lev = regcm_djf[0]*100
	regcm_mam_lev = regcm_mam[0]*100
	regcm_jja_lev = regcm_jja[0]*100
	regcm_son_lev = regcm_son[0]*100

elif var == 'cli':
	dict_var = {'cli': ['ciwc']}
	
	regcm_djf = import_grid(var, domain, 'RegCM5', 'DJF')
	regcm_mam = import_grid(var, domain, 'RegCM5', 'MAM')
	regcm_jja = import_grid(var, domain, 'RegCM5', 'JJA')
	regcm_son = import_grid(var, domain, 'RegCM5', 'SON')
		
	era5_djf = import_grid(dict_var[var][0], domain, 'ERA5', 'DJF')
	era5_mam = import_grid(dict_var[var][0], domain, 'ERA5', 'MAM')
	era5_jja = import_grid(dict_var[var][0], domain, 'ERA5', 'JJA')
	era5_son = import_grid(dict_var[var][0], domain, 'ERA5', 'SON')
	
	era5_djf_lev = era5_djf[0]*1000000
	era5_mam_lev = era5_mam[0]*1000000
	era5_jja_lev = era5_jja[0]*1000000
	era5_son_lev = era5_son[0]*1000000

	regcm_djf_lev = regcm_djf[0]*1000000
	regcm_mam_lev = regcm_mam[0]*1000000
	regcm_jja_lev = regcm_jja[0]*1000000
	regcm_son_lev = regcm_son[0]*1000000
	
else:
	dict_var = {'clw': ['clwc']}
	
	regcm_djf = import_grid(var, domain, 'RegCM5', 'DJF')
	regcm_mam = import_grid(var, domain, 'RegCM5', 'MAM')
	regcm_jja = import_grid(var, domain, 'RegCM5', 'JJA')
	regcm_son = import_grid(var, domain, 'RegCM5', 'SON')
		
	era5_djf = import_grid(dict_var[var][0], domain, 'ERA5', 'DJF')
	era5_mam = import_grid(dict_var[var][0], domain, 'ERA5', 'MAM')
	era5_jja = import_grid(dict_var[var][0], domain, 'ERA5', 'JJA')
	era5_son = import_grid(dict_var[var][0], domain, 'ERA5', 'SON')

	era5_djf_lev = era5_djf[0]*1000000
	era5_mam_lev = era5_mam[0]*1000000
	era5_jja_lev = era5_jja[0]*1000000
	era5_son_lev = era5_son[0]*1000000

	regcm_djf_lev = regcm_djf[0]*1000000
	regcm_mam_lev = regcm_mam[0]*1000000
	regcm_jja_lev = regcm_jja[0]*1000000
	regcm_son_lev = regcm_son[0]*1000000
	
# Plot figure   
fig = plt.figure(figsize=(10, 6))
font_size = 8

dict_plot = {
'cl': ['Cloud fraction (%)', 0, 20, np.arange(0, 24, 4)],
'cli': ['Cloud ice (mg kg⁻¹)', 0, 25, np.arange(0, 30, 5)],
'clw': ['Cloud liquid water (mg kg⁻¹)', 0, 25, np.arange(0, 30, 5)]
}

levels_i = (1000,975,950,925,900,875,850,825,800,775,750,700,650,600,550,500,450,400,350,300,250,225,200,175,150,125,100,70,50,30,20,10,7,5,3,2,1)
levels_ii = (1000,925,850,700,600,500,400,300,250,200,150,100)

ax = fig.add_subplot(1, 4, 1)
ax.plot(era5_djf_lev[::-1], levels_i, color='black', label='ERA5', linewidth=1, linestyle='--')
plt.plot(regcm_djf_lev, levels_ii, color='red', label='RegCM5', linewidth=1)
plt.title(u'(a) DJF', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')
plt.ylabel('Level pressure (hPa)', fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()
plt.legend(loc=1, ncol=1, fontsize=font_size)

ax = fig.add_subplot(1, 4, 2)
plt.plot(era5_mam_lev[::-1], levels_i, color='black', label='ERA5', linewidth=1, linestyle='--')
plt.plot(regcm_mam_lev, levels_ii, color='red', label='RegCM5', linewidth=1)
plt.title(u'(b) MAM', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.grid(linestyle='--')
plt.setp(ax.get_yticklabels(), visible=False)
plt.gca().invert_yaxis()

ax = fig.add_subplot(1, 4, 3)
plt.plot(era5_jja_lev[::-1], levels_i, color='black', label='ERA5', linewidth=1, linestyle='--')
plt.plot(regcm_jja_lev, levels_ii, color='red', label='RegCM5', linewidth=1)
plt.title(u'(c) JJA', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.grid(linestyle='--')
plt.setp(ax.get_yticklabels(), visible=False)
plt.gca().invert_yaxis()

ax = fig.add_subplot(1, 4, 4)
plt.plot(era5_son_lev[::-1], levels_i, color='black', label='ERA5', linewidth=1, linestyle='--')
plt.plot(regcm_son_lev, levels_ii, color='red', label='RegCM5', linewidth=1)
plt.title(u'(d) SON', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.grid(linestyle='--')
plt.setp(ax.get_yticklabels(), visible=False)
plt.gca().invert_yaxis()

# Path out to save figure
path_out = '{0}/user/mdasilva/sam_3km/figs'.format(path)
name_out = 'pyplt_vertical_profile_{0}_{1}_RegCM5_2018-2021.png'.format(var, domain)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
