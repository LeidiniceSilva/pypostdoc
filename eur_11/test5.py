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

var = 'cl'
domain = 'EUR-11'
path='/marconi/home/userexternal/mdasilva'


def import_obs(param):

	arq   = '{0}/user/mdasilva/EUR-11/post_evaluate/obs/{1}_{2}_FPS_ERA5_20000101_lonlat.nc'.format(path, param, domain)	     
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	value = var[:][:,:,:,:]
	mean = np.nanmean(np.nanmean(value, axis=2), axis=2)
	
	return mean


def import_rcm(param, experiment):

	arq   = '{0}/user/mdasilva/EUR-11/post_evaluate/rcm/{1}/{2}_{3}_FPS_RegCM5_20000101_lonlat.nc'.format(path, experiment, param, domain)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	value = var[:][:,:,:,:]
	mean = np.nanmean(np.nanmean(value, axis=2), axis=2)
	
	return mean
	
		
# Import model and obs dataset 
if var == 'cl':
	dict_var = {'cl': ['cc']}

	era5_jan_exp = import_obs(dict_var[var][0])
	era5_jan_lev = era5_jan_exp[0]*100

	regcm5_jan_exp1 = import_rcm(var, 'wdm7-Europe_v1')
	regcm5_jan_lev1 = regcm5_jan_exp1[0]*100
	
	regcm5_jan_exp2 = import_rcm(var, 'wdm7-Europe_v2')
	regcm5_jan_lev2 = regcm5_jan_exp2[0]*100

	regcm5_jan_exp3 = import_rcm(var, 'wdm7-Europe_v3')
	regcm5_jan_lev3 = regcm5_jan_exp3[0]*100

	regcm5_jan_exp4 = import_rcm(var, 'wdm7-Europe_v4')
	regcm5_jan_lev4 = regcm5_jan_exp4[0]*100

elif var == 'cli':
	dict_var = {'cli': ['ciwc']}
	
	era5_jan_exp = import_obs(dict_var[var][0])
	era5_jan_lev = era5_jan_exp[0]*1000000

	regcm5_jan_exp1 = import_rcm(var, 'wdm7-Europe_v1')
	regcm5_jan_lev1 = regcm5_jan_exp1[0]*1000000

	regcm5_jan_exp2 = import_rcm(var, 'wdm7-Europe_v2')
	regcm5_jan_lev2 = regcm5_jan_exp2[0]*1000000

	regcm5_jan_exp3 = import_rcm(var, 'wdm7-Europe_v3')
	regcm5_jan_lev3 = regcm5_jan_exp3[0]*1000000
	
	regcm5_jan_exp4 = import_rcm(var, 'wdm7-Europe_v4')
	regcm5_jan_lev4 = regcm5_jan_exp4[0]*1000000
		
else:
	dict_var = {'clw': ['clwc']}
	
	era5_jan_exp = import_obs(dict_var[var][0])
	era5_jan_lev = era5_jan_exp[0]*1000000

	regcm5_jan_exp1 = import_rcm(var, 'wdm7-Europe_v1')
	regcm5_jan_lev1 = regcm5_jan_exp1[0]*1000000

	regcm5_jan_exp2 = import_rcm(var, 'wdm7-Europe_v2')
	regcm5_jan_lev2 = regcm5_jan_exp2[0]*1000000
	
	regcm5_jan_exp3 = import_rcm(var, 'wdm7-Europe_v3')
	regcm5_jan_lev3 = regcm5_jan_exp3[0]*1000000
	
	regcm5_jan_exp4 = import_rcm(var, 'wdm7-Europe_v4')
	regcm5_jan_lev4 = regcm5_jan_exp4[0]*1000000


# Plot figure   
fig = plt.figure(figsize=(10, 6))
font_size = 8

dict_plot = {
'cl': ['Cloud fraction (%)', 0, 40, np.arange(0, 44, 4)],
'cli': ['Cloud ice (mg kg$^-$$^1$)', 0, 5, np.arange(0, 5.5, 0.5)],
'clw': ['Cloud liquid water (mg kg$^-$$^1$)', 0, 100, np.arange(0, 110, 10)]
}

levels_i = (1000,975,950,925,900,875,850,825,800,775,750,700,650,600,550,500,450,400,350,300,250,225,200,175,150,125,100,70,50,30,20,10,7,5,3,2,1)
levels_ii = (1000,925,850,700,600,500,400,300,250,200,150,100)

ax = fig.add_subplot(2, 2, 1)
plt.plot(era5_jan_lev[::-1], levels_i, color='black', label='ERA5', linewidth=1, linestyle='--')
plt.plot(regcm5_jan_lev1, levels_ii, color='red', label='RegCM5', linewidth=1)
plt.title(u'(a) WDM7_v1(ctrl)', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Level pressure (hPa)', fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()
plt.legend(loc=1, ncol=1, fontsize=font_size)

ax = fig.add_subplot(2, 2, 2)
ax.plot(era5_jan_lev[::-1], levels_i, color='black', label='ERA5', linewidth=1, linestyle='--')
plt.plot(regcm5_jan_lev2, levels_ii, color='red', label='RegCM5', linewidth=1)
plt.title(u'(b) WDM7_v2(more ccn)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.grid(linestyle='--')
plt.setp(ax.get_yticklabels(), visible=False)
plt.gca().invert_yaxis()

ax = fig.add_subplot(2, 2, 3)
ax.plot(era5_jan_lev[::-1], levels_i, color='black', label='ERA5', linewidth=1, linestyle='--')
plt.plot(regcm5_jan_lev3, levels_ii, color='red', label='RegCM5', linewidth=1)
plt.title(u'(c) WDM7_v3(less ccn)', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Level pressure (hPa)', fontsize=font_size, fontweight='bold')
plt.xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()

ax = fig.add_subplot(2, 2, 4)
ax.plot(era5_jan_lev[::-1], levels_i, color='black', label='ERA5', linewidth=1, linestyle='--')
plt.plot(regcm5_jan_lev4, levels_ii, color='red', label='RegCM5', linewidth=1)
plt.title(u'(d) WDM7_v4(less ccn2)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(dict_plot[var][0], fontsize=font_size, fontweight='bold')
plt.xlim(dict_plot[var][1], dict_plot[var][2])
plt.ylim(0,1000)
plt.xticks(dict_plot[var][3], fontsize=font_size)
plt.grid(linestyle='--')
plt.setp(ax.get_yticklabels(), visible=False)
plt.gca().invert_yaxis()

# Path out to save figure
path_out = '{0}/user/mdasilva/EUR-11/figs'.format(path)
name_out = 'pyplt_vertical_profile_{0}_{1}_RegCM5_WDM7_v1-v2-v3-v4_20000101.png'.format(var, domain)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
