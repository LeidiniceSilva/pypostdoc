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

var = 'nc'
domain = 'EUR-11'
dt = '2000-2001'
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11'

def import_rcm(param, dataset, season):

	arq   = '{0}/postproc/rcm/{1}_{2}_FPS_{3}_{4}_{5}_lonlat.nc'.format(path, param, domain, dataset, season, dt)
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	value = var[:][:,:,:,:]
	mean = np.nanmean(np.nanmean(value, axis=2), axis=2)
	
	return mean
	
	
# Import model and obs dataset
ncc_wdm7_djf = import_rcm('ncc', 'WDM7-Europe', 'DJF')
ncc_wdm7_mam = import_rcm('ncc', 'WDM7-Europe', 'MAM')
ncc_wdm7_jja = import_rcm('ncc', 'WDM7-Europe', 'JJA')
ncc_wdm7_son = import_rcm('ncc', 'WDM7-Europe', 'SON')

ncn_wdm7_djf = import_rcm('ncn', 'WDM7-Europe', 'DJF')
ncn_wdm7_mam = import_rcm('ncn', 'WDM7-Europe', 'MAM')
ncn_wdm7_jja = import_rcm('ncn', 'WDM7-Europe', 'JJA')
ncn_wdm7_son = import_rcm('ncn', 'WDM7-Europe', 'SON')

ncr_wdm7_djf = import_rcm('ncr', 'WDM7-Europe', 'DJF')
ncr_wdm7_mam = import_rcm('ncr', 'WDM7-Europe', 'MAM')
ncr_wdm7_jja = import_rcm('ncr', 'WDM7-Europe', 'JJA')
ncr_wdm7_son = import_rcm('ncr', 'WDM7-Europe', 'SON')

# Plot figure  
fig = plt.figure(figsize=(12, 8))
font_size = 8

levels_ii = (1000,925,850,700,600,500,400,300,250,200,150,100)

ax = fig.add_subplot(3, 4, 1)
plt.plot(ncc_wdm7_djf[0], levels_ii, marker='o', markersize=4, linestyle='-', color='green', alpha=0.70, label='NCC', linewidth=1)
plt.title(u'(a) DJF', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Level pressure (hPa)', fontsize=font_size, fontweight='bold')
plt.ylim(0,1000)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()
plt.legend(loc=1, ncol=1, fontsize=font_size)

ax = fig.add_subplot(3, 4, 2)
plt.plot(ncc_wdm7_mam[0], levels_ii, marker='o', markersize=4, linestyle='-', color='green', alpha=0.70, label='NCC', linewidth=1)
plt.title(u'(b) MAM', loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(0,1000)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()
plt.legend(loc=1, ncol=1, fontsize=font_size)

ax = fig.add_subplot(3, 4, 3)
plt.plot(ncc_wdm7_jja[0], levels_ii, marker='o', markersize=4, linestyle='-', color='green', alpha=0.70, label='NCC', linewidth=1)
plt.title(u'(c) JJA', loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(0,1000)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()
plt.legend(loc=1, ncol=1, fontsize=font_size)

ax = fig.add_subplot(3, 4, 4)
plt.plot(ncc_wdm7_son[0], levels_ii, marker='o', markersize=4, linestyle='-', color='green', alpha=0.70, label='NCC', linewidth=1)
plt.title(u'(d) SON', loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(0,1000)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()
plt.legend(loc=1, ncol=1, fontsize=font_size)

ax = fig.add_subplot(3, 4, 5)
plt.plot(ncn_wdm7_djf[0], levels_ii, marker='^', markersize=4, linestyle='-', color='green', alpha=0.70, label='NCN', linewidth=1)
plt.title(u'(e) DJF', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Level pressure (hPa)', fontsize=font_size, fontweight='bold')
plt.ylim(0,1000)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()
plt.legend(loc=1, ncol=1, fontsize=font_size)

ax = fig.add_subplot(3, 4, 6)
plt.plot(ncn_wdm7_mam[0], levels_ii, marker='^', markersize=4, linestyle='-', color='green', alpha=0.70, label='NCN', linewidth=1)
plt.title(u'(f) MAM', loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(0,1000)
plt.yticks(fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()
plt.legend(loc=1, ncol=1, fontsize=font_size)

ax = fig.add_subplot(3, 4, 7)
plt.plot(ncn_wdm7_jja[0], levels_ii, marker='^', markersize=4, linestyle='-', color='green', alpha=0.70, label='NCN', linewidth=1)
plt.title(u'(g) JJA', loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(0,1000)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()
plt.legend(loc=1, ncol=1, fontsize=font_size)

ax = fig.add_subplot(3, 4, 8)
plt.plot(ncn_wdm7_son[0], levels_ii, marker='^', markersize=4, linestyle='-', color='green', alpha=0.70, label='NCN', linewidth=1)
plt.title(u'(h) SON', loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(0,1000)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()
plt.legend(loc=1, ncol=1, fontsize=font_size)

ax = fig.add_subplot(3, 4, 9)
plt.plot(ncr_wdm7_djf[0], levels_ii, marker='s', markersize=4, linestyle='-', color='green', alpha=0.70, label='NCR', linewidth=1)
plt.title(u'(i) DJF', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Level pressure (hPa)', fontsize=font_size, fontweight='bold')
plt.xlabel('Number concentration WDM7 (m$^-$$^3$)', fontsize=font_size, fontweight='bold')
plt.ylim(0,1000)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()
plt.legend(loc=1, ncol=1, fontsize=font_size)

ax = fig.add_subplot(3, 4, 10)
plt.plot(ncr_wdm7_mam[0], levels_ii, marker='s', markersize=4, linestyle='-', color='green', alpha=0.70, label='NCR', linewidth=1)
plt.title(u'(j) MAM', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('Number concentration WDM7 (m$^-$$^3$)', fontsize=font_size, fontweight='bold')
plt.ylim(0,1000)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()
plt.legend(loc=1, ncol=1, fontsize=font_size)

ax = fig.add_subplot(3, 4, 11)
plt.plot(ncr_wdm7_jja[0], levels_ii, marker='s', markersize=4, linestyle='-', color='green', alpha=0.70, label='NCR', linewidth=1)
plt.title(u'(k) JJA', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('Number concentration WDM7 (m$^-$$^3$)', fontsize=font_size, fontweight='bold')
plt.ylim(0,1000)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()
plt.legend(loc=1, ncol=1, fontsize=font_size)

ax = fig.add_subplot(3, 4, 12)
plt.plot(ncr_wdm7_son[0], levels_ii, marker='s', markersize=4, linestyle='-', color='green', alpha=0.70, label='NCR', linewidth=1)
plt.title(u'(l) SON', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('Number concentration WDM7 (m$^-$$^3$)', fontsize=font_size, fontweight='bold')
plt.ylim(0,1000)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.grid(linestyle='--')
plt.gca().invert_yaxis()
plt.legend(loc=1, ncol=1, fontsize=font_size)

# Path out to save figure
path_out = '{0}/figs/totc'.format(path)
name_out = 'pyplt_graph_vertical_profile_{0}_{1}_RegCM5_{2}.png'.format(var, domain, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
