# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Sep 22, 2025"
__description__ = "This script plots composite around cyclone center"


import os
import glob
import netCDF4
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeat
from dict_inmet_stations import inmet

from matplotlib.patches import Circle
from scipy.interpolate import griddata
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/TRACK-CYCLONE/CPRCM-TF'
font_size = 10
day=29


def import_ws():
	
	yy, xx, mean = [], [], [] 
	for station in inmet.keys():
		print(station, inmet[station][0])

		file_path = '{0}/data/inmet/dados_{1}_D_2009-01-01_2009-12-31.csv'.format(path, inmet[station][0])
		if not os.path.exists(file_path):
                       continue

		df = pd.read_csv(file_path, skiprows=11, encoding="ISO-8859-1", decimal=",", delimiter=';', header=None, usecols=[0,1], names=["Data", "Precipitacao"])
		df["Precipitacao"] = pd.to_numeric(df["Precipitacao"], errors='coerce')
		row = df[df["Data"] == '2009-01-29']

		if row.empty:
			continue

		yy.append(inmet[station][2])
		xx.append(inmet[station][3])
		mean.append(float(row["Precipitacao"].iloc[0]))
		
	return yy, xx, mean


def import_data(param, dataset):
	
	arq = '{0}/data/Jan2009/{1}_Jan20-29_{2}_lonlat.nc'.format(path, param, dataset)
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:]    
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]

	if dataset == 'ERA5':
		mean = var[9, :, :] * 1000.0
	elif dataset == 'RegCM5':
		mean = var[9, :, :] * 86400
	else:
		mean = var[9, :, :] 
							
	return lat, lon, mean

# Import data
yy_ws, xx_ws, ws = import_ws()
lat, lon, pr_era5 = import_data('tp', 'ERA5')
lat, lon, pr_cmorph = import_data('cmorph', 'CMORPH')
lat, lon, pr_regcm5 = import_data('pr', 'RegCM5')
lat, lon, pr_wrf415 = import_data('PREC_ACC_NC', 'WRF415')

print(pr_era5.shape)
print(np.min(pr_era5), np.max(pr_era5))
print()
print(pr_cmorph.shape)
print(np.min(pr_cmorph), np.max(pr_cmorph))
print()
print(pr_regcm5.shape)
print(np.min(pr_regcm5), np.max(pr_regcm5))
print()
print(pr_wrf415.shape)
print(np.min(pr_wrf415), np.max(pr_wrf415))

def configure_subplot(ax):

	ax.set_extent([-66, -38, -36, -14], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(-66,-38,4), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(-36,-14,4), crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.set_xlabel('Longitude', fontsize=font_size, fontweight='bold')
	ax.set_ylabel('Latitude',  fontsize=font_size, fontweight='bold')
	ax.grid(c='k', ls='--', alpha=0.5)
	ax.add_feature(cfeat.BORDERS, linewidth=2, color='gray')
	ax.coastlines(linewidth=2, color='gray')	


# Plot figure
fig, axes = plt.subplots(2,2, figsize=(12, 10), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2), (ax3, ax4) = axes
levels_pr = np.arange(0, 105, 5)

cf1 = ax1.contourf(lon, lat, pr_era5, levels=levels_pr, cmap='Blues', extend='max')
sc1 = ax1.scatter(xx_ws, yy_ws, 12, ws, cmap='Blues', edgecolors='black', linewidth=0.5, marker='o', vmin=0, vmax=75) 
ax1.set_title(u'(a) ERA5+INMET', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

cf2 = ax2.contourf(lon, lat, pr_cmorph, levels=levels_pr, cmap='Blues', extend='max')
ax2.set_title(u'(b) CMORPH', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

cf3 = ax3.contourf(lon, lat, pr_regcm5, levels=levels_pr, cmap='Blues', extend='max')
ax3.set_title(u'(c) RegCM5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

cf4 = ax4.contourf(lon, lat, pr_wrf415, levels=levels_pr, cmap='Blues', extend='max')
ax4.set_title(u'(d) WRF415', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax4)

cbar = plt.colorbar(cf4, cax=fig.add_axes([0.91, 0.33, 0.015, 0.33]))
cbar.set_label('Precipitation (mm/d)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Save figura
path_out = '{0}/figs/'.format(path)
name_out = 'pyplt_maps_Subtropical_cyclone_29Jan2009.png'.format(day)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


