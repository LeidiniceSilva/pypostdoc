# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot map of genesis density"

import os
import numpy as np
import xarray as xr
import scipy as sc
import scipy.ndimage
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeat

from scipy import signal, misc
from matplotlib.patches import Rectangle
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/TRACK-CYCLONE/Test-Mi/ERA5/S2R-Vortrack/SAM'


def import_dataset_1(dataset):

	df = pd.read_csv('{0}/postproc/genesis_density_{1}.txt'.format(path, dataset), sep="\s+", engine='python')
	df.columns = df.columns.str.strip()

	if 'data' not in df.columns:
		raise KeyError("Column 'data' not found in the dataset.")

	df['data'] = pd.to_datetime(df['data'], format='%Y%m%d%H',errors='coerce')
	df = df[(df['data'].dt.year >= 2000) & (df['data'].dt.year <= 2009)]

	return df


def import_dataset_2(dataset):

	df = pd.read_csv('{0}/postproc/genesis_density-track_{1}.txt'.format(path, dataset), sep="\s+", engine='python')
	df.columns = df.columns.str.strip()

	if 'data' not in df.columns:
		raise KeyError("Column 'data' not found in the dataset.")

	df['data'] = pd.to_datetime(df['data'], format='%Y%m%d%H',errors='coerce')
	df = df[(df['data'].dt.year >= 2000) & (df['data'].dt.year <= 2009)]

	return df

	
def density_ciclones(df):

	a = 6370 #raio da terra
	a2 = a * a
	lon1 = -93.0
	lon2 = 0.5
	lat1 = -60
	lat2 = -5
	dlon = 3
	dlat = 3
	lats = np.arange(lat1,lat2,dlat) #criando array #lat ficticia
	lons = np.arange(lon1,lon2,dlon) #criando array #lon ficticia
	nlat = lats.size
	nlon = lons.size
	entrada = np.array(df)
	nt=entrada[:,1].size
	rlats = lats * np.pi / 180
	drlon = dlon * np.pi / 180
	den = np.zeros((nlon,nlat))
	
	for t in range(nt):
		y0=entrada[t,1]
		y=np.absolute(lats-y0)
		latloc=np.argmin(y)
		x0 = entrada[t,2]
		x=np.absolute(lons-x0)
		lonloc=np.argmin(x)
		
		for j in range(nlat):
			for i in range(nlon):
				if(i==lonloc and j==latloc):
					den[i,j]= den[i,j] + 1
	
	aread = drlon * a2 * (np.sin(rlats[1:nlat-1]) - np.sin(rlats[0:nlat-2]))
                     
       
	for k in range(nlat-2):
		den[:,k] = (den[:,k]/aread[k])
     
	den = xr.DataArray(den, dims=['lon','lat'], coords={'lat': lats,'lon': lons})                   
        
	return den


def configure_subplot(ax):

	states_provinces = cfeat.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')

	ax.set_extent([-88, -31, -55, -20.5], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(-88,-31,5), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(-55,-20.5,5), crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.add_feature(cfeat.BORDERS)
	ax.add_feature(states_provinces, edgecolor='0.25')
	ax.coastlines()


# Import model and obs dataset    
era5_dateset_dens = import_dataset_1('ERA5')
dens_era5 = density_ciclones(era5_dateset_dens)

era5_dateset_track = import_dataset_2('ERA5')
track_era5 = density_ciclones(era5_dateset_track)

# Plot figure
fig, axes = plt.subplots(1,1, figsize=(10, 10), subplot_kw={"projection": ccrs.PlateCarree()})
ax1 = axes

font_size = 10
colorb_1 = [10,20,30,40,50]
colorb = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200]

c1 = ax1.contour(scipy.ndimage.zoom(dens_era5.lon,3), scipy.ndimage.zoom(dens_era5.lat,3), scipy.ndimage.zoom(dens_era5,3).T/10*1e6, colorb_1, linewidths=1.25, colors="cyan", transform=ccrs.PlateCarree())
ax1.clabel(c1, colors='red', inline=True, inline_spacing=1, fontsize=8)
cf1 = ax1.contourf(scipy.ndimage.zoom(track_era5.lon,3), scipy.ndimage.zoom(track_era5.lat,3), scipy.ndimage.zoom(track_era5,3).T/10*1e6, colorb, transform=ccrs.PlateCarree(), extend='both', cmap='gnuplot2_r')
ax1.set_title('(a) ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_ylabel('Latitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

cb = plt.colorbar(cf1, ticks=colorb)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_maps_density-tracking_ERA5_SAM_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
