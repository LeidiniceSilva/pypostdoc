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

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/TRACK-CYCLONE/CORDEX-TF'


def import_dataset(dataset, domain):

	df = pd.read_csv('{0}/{1}/S2R-Vortrack/{2}/genesis_density_{1}_{2}.txt'.format(path, dataset, domain), sep="\s+", engine='python')
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
	lat1 = -56
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

	ax.set_extent([-80, -30, -56, -16], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(-80,-30,10), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(-56,-16,5), crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.add_feature(cfeat.BORDERS)
	ax.add_feature(states_provinces, edgecolor='0.25')
	ax.coastlines()


# Import model and obs dataset    
dateset_obs = import_dataset('ERA5', 'SAM')
dateset_gcm1 = import_dataset('EC-Earth3-Veg', 'SAM')
dateset_gcm2 = import_dataset('MPI-ESM1-2-HR', 'SAM')
dateset_gcm3 = import_dataset('CNRM-ESM2-1', 'SAM')

dens_obs = density_ciclones(dateset_obs)
dens_gcm1 = density_ciclones(dateset_gcm1)
dens_gcm2 = density_ciclones(dateset_gcm2)
dens_gcm3 = density_ciclones(dateset_gcm3)

# Plot figure
fig, axes = plt.subplots(2,2, figsize=(8, 6), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2), (ax3, ax4) = axes

font_size = 10
colorb = [1,5,10,20,30,40,50,60,70,80,90,100]

cf1 = ax1.contourf(scipy.ndimage.zoom(dens_obs.lon,3), scipy.ndimage.zoom(dens_obs.lat,3), scipy.ndimage.zoom(dens_obs,3).T/8*1e6, colorb, transform=ccrs.PlateCarree(), extend='both', cmap='gnuplot2_r')
ax1.set_title('(a)', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_ylabel('Latitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

cf2 = ax2.contourf(scipy.ndimage.zoom(dens_gcm1.lon,3), scipy.ndimage.zoom(dens_gcm1.lat,3), scipy.ndimage.zoom(dens_gcm1,3).T/8*1e6, colorb, transform=ccrs.PlateCarree(), extend='both', cmap='gnuplot2_r')
ax2.set_title('(b)', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

cf3 = ax3.contourf(scipy.ndimage.zoom(dens_gcm2.lon,3), scipy.ndimage.zoom(dens_gcm2.lat,3), scipy.ndimage.zoom(dens_gcm2,3).T/8*1e6, colorb, transform=ccrs.PlateCarree(), extend='both', cmap='gnuplot2_r')
ax3.set_title('(c)', loc='left', fontsize=font_size, fontweight='bold')
ax3.set_xlabel('Longitude', fontsize=font_size, fontweight='bold')
ax3.set_ylabel('Latitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

cf4 = ax4.contourf(scipy.ndimage.zoom(dens_gcm3.lon,3), scipy.ndimage.zoom(dens_gcm3.lat,3), scipy.ndimage.zoom(dens_gcm3,3).T/8*1e6, colorb, transform=ccrs.PlateCarree(), extend='both', cmap='gnuplot2_r')
ax4.set_title('(d)', loc='left', fontsize=font_size, fontweight='bold')
ax4.set_xlabel('Longitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax4)

cb = plt.colorbar(cf3, ticks=colorb, cax=fig.add_axes([0.92, 0.3, 0.015, 0.4]))

# Path out to save figure
path_out = '{0}/figs/S2R-Vortrack'.format(path)
name_out = 'pyplt_maps_density-tracking_SAM_2000-2009.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

