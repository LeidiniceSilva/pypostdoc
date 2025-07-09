
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

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/SAM-3km'


def import_dataset(dataset):

	df = pd.read_csv('{0}/postproc/cyclone/{1}/genesis_density_{1}.txt'.format(path, dataset), sep="\s+", engine='python')
	df.columns = df.columns.str.strip()

	if 'data' not in df.columns:
		raise KeyError("Column 'data' not found in the dataset.")

	df['data'] = pd.to_datetime(df['data'], format='%Y%m%d%H',errors='coerce')
	df = df[(df['data'].dt.year >= 2018) & (df['data'].dt.year <= 2021)]

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

	ax.set_extent([-76, -38.5, -34.5, -15], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(-76,-38.5,5), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(-34.5,-15,5), crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.add_feature(cfeat.BORDERS)
	ax.add_feature(states_provinces, edgecolor='0.25')
	ax.coastlines()


# Import model and obs dataset    
era5_dateset = import_dataset('ERA5')
regcm5_dateset = import_dataset('RegCM5')
wrf415_dateset = import_dataset('WRF415')

dens_era5 = density_ciclones(era5_dateset)
dens_regcm5 = density_ciclones(regcm5_dateset)
dens_wrf415 = density_ciclones(wrf415_dateset)

# Plot figure
fig, axes = plt.subplots(2,2, figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2), (ax3, ax4) = axes
fig.delaxes(ax4)

font_size = 10
colorb = [1,3,5,7,9,12,15,20]

cf1 = ax1.contourf(scipy.ndimage.zoom(dens_era5.lon,3), scipy.ndimage.zoom(dens_era5.lat,3), scipy.ndimage.zoom(dens_era5,3).T/5*1e6, colorb, transform=ccrs.PlateCarree(), extend='both', cmap='gnuplot2_r')
ax1.set_title('(a) ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_ylabel('Latitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

cf2 = ax2.contourf(scipy.ndimage.zoom(dens_regcm5.lon,3), scipy.ndimage.zoom(dens_regcm5.lat,3), scipy.ndimage.zoom(dens_regcm5,3).T/5*1e6, colorb, transform=ccrs.PlateCarree(), extend='both', cmap='gnuplot2_r')
ax2.set_title('(b) RegCM5', loc='left', fontsize=font_size, fontweight='bold')
ax2.set_xlabel('Longitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

cf3 = ax3.contourf(scipy.ndimage.zoom(dens_wrf415.lon,3), scipy.ndimage.zoom(dens_wrf415.lat,3), scipy.ndimage.zoom(dens_wrf415,3).T/5*1e6, colorb, transform=ccrs.PlateCarree(), extend='both', cmap='gnuplot2_r')
ax3.set_title('(c) WRF415', loc='left', fontsize=font_size, fontweight='bold')
ax3.set_xlabel('Longitude', fontsize=font_size, fontweight='bold')
ax3.set_ylabel('Latitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

cb = plt.colorbar(cf3, ticks=colorb, cax=fig.add_axes([0.92, 0.3, 0.015, 0.4]))

# Path out to save figure
path_out = '{0}/figs/cyclone'.format(path)
name_out = 'pyplt_maps_density_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

