# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot map of tracking density"

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
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

font_size = 10
dataset = 'ERA5'
path = '/marconi/home/userexternal/mdasilva/user/mdasilva/SAM-3km'


def density_ciclones(df):

    a = 6370 #raio da terra
    a2 = a * a
    lon1 = -76    
    lon2 = -38.5   
    lat1 = -34.5  
    lat2 = -15      
    dlon = 3 
    dlat = 3
    lats = np.arange(lat1,lat2,dlat) #criando array #lat ficticia
    lons = np.arange(lon1,lon2,dlon) #criando array #lon ficticia
    nlat = lats.size
    nlon = lons.size
    entrada = np.array(df)   
    nt=entrada[:,1].size
    
    rlats = lats * np.pi / 180         #latitudes em radianos
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
     
    den = xr.DataArray(den, dims=['lon','lat'],
                            coords={'lat': lats,'lon': lons})                   
        
    return den


# Import model and obs dataset    
df = pd.read_csv('{0}/post_cyclone/ECyclone/ECyclone_{1}/tracking_density_{1}.txt'.format(path, dataset), sep=" ")
df['data'] = pd.to_datetime(df['data'], format='%Y%m%d%H',errors='coerce')

df = df[(df['data'].dt.year >= 2018) & (df['data'].dt.year <= 2021)]
dens = density_ciclones(df)

# Plot figure
fig, ax = plt.subplots(figsize=(6,4), subplot_kw={'projection': ccrs.PlateCarree()})

xticks = np.arange(-76,-38.5,5)
yticks = np.arange(-34.5,-15,5)

colorb = [10,20,30,40,60,80,100,120,140]
levels = [10,20,30,40,60,80,100,120,140]

#[.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.2,1.4,1.6,1.8,2,2.2,2.4]

ax.coastlines()
ax.set_xticks(xticks, crs=ccrs.PlateCarree())
ax.set_yticks(yticks, crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
ax.grid(c='k', ls='--', alpha=0.3)
ax.add_feature(cfeat.BORDERS)

states_provinces = cfeat.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')

ax.add_feature(states_provinces, edgecolor='0.25')
ax.set_xlabel('Longitude',fontsize=font_size, fontweight='bold')
ax.set_ylabel('Latitude',fontsize=font_size, fontweight='bold')
ax.set_title('b) {0}'.format(dataset), loc='left', fontsize=font_size, fontweight='bold')

cf = ax.contourf(scipy.ndimage.zoom(dens.lon,3),scipy.ndimage.zoom(dens.lat,3),scipy.ndimage.zoom(dens,3).T/4*1e2,colorb,transform=ccrs.PlateCarree(),extend='max',cmap='rainbow')
cb = plt.colorbar(cf,ticks=colorb,shrink=0.90,pad=0.07)

#cf1 = ax.contour(scipy.ndimage.zoom(dens.lon,3),scipy.ndimage.zoom(dens.lat,3),scipy.ndimage.zoom(dens,3).T/4*1e6,[2.5,5.5,8.5,11.5,14.5],linewidths=1.5, linestyles = '-',colors = 'turquoise')
#plt.clabel(cf1,[2.5,5.5,8.5,11.5,14.5], colors='red', fmt='%3i', fontsize=font_size)

plt.xticks(visible=True,fontsize=font_size)
plt.yticks(visible=True,fontsize=font_size)

# Path out to save figure
path_out = '{0}/figs/cyclone'.format(path)
name_out = 'pyplt_genesis_density_EC_{0}_SAM-3km_2018-2021.png'.format(dataset)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
