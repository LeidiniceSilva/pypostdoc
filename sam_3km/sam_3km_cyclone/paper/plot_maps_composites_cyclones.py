# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot composites"


import os,sys
import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.crs as ccrs
import metpy as m
import metpy.calc as mpcalc
import pandas as pd

from metpy.units import units
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

freq='day'
path='/marconi/home/userexternal/mdasilva'


def import_obs(param, dataset):

	arq = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_cyclone/ECyclone/ECyclone_ERA5/'.format(path) + 'composite_cyclones_ERA5_msl_2018-2021.nc'.format(param), engine='netcdf4', decode_times=False)
	idx = arq.variables['time']
	lon = arq.variables['longitude']
	lat = arq.variables['latitude']
	var = arq.variables[param] / 100
	nc_m, ny_m, nx_m = var.shape
	
	return var, nc_m, ny_m, nx_m

    	
def radial_line(circle, lonc, latc, radius, nsec, a):

	angle = int(360/nsec)
	xyTicker = circle.vertices[::angle]
	
	for i in range(0, xyTicker.shape[0]):
		ax[a].plot((xyTicker[i][0],lonc), (xyTicker[i][1],latc), linewidth='0.5', linestyle='--',color='darkgray')


def get_mask(x, y, circle):

	yy, xx = np.meshgrid(y, x)
	path = mpath.Path(circle.vertices)
	points = np.transpose((xx.ravel(), yy.ravel()))
	mask = path.contains_points( points )
	mask = mask.reshape(xx.shape)
	mask = np.invert(mask)
	
	return mask

	
# Import model and obs dataset 
era5_, time_, lat_, lon_ = import_obs('msl', 'ERA5')

lat = '{0}'.format(lat_)
lon = '{0}'.format(lon_)

# Calculate mean
mean_msl = np.nanmean(era5_, 0)
	
# Radial plot
latc = lat[49] + 0.125
lonc = lon[49] + 0.125

theta = np.linspace(0, 2*np.pi, 360)
center, radius = [lonc, latc], 20.
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

rMask = get_mask(lon, lat, circle)
mean_msl[rMask] = np.nan
    
# Plot figure
fig, ax = plt.subplots(nrows=1,ncols=1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(4,4))
ax=ax.flatten() # axs is a 2 dimensional array of `GeoAxes`.  We will flatten it into a 1-D array

sl=10
mpl.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
mpl.rc('xtick', labelsize=sl); mpl.rc('ytick', labelsize=sl)

clev = np.arange(0, 14, 1)
clevARG = [2,4,6,8,10,12,14,16,17,18]
clevURU = [2,4,6,8,10,15,20,25,30,35,40]
clevSBR = [2,4,6,8,10,12,14,16,18,20,25,30,35]
cleva = [2,4,6,8,10,15,20,25,30,35,40]
			 
llev = np.arange(np.round(np.nanmin(mean_msl)-2), np.round(np.nanmax(mean_msl)+2), 2)
llev1 = np.arange(np.round(np.nanmin(mean_msl)-2), np.round(np.nanmax(mean_msl)+2), 2)

cs = ax.contour(lon, lat, mean_msl, llev, colors='k', linewidths=0.5, linestyles='-')
ax.set_title('a) ERA5', fontsize=14)
plt.clabel(cs1, llev1,fontsize=7, inline=1, inline_spacing=4, fmt='%i', use_clabeltext=True,colors='r')

# Radial grid and axes
patch = mpatches.PathPatch(circle, facecolor='none')
ax.add_patch(patch)    
ax.axis('off')
radial_line(circle,lonc, latc, radius,6,i) 

#Adjust the location of the subplots on the page to make room for the colorbar
fig.subplots_adjust(bottom=0.20, top=0.9, left=0.01, right=0.99, wspace=0.08, hspace=0.2)

# Add a colorbar axis at the bottom of the graph
cbar_ax = fig.add_axes([0.095, 0.2, 0.80, 0.055])

# Draw the colorbar
cbar=plt.colorbar(cs, ticks=cleva, cax=cbar_ax, orientation='horizontal', pad=0.1)
cbar.ax.tick_params(labelsize=14) 

# Path out to save figure
path_out = '{0}/user/mdasilva/SAM-3km/figs/cyclone/paper'.format(path)
name_out = 'pyplt_maps_composites_radial_mslp_CP-RCM_SAM-3km_2018-2021.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
