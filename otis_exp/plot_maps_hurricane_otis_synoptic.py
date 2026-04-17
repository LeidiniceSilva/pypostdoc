# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "May 14, 2025"
__description__ = "This script plot study area"

import os
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Paths
path = "/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/era5/"

# Domain
lon_min, lon_max = -114, -80
lat_min, lat_max = 2, 24


def load_var(file, varname, time_sel=None):
    ds = xr.open_dataset(path + file)

    # corrigir longitude
    if ds.longitude.max() > 180:
        ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180)).sortby("longitude")

    # selecionar tempo (se fornecido)
    if time_sel is not None:
        ds = ds.sel(valid_time=time_sel)

    return ds[varname], ds["longitude"], ds["latitude"]


def setup(ax):

    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.set_xticks(np.arange(lon_min, lon_max, 4), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(lat_min, lat_max, 4), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(font_size)
   
    ax.grid(c='gray', ls='--')
    ax.coastlines(linewidth=0.75)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)


# Import vars
mslp, lon, lat = load_var("msl_ERA5_Otis_6hr_23-24Oct2023.nc", "msl", "2023-10-23T06:00")
sst, _, _      = load_var("sst_ERA5_Otis_6hr_23-24Oct2023.nc", "sst", "2023-10-23T06:00")
u850, _, _     = load_var("u_850hPa_ERA5_Otis_6hr_23-24Oct2023.nc", "u", "2023-10-23T06:00")
v850, _, _    = load_var("v_850hPa_ERA5_Otis_6hr_23-24Oct2023.nc", "v", "2023-10-23T06:00")
u200, _, _     = load_var("u_200hPa_ERA5_Otis_6hr_23-24Oct2023.nc", "u", "2023-10-23T06:00")
v200, _, _     = load_var("v_200hPa_ERA5_Otis_6hr_23-24Oct2023.nc", "v", "2023-10-23T06:00")
u10, _, _      = load_var("u10_ERA5_Otis_6hr_23-24Oct2023.nc", "u10", "2023-10-23T06:00")
v10, _, _      = load_var("v10_ERA5_Otis_6hr_23-24Oct2023.nc", "v10", "2023-10-23T06:00")
q, _, _        = load_var("q_850hPa_ERA5_Otis_6hr_23-24Oct2023.nc", "q", "2023-10-23T06:00")
tp, _, _       = load_var("tp_ERA5_Otis_day_Oct2023.nc", "tp")

mslp = mslp / 100.0

du = u200[0] - u850[0]
dv = v200[0] - v850[0]
shear = np.sqrt(du**2 + dv**2)

wspd = np.sqrt(u10**2 + v10**2)

q = q[0] * 1000.0

tp = tp.sel(valid_time=slice("2023-10-24", "2023-10-26"))
tp = tp.sum(dim="valid_time") * 1000.0  

# Plot figure
plt.figure(figsize=(10, 6))
font_size = 8
step = 5

sst_levels   = np.arange(298, 306, 0.5)  # K
mslp_levels  = np.arange(998, 1018, 1)   # hPa
q_levels     = np.arange(0, 18, 1)       # g/kg
shear_levels = np.arange(0, 26, 1)       # m/s
prec_levels  = np.arange(0, 320, 20)     # mm/h

# MSLP
plt.subplot(2, 2, 1, projection=ccrs.PlateCarree())
ax = plt.gca()
setup(ax)
cf = plt.contourf(lon, lat, sst, levels=sst_levels, cmap='bwr', transform=ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=10)
ct = plt.contour(lon, lat, mslp, levels=np.arange(998, 1015, 0.5), linewidths=0.5, colors='black', transform=ccrs.PlateCarree())
plt.title("(a)", loc='left', fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(cf, ax=ax, pad=0.01, fraction=0.03)
cbar.set_label("SST (K)", fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Wind shear
plt.subplot(2, 2, 2, projection=ccrs.PlateCarree())
ax = plt.gca()
setup(ax)
cf = plt.contourf(lon, lat, shear, levels=shear_levels, cmap='gist_ncar_r', transform=ccrs.PlateCarree())
qv = plt.quiver(lon[::step], lat[::step], du[::step, ::step], dv[::step, ::step], scale=200, width=0.002, transform=ccrs.PlateCarree(), color='gray')
plt.title("(b)", loc='left', fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(cf, ax=ax, pad=0.01, fraction=0.03)
cbar.set_label("Wind shear (m s⁻¹)", fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Specific humidity 
plt.subplot(2, 2, 3, projection=ccrs.PlateCarree())
ax = plt.gca()
setup(ax)
cf = plt.contourf(lon, lat, q, levels=q_levels, cmap='rainbow_r', transform=ccrs.PlateCarree())
st = plt.streamplot(lon, lat, u10[:,:], v10[:,:], arrowsize=1, arrowstyle='->', color='white', density=1, linewidth=0.5)
plt.title("(c)", loc='left', fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(cf, ax=ax, pad=0.01, fraction=0.03)
cbar.set_label("Q (g kg⁻¹)", fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Total precipitation 
plt.subplot(2, 2, 4, projection=ccrs.PlateCarree())
ax = plt.gca()
setup(ax)
cf = plt.contourf(lon, lat, tp, levels=prec_levels, cmap='Blues', transform=ccrs.PlateCarree())
ct = plt.contour(lon, lat, tp, levels=prec_levels, linewidths=0.5, colors='gray', transform=ccrs.PlateCarree())
plt.quiverkey(qv, 0.85, 1.05, 10, '10 m/s', labelpos="E")
plt.title("(d)", loc='left', fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(cf, ax=ax, pad=0.01, fraction=0.03)
cbar.set_label("Precipitation (mm d⁻¹)", fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/figs/'
name_out = 'pyplt_Hurricane_Otis_synotic.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()



