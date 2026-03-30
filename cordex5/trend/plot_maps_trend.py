# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 02, 2026"
__description__ = "This script plot trend maps"

import os
import netCDF4
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeat

from function_trend import calculate_trend
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

parser = argparse.ArgumentParser(description='Plot trend maps')
parser.add_argument('--var', required=True, help='Variable name')
parser.add_argument('--domain', required=True, help='Domain name')
parser.add_argument('--idt', required=True, help='Initial year')
parser.add_argument('--fdt', required=True, help='Final year')
args = parser.parse_args()

# Configuration
var = args.var
domain = args.domain
idt = args.idt
fdt = args.fdt
dt = f'{idt}-{fdt}'
font_size = 10

path = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/postproc/trend'

# Domain extent
if domain == 'AUS-12': 
    DOMAIN_EXTENT = [110, 180, -48, 4]
elif domain == 'EUR-12': 
    DOMAIN_EXTENT = [-42, 61, 20, 70]
elif domain == 'NAM-12': 
    DOMAIN_EXTENT = [-171, -22, 11, 75]
else: 
    DOMAIN_EXTENT = [-105, -15, -57, 18]

# Variable name 
VAR_MAP = {
    'pr': {
        'cru': 'pre',
        'era5': 'tp',
        'NorESM2-MM': 'pr',
        'MPI-ESM1-2-HR': 'pr',
        'EC-Earth3-Veg': 'pr',
        'units': 'mm yr⁻¹'
    },
    'tas': {
        'cru': 'tmp',
        'era5': 't2m',
        'NorESM2-MM': 'tas',
        'MPI-ESM1-2-HR': 'tas',
        'EC-Earth3-Veg': 'tas',
        'units': '°C yr⁻¹'
    },
    'tasmax': {
        'cru': 'tmx',
        'era5': 'mx2t',
        'NorESM2-MM': 'tasmax',
        'MPI-ESM1-2-HR': 'tasmax',
        'EC-Earth3-Veg': 'tasmax',
        'units': '°C yr⁻¹'
    },
    'tasmin': {
        'cru': 'tmn',
        'era5': 'mn2t',
        'NorESM2-MM': 'tasmin',
        'MPI-ESM1-2-HR': 'tasmin',
        'EC-Earth3-Veg': 'tasmin',
        'units': '°C yr⁻¹'
    }
}

# Plot configuration
PR_LEVELS = { 
    'AUS-12': np.arange(-20,  22, 2),
    'EUR-12': np.arange(-10,  11, 1),
    'NAM-12': np.arange(-10,  11, 1),
    'SAM-12': np.arange(-20,  22, 2)
}

PLOT_CONFIG = {
    'pr': {
        'title': 'Precipitation trend 1970-2014 (mm yr⁻¹) ',
        'levels': PR_LEVELS,
        'cmap': cm.BrBG,
        'sig_label': 'Dots: p < 0.05'
    },
    'tas': {
        'title': 'Air temperature trend 1970-2014 (°C yr⁻¹)',
        'levels': np.arange(-0.2, 0.21, 0.01),
        'cmap': cm.bwr,
        'sig_label': 'Dots: p < 0.05'
    },
    'tasmax': {
        'title': 'Maximum air temperature trend 1970-2014 (°C yr⁻¹)',
        'levels': np.arange(-0.2, 0.21, 0.01),
        'cmap': cm.bwr,
        'sig_label': 'Dots: p < 0.05'
    },
    'tasmin': {
        'title': 'Minimum air temperature trend 1970-2014 (°C yr⁻¹)',
        'levels': np.arange(-0.2, 0.21, 0.01),
        'cmap': cm.bwr,
        'sig_label': 'Dots: p < 0.05'
    }
}


def import_regular_dataset(param, dataset):

    if dataset == 'CRU':
        var_name = VAR_MAP[param]['cru']
        lat_name, lon_name = 'lat', 'lon'
        lon_is_0360 = False
    elif dataset == 'ERA5':
        var_name = VAR_MAP[param]['era5']
        lat_name, lon_name = 'latitude', 'longitude'
        lon_is_0360 = True
    elif dataset == 'EC-Earth3-Veg':
        var_name = VAR_MAP[param]['EC-Earth3-Veg']
        lat_name, lon_name = 'lat', 'lon'
        lon_is_0360 = True
    elif dataset == 'MPI-ESM1-2-HR':
        var_name = VAR_MAP[param]['MPI-ESM1-2-HR']
        lat_name, lon_name = 'lat', 'lon'
        lon_is_0360 = True
    elif dataset == 'NorESM2-MM':
        var_name = VAR_MAP[param]['NorESM2-MM']
        lat_name, lon_name = 'lat', 'lon'
        lon_is_0360 = True
    else:
        raise ValueError(f"Unknown dataset: {dataset}")

    filename = f"{path}/global/{var_name}_{dataset}_year_{dt}.nc"
    print(f"Reading {dataset}: {filename}")

    with netCDF4.Dataset(filename) as nc:
        data = nc.variables[var_name][:]
        lat = nc.variables[lat_name][:]
        lon = nc.variables[lon_name][:]

        # Fix longitude if needed 
        if lon_is_0360 and lon.max() > 180:
            lon = ((lon + 180) % 360) - 180
            idx = np.argsort(lon)
            lon = lon[idx]

            # Assume data is (..., lat, lon)
            data = data[..., idx]

        if lat.ndim == 1 and lon.ndim == 1:
            lon, lat = np.meshgrid(lon, lat)

        print(f"  Data shape: {data.shape}")
        print(f"  Lon range: [{lon.min():.1f}, {lon.max():.1f}]")
        print(f"  Lat range: [{lat.min():.1f}, {lat.max():.1f}]")

    return lat, lon, data


def import_regcm_dataset(param, driven):

    dataset = f"{domain}_{driven}"
    filename = f"{path}/{domain}/{param}_{dataset}_year_{dt}.nc"
    print(f"Reading RegCM5: {filename}")
    
    with netCDF4.Dataset(filename) as nc:
        data = nc.variables[param][:]

        if data.ndim == 4 and data.shape[1] == 1:
            data = data[:, 0, :, :]
        
        if 'lat' in nc.variables and 'lon' in nc.variables:
            lat = nc.variables['lat'][:]
            lon = nc.variables['lon'][:]
            print("  2D lat/lon coordinates")
        elif 'xlat' in nc.variables and 'xlon' in nc.variables:
            lat = nc.variables['xlat'][:]
            lon = nc.variables['xlon'][:]
            print("  2D lat/lon coordinates")
        elif 'rlat' in nc.variables and 'rlon' in nc.variables:
            rlat = nc.variables['rlat'][:]
            rlon = nc.variables['rlon'][:]
            lon, lat = np.meshgrid(rlon, rlat)
            print("  Rotated coordinates (rlat/rlon)")
        else:
            lat_names = ['lat', 'latitude', 'y', 'nav_lat', 'XLAT']
            lon_names = ['lon', 'longitude', 'x', 'nav_lon', 'XLONG']
            
            lat_var = next((n for n in lat_names if n in nc.variables), None)
            lon_var = next((n for n in lon_names if n in nc.variables), None)
            
            if lat_var and lon_var:
                lat = nc.variables[lat_var][:]
                lon = nc.variables[lon_var][:]
                print(f"  Found coordinates: {lat_var}, {lon_var}")
            else:
                raise ValueError(f"Cannot find coordinates in {filename}")
        
        if len(lat.shape) == 2 and len(lon.shape) == 2:
            print("  2D coordinates detected")
        elif len(lat.shape) == 1 and len(lon.shape) == 1:
            lon, lat = np.meshgrid(lon, lat)
        
        print(f"  Data shape: {data.shape}")
        print(f"  Lon shape: {lon.shape}, range: [{lon.min():.1f}, {lon.max():.1f}]")
        print(f"  Lat shape: {lat.shape}, range: [{lat.min():.1f}, {lat.max():.1f}]")

        if domain == 'AUS-12':
            if lon.min() < 0:
               lon = lon % 360
        else:
            'Not needed'
        
        return lat, lon, data


def convert_rotated_to_regular(lon_rot, lat_rot):

    # Rotated pole parameters from your file
    rot_pole_lat = -60.31    # South pole latitude
    rot_pole_lon = 321.38    # South pole longitude
    
    # Convert to north pole convention for cartopy
    north_pole_lat = -rot_pole_lat         # 60.31°N
    north_pole_lon = rot_pole_lon - 180.0  # 141.38°E
    
    print(f"  Converting ({north_pole_lat:.2f}°N, {north_pole_lon:.2f}°E)")
    
    # Create coordinate systems
    rotated_crs = ccrs.RotatedPole(pole_latitude=north_pole_lat,
                                   pole_longitude=north_pole_lon)
    regular_crs = ccrs.PlateCarree()
    
    # Transform coordinates
    points = rotated_crs.transform_points(regular_crs, lon_rot, lat_rot)
    lon_reg = points[..., 0]
    lat_reg = points[..., 1]
    
    # Ensure longitude is in 0-360 range
    lon_reg = lon_reg % 360
    
    print(f"  Converted lon range: [{lon_reg.min():.1f}, {lon_reg.max():.1f}]")
    print(f"  Converted lat range: [{lat_reg.min():.1f}, {lat_reg.max():.1f}]")
    
    return lon_reg, lat_reg


def configure_plot(ax):

    ax.set_extent(DOMAIN_EXTENT, crs=ccrs.PlateCarree())
    ax.set_xticks(np.arange(DOMAIN_EXTENT[0], DOMAIN_EXTENT[1]+1, 20), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(DOMAIN_EXTENT[2], DOMAIN_EXTENT[3]+1, 10), crs=ccrs.PlateCarree())
    
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
   
    ax.add_feature(cfeat.COASTLINE, linewidth=0.6)
    ax.add_feature(cfeat.BORDERS, linewidth=0.5, linestyle='-', alpha=0.5)
    ax.add_feature(cfeat.LAND, facecolor='lightgray', alpha=0.2)
    ax.gridlines(draw_labels=False, linewidth=0.5, linestyle='--', color='gray', alpha=0.5)


def main():
    """Processing function."""
    print(f"Processing variable: {var}")
    print(f"Time period: {dt}")
    print(f"Domain: {domain}")
    print("-" * 50)
    
    # Import data
    print("\n1. Importing data...")
    lat_cru,   lon_cru,   cru_data   = import_regular_dataset(var, 'CRU')
    lat_era5,  lon_era5,  era5_data  = import_regular_dataset(var, 'ERA5')
    lat_regcm, lon_regcm, regcm_data = import_regcm_dataset(var, 'ERA5_RegCM5')
    lat_rcm1, lon_rcm1, rcm1_data = import_regular_dataset(var, 'EC-Earth3-Veg')
    lat_rcm2, lon_rcm2, rcm2_data = import_regular_dataset(var, 'MPI-ESM1-2-HR')
    lat_rcm3, lon_rcm3, rcm3_data = import_regular_dataset(var, 'NorESM2-MM')
    lat_gcm1, lon_gcm1, gcm1_data = import_regular_dataset(var, 'EC-Earth3-Veg')
    lat_gcm2, lon_gcm2, gcm2_data = import_regular_dataset(var, 'MPI-ESM1-2-HR')
    lat_gcm3, lon_gcm3, gcm3_data = import_regular_dataset(var, 'NorESM2-MM')
    
    # Convert RegCM coordinates if needed
    print("\n2. Processing RegCM coordinates...")
    
    if domain == 'AUS-12':
        if (lon_regcm.min() < -100 or lon_regcm.max() > 300 or 
        lat_regcm.min() < -100 or lat_regcm.max() > 100):
            print("  Detected rotated coordinates, converting to regular lat/lon...")
            lon_regcm, lat_regcm = convert_rotated_to_regular(lon_regcm, lat_regcm)
        else:
            print("  Regular coordinates detected")
    
	# Mask data outside Australia domain
        mask = ((lon_regcm >= DOMAIN_EXTENT[0]) & (lon_regcm <= DOMAIN_EXTENT[1]) &
                (lat_regcm >= DOMAIN_EXTENT[2]) & (lat_regcm <= DOMAIN_EXTENT[3]))
        regcm_data = np.where(mask[None, :, :], regcm_data, np.nan)
    else:
        print("All the other domains uses regular lat/lon coordinates")
    
    # Calculate trends
    print("\n3. Calculating trends...")
    cru_trend, cru_sig = calculate_trend(cru_data, alpha=0.05)
    print(f"  CRU trend range: [{np.nanmin(cru_trend):.3f}, {np.nanmax(cru_trend):.3f}]")
    era5_trend, era5_sig = calculate_trend(era5_data, alpha=0.05)
    print(f"  ERA5 trend range: [{np.nanmin(era5_trend):.3f}, {np.nanmax(era5_trend):.3f}]")
    regcm_trend, regcm_sig = calculate_trend(regcm_data, alpha=0.05)
    print(f"  RegCM-ERA5 trend range: [{np.nanmin(regcm_trend):.3f}, {np.nanmax(regcm_trend):.3f}]")
    regcm_gcm1_trend, regcm_gcm1_sig = calculate_trend(rcm1_data, alpha=0.05)
    print(f"  RegCM-EC-Earth3-Veg trend range: [{np.nanmin(regcm_gcm1_trend):.3f}, {np.nanmax(regcm_gcm1_trend):.3f}]")
    regcm_gcm2_trend, regcm_gcm2_sig = calculate_trend(rcm2_data, alpha=0.05)
    print(f"  RegCM-MPI-ESM1-2-HR trend range: [{np.nanmin(regcm_gcm2_trend):.3f}, {np.nanmax(regcm_gcm2_trend):.3f}]")
    regcm_gcm3_trend, regcm_gcm3_sig = calculate_trend(rcm3_data, alpha=0.05)
    print(f"  RegCM-NorESM2-MM trend range: [{np.nanmin(regcm_gcm3_trend):.3f}, {np.nanmax(regcm_gcm3_trend):.3f}]")
    gcm1_trend, gcm1_sig = calculate_trend(gcm1_data, alpha=0.05)
    print(f"  EC-Earth3-Veg trend range: [{np.nanmin(gcm1_trend):.3f}, {np.nanmax(gcm1_trend):.3f}]")
    gcm2_trend, gcm2_sig = calculate_trend(gcm2_data, alpha=0.05)
    print(f"  MPI-ESM1-2-HR trend range: [{np.nanmin(gcm2_trend):.3f}, {np.nanmax(gcm2_trend):.3f}]")
    gcm3_trend, gcm3_sig = calculate_trend(gcm3_data, alpha=0.05)
    print(f"  NorESM2-MM trend range: [{np.nanmin(gcm3_trend):.3f}, {np.nanmax(gcm3_trend):.3f}]")
      
    # Create figure
    print("\n4. Creating plot...")

    if domain == 'AUS-12': 
        fig, axes = plt.subplots(3, 3, figsize=(12, 13), subplot_kw={'projection': ccrs.PlateCarree()})
    elif domain == 'EUR-12': 
        fig, axes = plt.subplots(3, 3, figsize=(14, 10), subplot_kw={'projection': ccrs.PlateCarree()})
    elif domain == 'NAM-12': 
        fig, axes = plt.subplots(3, 3, figsize=(15, 10), subplot_kw={'projection': ccrs.PlateCarree()})
    else: 
        fig, axes = plt.subplots(3, 3, figsize=(14, 12), subplot_kw={'projection': ccrs.PlateCarree()})
    
    datasets = [('CRU', lon_cru, lat_cru, cru_trend, cru_sig),
    ('ERA5', lon_era5, lat_era5, era5_trend, era5_sig),
    ('RegCM5–ERA5', lon_regcm, lat_regcm, regcm_trend, regcm_sig),
    ('RegCM5–EC-Earth3-Veg', lon_rcm1, lat_rcm1, regcm_gcm1_trend, regcm_gcm1_sig),
    ('RegCM5–MPI-ESM1-2-HR', lon_rcm2, lat_rcm2, regcm_gcm2_trend, regcm_gcm2_sig),
    ('RegCM5–NorESM2-MM', lon_rcm3, lat_rcm3, regcm_gcm3_trend, regcm_gcm3_sig),
    ('EC-Earth3-Veg', lon_gcm1, lat_gcm1, gcm1_trend, gcm1_sig),
    ('MPI-ESM1-2-HR', lon_gcm2, lat_gcm2, gcm2_trend, gcm2_sig),
    ('NorESM2-MM', lon_gcm3, lat_gcm3, gcm3_trend, gcm3_sig)]

    axes = axes.flatten() 
    for idx, (ax, (name, lon, lat, trend, sig)) in enumerate(zip(axes, datasets)):

        levels_cfg = PLOT_CONFIG[var]['levels']
        if var == 'pr':
            levels = levels_cfg.get(domain)
        else:
            levels = levels_cfg

        # Plot trend
        im = ax.contourf(lon, lat, trend, levels=levels, cmap=PLOT_CONFIG[var]['cmap'], extend='both', transform=ccrs.PlateCarree())
        
        # Plot significance dots
        if np.any(sig):
            sig_binary = np.where(sig, 1.0, 0.0)
            ax.contourf(lon, lat, sig_binary, levels=[0.5, 1.5], hatches=['...'], colors='none', transform=ccrs.PlateCarree(), alpha=0.0) 
        
        # Configure plot
        configure_plot(ax)
        label = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)'][idx]
        ax.set_title(f'{label} {name}', loc='left', fontsize=font_size, fontweight='bold')
    
    cbar = fig.colorbar(im, ax=axes, orientation='horizontal', pad=0.07, aspect=40, shrink=0.8)
    cbar.set_label(f"{PLOT_CONFIG[var]['title']} {PLOT_CONFIG[var]['sig_label']}", fontsize=font_size, fontweight='bold')
    cbar.ax.tick_params(labelsize=font_size)
    
    # Save figure 
    path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/figs/trend'
    os.makedirs(path_out, exist_ok=True)
    output_file = f'{path_out}/pyplt_maps_trend_{var}_{domain}_RegCM5_{dt}.png'
    print(f"\n4. Saving figure: {output_file}")
    plt.savefig(output_file, dpi=400, bbox_inches='tight')
    print("Done!")
    plt.show()

if __name__ == '__main__':
    main()

