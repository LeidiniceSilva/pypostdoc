# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 02, 2026"
__description__ = "This script plot trend maps"

import os
import numpy as np
import argparse
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from scipy.stats import t as student_t
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
if domain == 'EURR-3': 
    DOMAIN_EXTENT = [-25, 38, 33, 65]
elif domain == 'CAR-4':
    DOMAIN_EXTENT = [-119, -58, 9, 36]
else: 
    DOMAIN_EXTENT = [-79, -35, -37, -12]

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
        'era5': 'tas',
        'NorESM2-MM': 'tas',
        'MPI-ESM1-2-HR': 'tas',
        'EC-Earth3-Veg': 'tas',
        'units': '°C yr⁻¹'
    },
    'tasmax': {
        'cru': 'tmx',
        'era5': 'tasmax',
        'NorESM2-MM': 'tasmax',
        'MPI-ESM1-2-HR': 'tasmax',
        'EC-Earth3-Veg': 'tasmax',
        'units': '°C yr⁻¹'
    },
    'tasmin': {
        'cru': 'tmn',
        'era5': 'tasmin',
        'NorESM2-MM': 'tasmin',
        'MPI-ESM1-2-HR': 'tasmin',
        'EC-Earth3-Veg': 'tasmin',
        'units': '°C yr⁻¹'
    }
}

# Plot configuration
PR_LEVELS = { 
    'EURR-3': np.arange(-30,  33, 3),
    'CAR-4': np.arange(-40,  44, 4),
    'CSAM-3': np.arange(-60,  66, 6)
    }

PLOT_CONFIG = {
    'pr': {
        'title': 'Precipitation trend 2000-2009 (mm yr⁻¹) ',
        'levels': PR_LEVELS,
        'cmap': cm.BrBG,
        'sig_label': 'Dots: p < 0.05'
    },
    'tas': {
        'title': 'Air temperature trend 2000-2009 (°C yr⁻¹)',
        'levels': np.arange(-0.2, 0.21, 0.01),
        'cmap': cm.bwr,
        'sig_label': 'Dots: p < 0.05'
    },
    'tasmax': {
        'title': 'Maximum air temperature trend 2000-2009 (°C yr⁻¹)',
        'levels': np.arange(-0.2, 0.21, 0.01),
        'cmap': cm.bwr,
        'sig_label': 'Dots: p < 0.05'
    },
    'tasmin': {
        'title': 'Minimum air temperature trend 2000-2009 (°C yr⁻¹)',
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
    else:
        raise ValueError(f"Unknown dataset: {dataset}")

    filename = f"{path}/global/{var_name}_{dataset}_year_{dt}.nc"
    print(f"Reading {dataset}: {filename}")

    with netCDF4.Dataset(filename) as nc:
        data = nc.variables[var_name][:]
        lat = nc.variables[lat_name][:]
        lon = nc.variables[lon_name][:]

        # --- Fix longitude if needed ---
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
        
        if 'lat' in nc.variables and 'lon' in nc.variables:
            lat = nc.variables['lat'][:]
            lon = nc.variables['lon'][:]
            print("  2D lat/lon coordinates")
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
        
        return lat, lon, data


def calculate_trend(data, alpha=0.05, per_decade=False):

    nt, ny, nx = data.shape
    t = np.arange(nt, dtype=float)
    
    trend = np.full((ny, nx), np.nan)
    sig_mask = np.full((ny, nx), False)
   
    for i in range(ny):
        for j in range(nx):
            y = data[:, i, j]
            valid = np.isfinite(y)
            
            if np.sum(valid) >= 7:  
                t_valid = t[valid]
                y_valid = y[valid]
                
                # Linear regression
                A = np.vstack([t_valid, np.ones_like(t_valid)]).T
                slope, intercept = np.linalg.lstsq(A, y_valid, rcond=None)[0]
                
                # Calculate trend 
                trend_val = slope * 10.0 if per_decade else slope
                trend[i, j] = trend_val
                
                # Calculate p-value
                y_pred = slope * t_valid + intercept
                residuals = y_valid - y_pred
                ss_res = np.sum(residuals**2)
                
                if len(t_valid) > 2 and ss_res > 0:
                    ss_tot = np.sum((y_valid - np.mean(y_valid))**2)
                    r_squared = 1 - (ss_res / ss_tot)
                    df = len(t_valid) - 2
                    
                    # t-statistic
                    t_stat = slope / np.sqrt(ss_res / (df * np.sum((t_valid - np.mean(t_valid))**2)))
                    p_value = 2 * (1 - student_t.cdf(np.abs(t_stat), df))
                    
                    if p_value < alpha:
                        sig_mask[i, j] = True
    
    return trend, sig_mask


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
    print(f"  RegCM-ERA5 trend range: [{np.nanmin(regcm_data):.3f}, {np.nanmax(regcm_data):.3f}]")
    
    # Convert RegCM coordinates if needed
    print("\n2. Processing RegCM coordinates...")
    
    # Calculate trends
    print("\n3. Calculating trends...")
    cru_trend, cru_sig = calculate_trend(cru_data, alpha=0.05)
    print(f"  CRU trend range: [{np.nanmin(cru_trend):.3f}, {np.nanmax(cru_trend):.3f}]")
    era5_trend, era5_sig = calculate_trend(era5_data, alpha=0.05)
    print(f"  ERA5 trend range: [{np.nanmin(era5_trend):.3f}, {np.nanmax(era5_trend):.3f}]")
    regcm_trend, regcm_sig = calculate_trend(regcm_data, alpha=0.05)
    print(f"  RegCM-ERA5 trend range: [{np.nanmin(regcm_trend):.3f}, {np.nanmax(regcm_trend):.3f}]")
          
    # Create figure
    print("\n4. Creating plot...")
    fig, axes = plt.subplots(1, 3, figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree()})

    datasets = [('CRU', lon_cru, lat_cru, cru_trend, cru_sig),
    ('ERA5', lon_era5, lat_era5, era5_trend, era5_sig),
    ('RegCM5–ERA5', lon_regcm, lat_regcm, regcm_trend, regcm_sig)]

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
            
            ax.contourf(lon, lat, sig_binary, 
                       levels=[0.5, 1.5], 
                       hatches=['...'], 
                       colors='none', 
                       transform=ccrs.PlateCarree(),
                       alpha=0.0) 
        
        # Configure plot
        configure_plot(ax)
        label = ['(a)', '(b)', '(c)'][idx]
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

