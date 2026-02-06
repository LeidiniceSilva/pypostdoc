# -*- coding: utf-8 -*-
"""
Script to plot trend maps for CRU, ERA5, and RegCM5 data.
Handles RegCM5 rotated pole coordinate system for AUS domain only.
"""

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jan 25, 2026"
__description__ = "This script plot trend maps with rotated pole support for AUS domain"

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
dataset = f'{domain}_RegCM5'
font_size = 8
path = '/leonardo_scratch/large/userexternal/fraffael/plots/trend/inputs'
#path = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/postproc/trend'

# Domain extent
if domain == 'AUS-12':
    DOMAIN_EXTENT = [110, 156, -45, -9]
elif domain == 'CSAM-3':
    DOMAIN_EXTENT = [-79, -35, -35, -10]
else:
    raise ValueError(f"Unknown domain: {domain}")

# Variable mapping
VAR_MAP = {
    'pr': {'cru': 'pre', 'era5': 'tp', 'units': 'mm yr⁻¹'},
    'tas': {'cru': 'tmp', 'era5': 't2m', 'units': '°C yr⁻¹'},
    'tasmax': {'cru': 'tmx', 'era5': 'tasmax', 'units': '°C yr⁻¹'},
    'tasmin': {'cru': 'tmn', 'era5': 'tasmin', 'units': '°C yr⁻¹'}
}

# Plot configuration
PLOT_CONFIG = {
    'pr': {
        'title': 'Precipitation trend (mm yr⁻¹)',
        'levels': np.arange(-6, 6.2, 0.2),
        'cmap': cm.BrBG,
        'sig_label': 'Dots: p < 0.05'
    },
    'tas': {
        'title': 'Air temperature trend (°C yr⁻¹)',
        'levels': np.arange(-0.3, 0.31, 0.01),
        'cmap': cm.bwr,
        'sig_label': 'Dots: p < 0.05'
    },
    'tasmax': {
        'title': 'Maximum air temperature trend (°C yr⁻¹)',
        'levels': np.arange(-0.3, 0.31, 0.01),
        'cmap': cm.bwr,
        'sig_label': 'Dots: p < 0.05'
    },
    'tasmin': {
        'title': 'Minimum air temperature trend (°C yr⁻¹)',
        'levels': np.arange(-0.3, 0.31, 0.01),
        'cmap': cm.bwr,
        'sig_label': 'Dots: p < 0.05'
    }
}


def import_regular_dataset(param, dataset):
    """Import CRU or ERA5 data (regular lat/lon grid)."""
    if dataset == 'CRU':
        var_name = VAR_MAP[param]['cru']
        lat_name, lon_name = 'lat', 'lon'
    elif dataset == 'ERA5':
        var_name = VAR_MAP[param]['era5']
        lat_name, lon_name = 'latitude', 'longitude'
    
    filename = f"{path}/{var_name}_{dataset}_year_{dt}.nc"
    print(f"Reading {dataset}: {filename}")
    
    with netCDF4.Dataset(filename) as nc:
        data = nc.variables[var_name][:]
        lat = nc.variables[lat_name][:]
        lon = nc.variables[lon_name][:]
        
        # Convert to 2D if 1D
        if len(lat.shape) == 1 and len(lon.shape) == 1:
            lon, lat = np.meshgrid(lon, lat)
        
        print(f"  Data shape: {data.shape}")
        print(f"  Lon range: [{lon.min():.1f}, {lon.max():.1f}]")
        print(f"  Lat range: [{lat.min():.1f}, {lat.max():.1f}]")
        
    return lat, lon, data


def import_regcm_dataset(param):
    """Import RegCM5 data with rotated pole coordinates."""
    filename = f"{path}/{param}_{dataset}_year_{dt}.nc"
    print(f"Reading RegCM5: {filename}")
    
    with netCDF4.Dataset(filename) as nc:
        data = nc.variables[param][:]
        
        # Try different coordinate variable names
        if 'lat' in nc.variables and 'lon' in nc.variables:
            lat = nc.variables['lat'][:]
            lon = nc.variables['lon'][:]
            print("  Using 2D lat/lon coordinates")
        elif 'rlat' in nc.variables and 'rlon' in nc.variables:
            # Rotated coordinates
            rlat = nc.variables['rlat'][:]
            rlon = nc.variables['rlon'][:]
            lon, lat = np.meshgrid(rlon, rlat)
            print("  Using rotated coordinates (rlat/rlon)")
        else:
            # Last resort: try common names
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
        
        # Check if we need to convert rotated coordinates
        if len(lat.shape) == 2 and len(lon.shape) == 2:
            # Already 2D coordinates
            print("  2D coordinates detected")
        elif len(lat.shape) == 1 and len(lon.shape) == 1:
            # Convert to 2D
            lon, lat = np.meshgrid(lon, lat)
        
        print(f"  Data shape: {data.shape}")
        print(f"  Lon shape: {lon.shape}, range: [{lon.min():.1f}, {lon.max():.1f}]")
        print(f"  Lat shape: {lat.shape}, range: [{lat.min():.1f}, {lat.max():.1f}]")
        
        # Handle longitude wrapping if needed
        if lon.min() < 0:
            lon = lon % 360
        
        return lat, lon, data


def convert_rotated_to_regular(lon_rot, lat_rot):
    """
    Convert rotated pole coordinates to regular lat/lon.
    Based on RegCM5's rotated pole at (60.31°S, 321.38°E) for AUS domain.
    """
    # Rotated pole parameters for AUS-12 domain
    rot_pole_lat = -60.31    # South pole latitude
    rot_pole_lon = 321.38    # South pole longitude
    
    # Convert to north pole convention for cartopy
    north_pole_lat = -rot_pole_lat      # 60.31°N
    north_pole_lon = rot_pole_lon - 180.0  # 141.38°E
    
    print(f"  Converting rotated pole: North pole at ({north_pole_lat:.2f}°N, {north_pole_lon:.2f}°E)")
    
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


def calculate_trend(data, alpha=0.05, per_decade=False):
    """Calculate linear trend and significance."""
    nt, ny, nx = data.shape
    t = np.arange(nt, dtype=float)
    
    trend = np.full((ny, nx), np.nan)
    sig_mask = np.full((ny, nx), False)
    
    for i in range(ny):
        for j in range(nx):
            y = data[:, i, j]
            valid = np.isfinite(y)
            
            if np.sum(valid) >= 10:  # Minimum points for trend
                t_valid = t[valid]
                y_valid = y[valid]
                
                # Linear regression
                A = np.vstack([t_valid, np.ones_like(t_valid)]).T
                slope, intercept = np.linalg.lstsq(A, y_valid, rcond=None)[0]
                
                # Calculate trend (per year by default)
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
    """Configure map plot settings."""
    ax.set_extent(DOMAIN_EXTENT, crs=ccrs.PlateCarree())
    
    # Set gridlines
    if domain == 'AUS-12':
        x_ticks = np.arange(110, 157, 10)
        y_ticks = np.arange(-45, -8, 5)
    else:  # CSAM-3
        x_ticks = np.arange(-80, -34, 10)
        y_ticks = np.arange(-35, -9, 5)
    
    ax.set_xticks(x_ticks, crs=ccrs.PlateCarree())
    ax.set_yticks(y_ticks, crs=ccrs.PlateCarree())
    
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    
    # Add features
    ax.add_feature(cfeat.COASTLINE, linewidth=0.6)
    ax.add_feature(cfeat.BORDERS, linewidth=0.5, linestyle='-', alpha=0.5)
    ax.add_feature(cfeat.LAND, facecolor='lightgray', alpha=0.2)
    
    # Grid
    ax.gridlines(draw_labels=False, linewidth=0.3, 
                 linestyle='--', color='gray', alpha=0.5)


def main():
    """Main processing function."""
    print(f"Processing variable: {var}")
    print(f"Time period: {dt}")
    print(f"Domain: {domain}")
    print("-" * 50)
    
    # Import data
    print("\n1. Importing data...")
    lat_cru, lon_cru, cru_data = import_regular_dataset(var, 'CRU')
    lat_era5, lon_era5, era5_data = import_regular_dataset(var, 'ERA5')
    lat_regcm, lon_regcm, regcm_data = import_regcm_dataset(var)
    
    # Process RegCM coordinates based on domain
    print("\n2. Processing RegCM coordinates...")
    
    if domain == 'AUS-12':
        # Convert rotated coordinates for AUS domain
        print("  Converting rotated coordinates to regular lat/lon for AUS domain...")
        lon_regcm, lat_regcm = convert_rotated_to_regular(lon_regcm, lat_regcm)
        
        # Adjust longitude to match domain extent [110, 156]
        # Convert from 0-360 to domain-specific range
        lon_regcm = np.where(lon_regcm > 180, lon_regcm - 360, lon_regcm)
    else:
        # CSAM-3: Coordinates are already regular lat/lon
        print("  CSAM domain uses regular lat/lon coordinates (no conversion needed)")
        
        # Adjust longitude to [-180, 180] if needed
        if lon_regcm.min() > 180:
            lon_regcm = np.where(lon_regcm > 180, lon_regcm - 360, lon_regcm)
    
    print(f"  Final coordinates:")
    print(f"    Lon range: [{lon_regcm.min():.2f}, {lon_regcm.max():.2f}]")
    print(f"    Lat range: [{lat_regcm.min():.2f}, {lat_regcm.max():.2f}]")
    
    # Mask data outside domain
    mask = ((lon_regcm >= DOMAIN_EXTENT[0]) & (lon_regcm <= DOMAIN_EXTENT[1]) &
            (lat_regcm >= DOMAIN_EXTENT[2]) & (lat_regcm <= DOMAIN_EXTENT[3]))
    regcm_data = np.where(mask[None, :, :], regcm_data, np.nan)
    
    # Calculate trends
    print("\n3. Calculating trends...")
    cru_trend, cru_sig = calculate_trend(cru_data, alpha=0.05)
    era5_trend, era5_sig = calculate_trend(era5_data, alpha=0.05)
    regcm_trend, regcm_sig = calculate_trend(regcm_data, alpha=0.05)
    
    print(f"  CRU trend range: [{np.nanmin(cru_trend):.3f}, {np.nanmax(cru_trend):.3f}]")
    print(f"  ERA5 trend range: [{np.nanmin(era5_trend):.3f}, {np.nanmax(era5_trend):.3f}]")
    print(f"  RegCM trend range: [{np.nanmin(regcm_trend):.3f}, {np.nanmax(regcm_trend):.3f}]")
    
    # Create figure
    print("\n4. Creating plot...")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5),
                            subplot_kw={'projection': ccrs.PlateCarree()})
    
    datasets = [
        ('CRU', lon_cru, lat_cru, cru_trend, cru_sig),
        ('ERA5', lon_era5, lat_era5, era5_trend, era5_sig),
        ('RegCM5', lon_regcm, lat_regcm, regcm_trend, regcm_sig)
    ]
    
    for idx, (ax, (name, lon, lat, trend, sig)) in enumerate(zip(axes, datasets)):
        # Plot trend - always use pcolormesh for RegCM5 due to curvilinear grid
        if name == 'RegCM5':
            im = ax.pcolormesh(lon, lat, trend,
                              cmap=PLOT_CONFIG[var]['cmap'],
                              vmin=PLOT_CONFIG[var]['levels'][0],
                              vmax=PLOT_CONFIG[var]['levels'][-1],
                              transform=ccrs.PlateCarree())
        else:
            # Use contourf for regular grids (CRU and ERA5)
            im = ax.contourf(lon, lat, trend,
                            levels=PLOT_CONFIG[var]['levels'],
                            cmap=PLOT_CONFIG[var]['cmap'],
                            extend='both',
                            transform=ccrs.PlateCarree())
        
        # Plot significance using contour hatching
        if np.any(sig):
            # Create a binary mask for significant points
            sig_binary = np.where(sig, 1.0, 0.0)
            
            # Use contourf with hatching to mark significant areas
            ax.contourf(lon, lat, sig_binary, 
                       levels=[0.5, 1.5], 
                       hatches=['...'], 
                       colors='none', 
                       transform=ccrs.PlateCarree(),
                       alpha=0.0)  # Set alpha to 0 to make fill transparent
        
        # Configure plot
        configure_plot(ax)
        
        # Title
        label = ['(a)', '(b)', '(c)'][idx]
        ax.set_title(f'{label} {name}', loc='left', 
                     fontsize=font_size+2, fontweight='bold')
    
    # Colorbar
    cbar = fig.colorbar(im, ax=axes, orientation='horizontal',
                       pad=0.07, aspect=40, shrink=0.8)
    cbar.set_label(f"{PLOT_CONFIG[var]['title']} ({PLOT_CONFIG[var]['sig_label']})",
                   fontsize=font_size+2, fontweight='bold')
    cbar.ax.tick_params(labelsize=font_size)
    
    # Save figure
    path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/figs/trend'
    os.makedirs(path_out, exist_ok=True)
    output_file = f'{path_out}/pyplt_maps_trend_{var}_{domain}_RegCM5_{dt}.png'
    print(f"\n5. Saving figure: {output_file}")
    plt.savefig(output_file, dpi=400, bbox_inches='tight')
    print("Done!")
    plt.show()

if __name__ == '__main__':
    main()
