# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jan 25, 2026"
__description__ = "This script plot trend maps with rotated pole support"

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
parser.add_argument('--var', required=True, help='Variable name (pr, tas, tasmax, tasmin)')
parser.add_argument('--domain', required=True, help='Domain name')
parser.add_argument('--idt', required=True, help='Initial year')
parser.add_argument('--fdt', required=True, help='Final year')
args = parser.parse_args()

# Argparse
var = args.var
domain = args.domain
idt = args.idt
fdt = args.fdt
dt = f'{idt}-{fdt}'
dataset = f'{domain}_RegCM5'

font_size = 8
path = '/leonardo_scratch/large/userexternal/fraffael/plots/trend/inputs'
#path = '/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/postproc/evaluate/trend'

# Domain extent for AUS-12
DOMAIN_EXTENT = [110, 156, -45, -9]  # [lon_min, lon_max, lat_min, lat_max]

# Domain extent for  CSAM-3
#DOMAIN_EXTENT = [-79, -35, -35, -10]  # [lon_min, lon_max, lat_min, lat_max]

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
        'levels': np.arange(-60, 62, 2),
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


def extract_rotated_pole_params(nc_file):
    """Extract rotated pole parameters from NetCDF file attributes."""
    # Try to get from global attributes
    if hasattr(nc_file, 'proj_params'):
        proj_str = nc_file.proj_params
    elif 'mapping' in nc_file.variables:
        mapping_var = nc_file.variables['mapping']
        if hasattr(mapping_var, 'grid_mapping_name'):
            if mapping_var.grid_mapping_name == 'rotated_latitude_longitude':
                proj_str = f"+proj=ob_tran +o_proj=longlat +o_lat_p={mapping_var.grid_north_pole_latitude} +o_lon_p={mapping_var.grid_north_pole_longitude} +R=6371229."
    else:
        # Try to find in variable attributes
        for var_name in nc_file.variables:
            var = nc_file.variables[var_name]
            if hasattr(var, 'grid_mapping'):
                mapping_var = nc_file.variables[var.grid_mapping]
                if hasattr(mapping_var, 'grid_mapping_name'):
                    if mapping_var.grid_mapping_name == 'rotated_latitude_longitude':
                        proj_str = f"+proj=ob_tran +o_proj=longlat +o_lat_p={mapping_var.grid_north_pole_latitude} +o_lon_p={mapping_var.grid_north_pole_longitude} +R=6371229."
                        break
    
    # Parse the projection string
    params = {}
    if 'proj_str' in locals():
        parts = proj_str.split()
        for part in parts:
            if '=' in part:
                key, value = part.split('=')
                params[key] = value
        print(f"  Found projection parameters: {params}")
    else:
        print("  Warning: Could not find projection parameters, using defaults")
        params = {'o_lat_p': '-70.60', 'o_lon_p': '123.94'}  # Your new parameters
    
    return float(params.get('o_lat_p', -70.60)), float(params.get('o_lon_p', 123.94))


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
        
        # Extract rotated pole parameters
        rot_pole_lat, rot_pole_lon = extract_rotated_pole_params(nc)
        
        # Try to get rotated coordinates first
        if 'rlat' in nc.variables and 'rlon' in nc.variables:
            rlat = nc.variables['rlat'][:]  # rotated latitude
            rlon = nc.variables['rlon'][:]  # rotated longitude
            lon_rot, lat_rot = np.meshgrid(rlon, rlat)
            print(f"  Using rlon/rlat coordinates: {rlon.shape}, {rlat.shape}")
            print(f"  Rotated lon range: [{rlon.min():.1f}, {rlon.max():.1f}]")
            print(f"  Rotated lat range: [{rlat.min():.1f}, {rlat.max():.1f}]")
            
            # Convert to regular coordinates
            lon_reg, lat_reg = convert_rotated_to_regular(lon_rot, lat_rot, 
                                                          rot_pole_lat, rot_pole_lon)
            
            return lat_reg, lon_reg, data, rot_pole_lat, rot_pole_lon
            
        elif 'lat' in nc.variables and 'lon' in nc.variables:
            # Already regular coordinates
            lat = nc.variables['lat'][:]
            lon = nc.variables['lon'][:]
            print("  Using regular lat/lon coordinates")
            
            return lat, lon, data, rot_pole_lat, rot_pole_lon
        
        else:
            raise ValueError(f"Cannot find coordinates in {filename}")


def convert_rotated_to_regular(lon_rot, lat_rot, rot_pole_lat, rot_pole_lon):
    """
    Convert rotated pole coordinates to regular lat/lon.
    """
    print(f"  Converting rotated pole: South pole at ({rot_pole_lat:.2f}°N, {rot_pole_lon:.2f}°E)")
    
    # Convert to north pole convention for cartopy
    north_pole_lat = -rot_pole_lat      # Convert south pole to north pole
    north_pole_lon = rot_pole_lon - 180.0  # Convert to -180 to 180 range
    
    print(f"  North pole at: ({north_pole_lat:.2f}°N, {north_pole_lon:.2f}°E)")
    
    # Create coordinate systems
    rotated_crs = ccrs.RotatedPole(pole_latitude=north_pole_lat,
                                   pole_longitude=north_pole_lon)
    regular_crs = ccrs.PlateCarree()
    
    # Transform coordinates
    print(f"  Transforming {lon_rot.size} points...")
    points = rotated_crs.transform_points(regular_crs, lon_rot, lat_rot)
    lon_reg = points[..., 0]
    lat_reg = points[..., 1]
    
    # Ensure longitude is in -180 to 180 range for South America
    lon_reg = np.where(lon_reg > 180, lon_reg - 360, lon_reg)
    
    print(f"  Converted lon range: [{lon_reg.min():.1f}, {lon_reg.max():.1f}]")
    print(f"  Converted lat range: [{lat_reg.min():.1f}, {lat_reg.max():.1f}]")
    
    return lon_reg, lat_reg


def calculate_trend(data, alpha=0.05, per_decade=False):
    """Calculate linear trend and significance."""
    nt, ny, nx = data.shape
    t = np.arange(nt, dtype=float)
    
    trend = np.full((ny, nx), np.nan)
    sig_mask = np.full((ny, nx), False)
    
    # Vectorized calculation where possible
    for i in range(ny):
        for j in range(nx):
            y = data[:, i, j]
            valid = np.isfinite(y)
            n_valid = np.sum(valid)
            
            if n_valid >= 10:  # Minimum points for trend
                t_valid = t[valid]
                y_valid = y[valid]
                
                # Linear regression
                t_mean = np.mean(t_valid)
                y_mean = np.mean(y_valid)
                
                # Calculate slope
                numerator = np.sum((t_valid - t_mean) * (y_valid - y_mean))
                denominator = np.sum((t_valid - t_mean) ** 2)
                
                if denominator > 0:
                    slope = numerator / denominator
                    
                    # Calculate trend (per year by default)
                    trend_val = slope * 10.0 if per_decade else slope
                    trend[i, j] = trend_val
                    
                    # Calculate p-value
                    y_pred = slope * t_valid + y_mean - slope * t_mean
                    residuals = y_valid - y_pred
                    ss_res = np.sum(residuals ** 2)
                    
                    if n_valid > 2 and ss_res > 0:
                        # Standard error of the slope
                        se_slope = np.sqrt(ss_res / ((n_valid - 2) * denominator))
                        
                        # t-statistic
                        t_stat = slope / se_slope
                        
                        # p-value (two-tailed)
                        p_value = 2 * (1 - student_t.cdf(np.abs(t_stat), df=n_valid - 2))
                        
                        if p_value < alpha:
                            sig_mask[i, j] = True
    
    return trend, sig_mask


def configure_plot(ax):
    """Configure map plot settings for South America."""
    ax.set_extent(DOMAIN_EXTENT, crs=ccrs.PlateCarree())
    
    # Set gridlines - adjust spacing for South America
    lon_ticks = np.arange(-80, -30, 10)
    lat_ticks = np.arange(-40, -5, 5)
    
    ax.set_xticks(lon_ticks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_ticks, crs=ccrs.PlateCarree())
    
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    
    # Add features
    ax.add_feature(cfeat.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeat.BORDERS, linewidth=0.6, linestyle='-', alpha=0.7)
    ax.add_feature(cfeat.LAND, facecolor='lightgray', alpha=0.1)
    ax.add_feature(cfeat.OCEAN, facecolor='lightblue', alpha=0.1)
    
    # Grid
    ax.gridlines(draw_labels=False, linewidth=0.4, 
                 linestyle='--', color='gray', alpha=0.5)


def main():
    """Main processing function."""
    print(f"Processing variable: {var}")
    print(f"Time period: {dt}")
    print(f"Domain: {domain}")
    print(f"Domain extent: {DOMAIN_EXTENT}")
    print("-" * 50)
    
    # Import data
    print("\n1. Importing data...")
    lat_cru, lon_cru, cru_data = import_regular_dataset(var, 'CRU')
    lat_era5, lon_era5, era5_data = import_regular_dataset(var, 'ERA5')
    
    # For RegCM, you might need to adjust the dataset name
    # If it's not 'AUS-12_RegCM5', update the import_regcm_dataset function
    lat_regcm, lon_regcm, regcm_data, rot_lat, rot_lon = import_regcm_dataset(var)
    
    print(f"\n2. Rotated pole parameters used:")
    print(f"   Rotated pole latitude: {rot_lat:.2f}")
    print(f"   Rotated pole longitude: {rot_lon:.2f}")
    
    # Mask data outside domain
    print("\n3. Masking data outside domain...")
    mask = ((lon_regcm >= DOMAIN_EXTENT[0]) & (lon_regcm <= DOMAIN_EXTENT[1]) &
            (lat_regcm >= DOMAIN_EXTENT[2]) & (lat_regcm <= DOMAIN_EXTENT[3]))
    regcm_data = np.where(mask[None, :, :], regcm_data, np.nan)
    
    # Calculate trends
    print("\n4. Calculating trends...")
    cru_trend, cru_sig = calculate_trend(cru_data, alpha=0.05)
    era5_trend, era5_sig = calculate_trend(era5_data, alpha=0.05)
    regcm_trend, regcm_sig = calculate_trend(regcm_data, alpha=0.05)
    
    print(f"  CRU trend range: [{np.nanmin(cru_trend):.3f}, {np.nanmax(cru_trend):.3f}]")
    print(f"  ERA5 trend range: [{np.nanmin(era5_trend):.3f}, {np.nanmax(era5_trend):.3f}]")
    print(f"  RegCM trend range: [{np.nanmin(regcm_trend):.3f}, {np.nanmax(regcm_trend):.3f}]")
    
    # Create figure
    print("\n5. Creating plot...")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5),
                            subplot_kw={'projection': ccrs.PlateCarree()})
    
    datasets = [
        ('CRU', lon_cru, lat_cru, cru_trend, cru_sig),
        ('ERA5', lon_era5, lat_era5, era5_trend, era5_sig),
        ('RegCM5', lon_regcm, lat_regcm, regcm_trend, regcm_sig)
    ]
    
    for idx, (ax, (name, lon, lat, trend, sig)) in enumerate(zip(axes, datasets)):
        # Plot trend
        if name == 'RegCM5':
            # Use pcolormesh for irregular grids
            im = ax.pcolormesh(lon, lat, trend,
                              cmap=PLOT_CONFIG[var]['cmap'],
                              vmin=PLOT_CONFIG[var]['levels'][0],
                              vmax=PLOT_CONFIG[var]['levels'][-1],
                              transform=ccrs.PlateCarree())
        else:
            # Use contourf for regular grids
            im = ax.contourf(lon, lat, trend,
                            levels=PLOT_CONFIG[var]['levels'],
                            cmap=PLOT_CONFIG[var]['cmap'],
                            extend='both',
                            transform=ccrs.PlateCarree())
        
        # Plot significance dots
        if np.any(sig):
            # For efficiency, sample points if there are too many
            sig_positions = np.where(sig)
            if len(sig_positions[0]) > 10000:
                # Random sample for display
                sample_idx = np.random.choice(len(sig_positions[0]), 10000, replace=False)
                sig_lon = lon[sig_positions[0][sample_idx], sig_positions[1][sample_idx]]
                sig_lat = lat[sig_positions[0][sample_idx], sig_positions[1][sample_idx]]
            else:
                sig_lon = lon[sig]
                sig_lat = lat[sig]
            
            ax.scatter(sig_lon, sig_lat,
                      s=1, color='black', alpha=0.3,
                      transform=ccrs.PlateCarree(),
                      marker='.', linewidths=0)
        
        # Configure plot
        configure_plot(ax)
        
        # Title
        label = ['(a)', '(b)', '(c)'][idx]
        ax.set_title(f'{label} {name}', loc='left', 
                     fontsize=font_size+2, fontweight='bold')
    
    # Colorbar
    cbar = fig.colorbar(im, ax=axes, orientation='horizontal',
                       pad=0.06, aspect=40, shrink=0.8)
    cbar.set_label(f"{PLOT_CONFIG[var]['title']} ({PLOT_CONFIG[var]['sig_label']})",
                   fontsize=font_size+2, fontweight='bold')
    cbar.ax.tick_params(labelsize=font_size)
        
    # Save figure
    output_file = f'pyplt_maps_trend_{var}_{domain}_RegCM5_{dt}.png'
    print(f"\n6. Saving figure: {output_file}")
    plt.savefig(output_file, dpi=400, bbox_inches='tight', facecolor='white')
    print("Done!")
    
    plt.show()

if __name__ == '__main__':
    main()
