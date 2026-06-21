# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 17, 2026"
__description__ = "This script plot Otis hurricane VTS"

import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.units import units
import warnings
warnings.filterwarnings('ignore')

# Domain type
DOMAIN = "large"
EXPS_ = 'exps_v3'
MAP_EXTENT = [-103, -92.5, 5, 17]  

# Experiment names
EXPERIMENTS = ["ctrl", "holt_r2", "holt_r3", "uw_r2", "uw_r3"]
HURRICANE_NAME = "Otis"

START_DATE = "2023-10-21"
END_DATE = "2023-10-25"

# Paths
DATA_DIR_MODELS = f"/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/exps/{EXPS_}/domain_{DOMAIN}_regridded/"
DATA_DIR_ERA5 = "/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/era5"
OUTPUT_PATH = f"/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/figs/{EXPS_}/"
os.makedirs(OUTPUT_PATH, exist_ok=True)


# Data loading function
def load_era5_data(data_dir, start_date=None, end_date=None):
    """Load ERA5 data from NetCDF files"""
    nc_files = sorted(glob.glob(os.path.join(data_dir, '*.nc')))
    
    if not nc_files:
        raise ValueError(f"No NetCDF files found in {data_dir}")
    
    print(f"Found {len(nc_files)} ERA5 files")
    
    # Load and merge all files
    datasets = [xr.open_dataset(file) for file in nc_files]
    ds_merged = xr.merge(datasets, compat='override')
    
    # Rename coordinates to standard names
    rename_dict = {
        "valid_time": "time",
        "latitude": "lat",
        "longitude": "lon",
        "pressure_level": "plev",
        "msl": "psl",
        "t": "ta",
        "q": "hus",
        "r": "rh"
    }
    
    for old_name, new_name in rename_dict.items():
        if old_name in ds_merged.variables:
            ds_merged = ds_merged.rename({old_name: new_name})
    
    # Fix longitude to [-180, 180] range
    if np.any(ds_merged.lon > 180):
        ds_merged.coords['lon'] = (ds_merged.lon + 180) % 360 - 180
        ds_merged = ds_merged.sortby('lon')
    
    # Select time range
    if start_date and end_date:
        ds_merged = ds_merged.sel(time=slice(start_date, end_date))
        print(f"Selected time range: {start_date} to {end_date}")
    
    return ds_merged


def load_regcm5_data(data_dir, experiment):
    """Load RegCM5 model data for a specific experiment"""
    data_path = os.path.join(data_dir, experiment)
    nc_files = sorted(glob.glob(os.path.join(data_path, '*.nc')))
    
    if not nc_files:
        raise ValueError(f"No NetCDF files found in {data_path}")
    
    # Load all datasets
    datasets = [xr.open_dataset(file) for file in nc_files]
    
    # Find intersection of time coordinates
    time_coords = [set(ds.time.values) for ds in datasets]
    common_times = sorted(time_coords[0].intersection(*time_coords[1:]))
    
    # Select only common time steps
    datasets_aligned = [ds.sel(time=common_times) for ds in datasets]
    
    # Merge datasets
    ds_merged = xr.merge(datasets_aligned, compat='override')
    
    # Fix longitude
    if np.any(ds_merged.lon > 180):
        ds_merged.coords['lon'] = (ds_merged.lon + 180) % 360 - 180
        ds_merged = ds_merged.sortby('lon')
    
    # Select time range
    ds_merged = ds_merged.sel(time=slice(START_DATE, END_DATE))
    
    return ds_merged


# Cyclone tracking function
def find_cyclone_center(ds_time, map_extent, is_era5=False):
    """Find cyclone center based on minimum sea level pressure"""
    lon_min, lon_max, lat_min, lat_max = map_extent
    
    # Handle latitude order
    if is_era5:
        # ERA5 latitude is descending (90 to -90)
        lat_slice = slice(lat_max, lat_min)
    else:
        # Model latitude is ascending
        lat_slice = slice(lat_min, lat_max)
    
    try:
        ds_domain = ds_time.sel(lat=lat_slice, lon=slice(lon_min, lon_max))
    except:
        return None
    
    if ds_domain.lat.size == 0 or ds_domain.lon.size == 0:
        return None
    
    # Get sea level pressure (convert from Pa to hPa if necessary)
    psl = ds_domain['psl']
    if hasattr(psl, 'units') and psl.units.lower() == 'pa':
        psl = psl / 100.0
    
    # Find minimum pressure location
    min_pressure = float(psl.min())
    
    if min_pressure > 1010:  # Threshold for cyclone
        return None
    
    min_point = psl.where(psl == psl.min(), drop=True)
    
    if min_point.lat.size == 0 or min_point.lon.size == 0:
        return None
    
    lat_center = float(min_point.lat.values[0])
    lon_center = float(min_point.lon.values[0])
    
    return {
        'lat': lat_center,
        'lon': lon_center,
        'min_pressure': min_pressure
    }


def track_cyclone(ds, map_extent, exp_name):
    """Track cyclone to find peak intensity"""
    track_data = []
    is_era5 = (exp_name == "ERA5")
    
    for t in ds.time.values:
        ds_time = ds.sel(time=t)
        center = find_cyclone_center(ds_time, map_extent, is_era5)
        
        if center is not None:
            track_data.append({
                'time': pd.to_datetime(t),
                'lat': center['lat'],
                'lon': center['lon'],
                'min_pressure': center['min_pressure']
            })
    
    if not track_data:
        print(f"Warning: No cyclone centers found for {exp_name}")
        return None
    
    df_track = pd.DataFrame(track_data)
    
    # Find peak intensity (minimum pressure)
    min_idx = df_track['min_pressure'].idxmin()
    peak = df_track.loc[min_idx]
    
    print(f"{exp_name:10} - Time: {peak['time']}, Lat: {peak['lat']:.2f}°N, Lon: {peak['lon']:.2f}°W, Pressure: {peak['min_pressure']:.1f} hPa")
    
    return peak


# Vertical thermal structure function
def compute_equivalent_potential_temperature(ds_slice, is_era5=False):
    """Compute equivalent potential temperature (theta_e)"""
    T = ds_slice['ta']  # Temperature in K
    p = ds_slice['plev']  # Pressure levels in hPa
    
    if is_era5:
        # ERA5 uses relative humidity
        rh = ds_slice['rh'].clip(min=0.0, max=100)
        T_units = T * units.kelvin
        p_units = p.broadcast_like(T) * units.hPa
        rh_units = rh * units.percent
        
        dewpoint = mpcalc.dewpoint_from_relative_humidity(T_units, rh_units)
    else:
        # Models use specific humidity
        q = ds_slice['hus']  # Specific humidity in kg/kg
        p_broadcasted = p.broadcast_like(T)
        
        T_units = T * units.kelvin
        p_units = p_broadcasted * units.hPa
        q_units = q * units('kg/kg')
        
        dewpoint = mpcalc.dewpoint_from_specific_humidity(p_units, T_units, q_units)
    
    # Compute equivalent potential temperature
    theta_e = mpcalc.equivalent_potential_temperature(p_units, T_units, dewpoint)
    
    return theta_e

def extract_cross_section(ds, center_time, center_lat, center_lon, is_era5=False, radius_deg=5):
    """Extract a zonal cross-section through the cyclone center"""
    # Select time
    ds_time = ds.sel(time=center_time, method='nearest')
    
    # Select longitude range around center
    lon_min = center_lon - radius_deg
    lon_max = center_lon + radius_deg
    
    # Select the closest latitude to center
    lat_closest = ds_time.lat.sel(lat=center_lat, method='nearest')
    
    # Extract the cross-section
    ds_cross = ds_time.sel(
        lon=slice(lon_min, lon_max),
        lat=lat_closest
    )
    
    # Compute signed distance from center (km)
    R = 6371.0
    center_lon_rad = np.radians(center_lon)
    center_lat_rad = np.radians(center_lat)
    lon_rad = np.radians(ds_cross.lon.values)
    
    # Haversine formula for east-west distance
    distances = R * np.cos(center_lat_rad) * (lon_rad - center_lon_rad)
    
    # Add distance as coordinate
    ds_cross = ds_cross.assign_coords(r=("lon", distances))
    
    # Remove any singleton dimensions
    ds_cross = ds_cross.squeeze()
    
    return ds_cross

# Plot function
def plot_all_vts_subplots(all_vts_data, output_file=None):
    """Create a single figure with subplots for all experiments"""
    
    # Define subplot layout (2 rows, 3 columns)
    fig, axes = plt.subplots(2, 3, figsize=(14, 10))
    axes = axes.flatten()
    
    # Common contour levels
    levels = np.arange(314, 390, 2)
    
    # List of labels for each subplot
    labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
    
    # Plot each experiment
    for idx, (exp_name, data) in enumerate(all_vts_data.items()):
        ax = axes[idx]
        
        # Get data
        Z = np.squeeze(data['theta_e'].values)
        r_coords = data['theta_e'].r.values
        plev_coords = data['theta_e'].plev.values
        
        # Create meshgrid
        X, Y = np.meshgrid(r_coords, plev_coords)
        
        # Create filled contour plot
        cf = ax.contourf(X, Y, Z, levels=levels, cmap='jet', extend='max')
        
        # Add contour lines
        cs = ax.contour(X, Y, Z, levels=levels[::3], colors='black', linewidths=0.5, alpha=0.5)
        ax.clabel(cs, inline=True, fontsize=7, fmt='%d')
        
        # Formatting
        ax.set_xlabel('Distance from Center (km)', fontsize=10, fontweight='bold')
        ax.set_ylabel('Pressure (hPa)', fontsize=10, fontweight='bold')
        ax.set_ylim(1000, 100)
        ax.axvline(x=0, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.grid(True, alpha=0.5, linestyle='--')
        
        # Title with experiment info and label
        time_str = pd.to_datetime(data['time']).strftime('%Y-%m-%d %H:%M')
        ax.set_title(f'{labels[idx]} {exp_name.upper()} ({data["lat"]:.1f}°N, {abs(data["lon"]):.1f}°W) {time_str}', fontsize=10, fontweight='bold')
    
    # Hide unused subplot (if any)
    for idx in range(len(all_vts_data), len(axes)):
        axes[idx].set_visible(False)
    
    # Add colorbar
    cbar_ax = fig.add_axes([0.999, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(cf, cbar_ax)
    cbar.set_label('Equivalent Potential Temperature (K)', fontsize=10, fontweight='bold')

    plt.tight_layout()    
    if output_file:
        plt.savefig(output_file, dpi=400, bbox_inches='tight')
        print(f"\nVTS comparison plot saved to {output_file}")
    
    plt.show()
    return fig


# Main code
def main():
    """Main function to compute and plot VTS for all experiments in one figure"""
    
    print("=" * 80)
    print(f"Computing Vertical Thermal Structure for Hurricane {HURRICANE_NAME}")
    print(f"Experiments: {EXPERIMENTS + ['ERA5']}")
    print("=" * 80)
    
    # Dictionary to store all results
    all_vts_data = {}
    
    # 1. Process ERA5
    print("\nLoading ERA5 data...")
    ds_era5 = load_era5_data(DATA_DIR_ERA5, START_DATE, END_DATE)
    
    print("Tracking ERA5 cyclone...")
    peak_era5 = track_cyclone(ds_era5, MAP_EXTENT, "ERA5")
    
    if peak_era5 is not None:
        print("Extracting ERA5 cross-section...")
        ds_cross_era5 = extract_cross_section(ds_era5, peak_era5['time'], peak_era5['lat'], 
                                              peak_era5['lon'], is_era5=True, radius_deg=5)
        
        print("Computing ERA5 theta_e...")
        theta_e_era5 = compute_equivalent_potential_temperature(ds_cross_era5, is_era5=True)
        
        all_vts_data["ERA5"] = {
            'theta_e': theta_e_era5,
            'lat': peak_era5['lat'],
            'lon': peak_era5['lon'],
            'time': peak_era5['time']
        }
    
    # 2. Process each model experiment
    for exp in EXPERIMENTS:
        print(f"\nLoading {exp} data...")
        ds_model = load_regcm5_data(DATA_DIR_MODELS, exp)
        
        print(f"Tracking {exp} cyclone...")
        peak_model = track_cyclone(ds_model, MAP_EXTENT, exp)
        
        if peak_model is not None:
            print(f"Extracting {exp} cross-section...")
            ds_cross_model = extract_cross_section(ds_model, peak_model['time'], peak_model['lat'],
                                                   peak_model['lon'], is_era5=False, radius_deg=5)
            
            print(f"Computing {exp} theta_e...")
            theta_e_model = compute_equivalent_potential_temperature(ds_cross_model, is_era5=False)
            
            all_vts_data[exp] = {
                'theta_e': theta_e_model,
                'lat': peak_model['lat'],
                'lon': peak_model['lon'],
                'time': peak_model['time']
            }
    
    # 3. Create single figure with all subplots
    if len(all_vts_data) > 0:
        print("\n" + "="*40)
        print("Creating single figure with all VTS subplots")
        print("="*40)
        
        output_file = os.path.join(OUTPUT_PATH, f'pyplt_Hurricane_Otis_{HURRICANE_NAME.lower()}_vts_{DOMAIN}.png')
        plot_all_vts_subplots(all_vts_data, output_file)
        
        # Summary
        print("\n" + "="*80)
        print("VTS Analysis Complete!")
        print("="*80)
        print("\nPeak Intensity Summary:")
        for exp, data in all_vts_data.items():
            print(f"  {exp.upper():10} - {data['time']} - ({data['lat']:.2f}°N, {data['lon']:.2f}°W)")
        print(f"\nFigure saved to: {output_file}")
        print("="*80)
    else:
        print("Error: No data was successfully processed!")

if __name__ == "__main__":
    main()
