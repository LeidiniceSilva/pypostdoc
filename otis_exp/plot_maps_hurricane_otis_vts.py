#!/usr/bin/env python3
# %%
__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 17, 2026"
__description__ = "Plot Otis hurricane VTS (Vertical Thermal Structure) - Combined figure with ERA5"

import os
import glob
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from metpy.units import units
import metpy.calc as mpcalc

domain = 'small'
OUTPUT_PATH = "/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/figs/"

# %%
# FUNCTIONS NECESSARY FOR VTS

def load_regcm5_multi_file(data_dir):
    """Load RegCM5 data from separate variable files and merge."""
    nc_files = sorted(glob.glob(os.path.join(data_dir, '*.nc')))
    if not nc_files:
        raise ValueError(f"No NetCDF files found in {data_dir}")
    
    datasets = [xr.open_dataset(file) for file in nc_files]
    time_coords = [set(ds.time.values) for ds in datasets]
    common_times = sorted(time_coords[0].intersection(*time_coords[1:]))
    datasets_aligned = [ds.sel(time=common_times) for ds in datasets]
    ds_merged = xr.merge(datasets_aligned, compat='override')
    
    return ds_merged

def load_era5_vts(era5_dir, map_extent):
    """Load ERA5 temperature and specific humidity and crop to hurricane domain."""
    # Load temperature
    t_file = os.path.join(era5_dir, 't_ERA5_Otis_6hr_Oct2023.nc')
    q_file = os.path.join(era5_dir, 'q_ERA5_Otis_6hr_Oct2023.nc')
    
    ds_t = xr.open_dataset(t_file)
    ds_q = xr.open_dataset(q_file)
    
    # Rename coordinates
    ds_t = ds_t.rename({
        'longitude': 'lon',
        'latitude': 'lat',
        'valid_time': 'time'
    })
    ds_q = ds_q.rename({
        'longitude': 'lon',
        'latitude': 'lat',
        'valid_time': 'time'
    })
    
    # Fix longitude to [-180, 180]
    if np.any(ds_t.lon > 180):
        ds_t.coords['lon'] = (ds_t.coords['lon'] + 180) % 360 - 180
        ds_t = ds_t.sortby(ds_t.lon)
    if np.any(ds_q.lon > 180):
        ds_q.coords['lon'] = (ds_q.coords['lon'] + 180) % 360 - 180
        ds_q = ds_q.sortby(ds_q.lon)
    
    # Crop to hurricane domain (add some margin)
    lon_min, lon_max, lat_min, lat_max = map_extent
    # Add 5 degrees margin for safety
    lon_min_crop = lon_min - 5
    lon_max_crop = lon_max + 5
    lat_min_crop = lat_min - 5
    lat_max_crop = lat_max + 5
    
    ds_t = ds_t.sel(
        lat=slice(lat_min_crop, lat_max_crop),
        lon=slice(lon_min_crop, lon_max_crop)
    )
    ds_q = ds_q.sel(
        lat=slice(lat_min_crop, lat_max_crop),
        lon=slice(lon_min_crop, lon_max_crop)
    )
    
    # Merge temperature and specific humidity
    ds_era = xr.merge([ds_t, ds_q])
    
    print(f"ERA5 cropped to domain: lon=[{lon_min_crop:.1f}, {lon_max_crop:.1f}], lat=[{lat_min_crop:.1f}, {lat_max_crop:.1f}]")
    print(f"ERA5 shape: time={len(ds_era.time)}, pressure_level={len(ds_era.pressure_level)}, lat={len(ds_era.lat)}, lon={len(ds_era.lon)}")
    
    return ds_era

def haversine_distance(lat1, lon1, lat2, lon2):
    """Compute great-circle distance (km) between two lat/lon points."""
    R = 6371.0
    phi1, phi2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlambda = np.radians(lon2 - lon1)
    a = np.sin(dphi/2)**2 + np.cos(phi1)*np.cos(phi2)*np.sin(dlambda/2)**2
    return 2 * R * np.arcsin(np.sqrt(a))

def find_cyclone_centers(ds_models, map_extent, psl_threshold=1010.0, max_jump_km=300):
    """Find cyclone centers for each experiment (including ERA5)."""
    cyclone_centers = {}
    minimum_datetime = {}
    
    for exp_name, ds in ds_models.items():
        print(f"Processing: {exp_name}")
        
        # Fix longitude if needed
        if np.any(ds.lon > 180):
            ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
            ds = ds.sortby(ds.lon)
        
        centers = []
        last_lat, last_lon = None, None
        
        for time_val in ds.time.values:
            try:
                ds_sel = ds.sel(
                    lat=slice(map_extent[2], map_extent[3]),
                    lon=slice(map_extent[0], map_extent[1]),
                    time=time_val
                )
            except Exception as e:
                print(f"  Warning: Could not select data for time {time_val}: {e}")
                continue
            
            # Get pressure variable
            if 'psl' in ds_sel:
                pressure = ds_sel['psl']
                if hasattr(pressure, 'attrs') and pressure.attrs.get('units', '').lower() == 'pa':
                    pressure = pressure / 100.0
                use_pressure = True
            elif 'msl' in ds_sel:
                pressure = ds_sel['msl']
                if hasattr(pressure, 'attrs') and pressure.attrs.get('units', '').lower() == 'pa':
                    pressure = pressure / 100.0
                use_pressure = True
            else:
                # For ERA5 without msl/psl, use temperature at 850hPa as proxy
                if 't' in ds_sel:
                    # Select temperature at 850hPa
                    if 'pressure_level' in ds_sel.dims:
                        t_850 = ds_sel['t'].sel(pressure_level=850, method='nearest')
                        pressure = t_850
                        use_pressure = False
                    else:
                        print(f"  Warning: No pressure_level dimension for {exp_name}")
                        continue
                else:
                    print(f"  Warning: No pressure variable found for {exp_name}")
                    continue
            
            if np.all(np.isnan(pressure.values)):
                continue
            
            # Find minimum pressure (or minimum temperature for ERA5 proxy)
            min_value = pressure.min().item()
            min_point = pressure.where(pressure == min_value, drop=True)
            
            if min_point.lat.size == 0 or min_point.lon.size == 0:
                continue
                
            lat_center = float(min_point.lat.values.flatten()[0])
            lon_center = float(min_point.lon.values.flatten()[0])
            
            valid = True
            if use_pressure:
                if min_value > psl_threshold:
                    valid = False
            if valid and last_lat is not None and last_lon is not None:
                dist = haversine_distance(last_lat, last_lon, lat_center, lon_center)
                if dist > max_jump_km:
                    valid = False
            
            if not valid:
                lat_center, lon_center, min_value = np.nan, np.nan, np.nan
            
            centers.append({
                'time': pd.to_datetime(time_val),
                'lat': lat_center,
                'lon': lon_center,
                'min_pressure': min_value if use_pressure else np.nan
            })
            
            if valid and not np.isnan(lat_center):
                last_lat, last_lon = lat_center, lon_center
        
        if len(centers) == 0:
            print(f"  Warning: No cyclone centers found for {exp_name}")
            cyclone_centers[exp_name] = pd.DataFrame()
            minimum_datetime[exp_name] = {'time': None, 'lat': np.nan, 'lon': np.nan, 'min_pressure': np.nan}
            continue
        
        df = pd.DataFrame(centers)
        cyclone_centers[exp_name] = df
        
        # Find index of minimum pressure (or minimum temperature for ERA5)
        if 'min_pressure' in df.columns and not df['min_pressure'].isna().all():
            # Filter out NaN values before finding idxmin
            valid_df = df[df['min_pressure'].notna()]
            if len(valid_df) > 0:
                min_idx = valid_df["min_pressure"].idxmin()
                minimum_datetime[exp_name] = df.iloc[min_idx]
            else:
                # Use the first valid center
                valid_df = df[df['lat'].notna()]
                if len(valid_df) > 0:
                    minimum_datetime[exp_name] = valid_df.iloc[len(valid_df)//2]
                else:
                    minimum_datetime[exp_name] = {'time': None, 'lat': np.nan, 'lon': np.nan, 'min_pressure': np.nan}
        else:
            # Use the middle point
            valid_df = df[df['lat'].notna()]
            if len(valid_df) > 0:
                min_idx = len(valid_df) // 2
                minimum_datetime[exp_name] = valid_df.iloc[min_idx]
            else:
                minimum_datetime[exp_name] = {'time': None, 'lat': np.nan, 'lon': np.nan, 'min_pressure': np.nan}
    
    return cyclone_centers, minimum_datetime

def compute_equiv_potential_temp_regcm5(da):
    """Compute equivalent potential temperature for RegCM5."""
    T = da['ta']  # K
    q = da['hus']  # RegCM5 specific humidity
    p = da['plev']
    p_broadcasted = p.broadcast_like(T)
    
    dpt = mpcalc.dewpoint_from_specific_humidity(p_broadcasted * units.hPa, q * units("kg/kg"))
    theta_e = mpcalc.equivalent_potential_temperature(p_broadcasted * units.hPa, T * units.K, dpt)
    return theta_e

def compute_equiv_potential_temp_era5(da):
    """Compute equivalent potential temperature for ERA5."""
    T = da['t']  # K
    q = da['q']  # ERA5 specific humidity (kg/kg)
    p = da['pressure_level']  # ERA5 pressure level
    p_broadcasted = p.broadcast_like(T)
    
    dpt = mpcalc.dewpoint_from_specific_humidity(p_broadcasted * units.hPa, q * units("kg/kg"))
    theta_e = mpcalc.equivalent_potential_temperature(p_broadcasted * units.hPa, T * units.K, dpt)
    return theta_e

def haversine_signed_1d_from_center(lat, lon):
    """Compute signed distance (km) from central point along a cross-section."""
    R = 6371.0
    idx_c = len(lat) // 2
    lat_c = lat[idx_c]
    lon_c = lon[idx_c]
    
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)
    lat_c_rad = np.radians(lat_c)
    lon_c_rad = np.radians(lon_c)
    
    delta_lat = lat_rad - lat_c_rad
    delta_lon = lon_rad - lon_c_rad
    
    dy = R * delta_lat
    dx = R * np.cos((lat_rad + lat_c_rad) / 2) * delta_lon
    signs = np.sign(dy + dx)
    radius_km = signs * np.sqrt(dx**2 + dy**2)
    
    return radius_km

def plot_vts_combined(ds_models, minimum_datetime, experiments, figname=None):
    """Plot VTS for all experiments in a single figure with subplots."""
    
    # Create list with all datasets (ERA5 first, then experiments)
    all_datasets = []
    
    # Add ERA5 if available
    if 'ERA5' in ds_models and 'ERA5' in minimum_datetime:
        center = minimum_datetime['ERA5']
        if center['time'] is not None and not np.isnan(center['lat']):
            all_datasets.append(('ERA5', ds_models['ERA5'], center))
    
    # Add experiments
    for exp in experiments:
        if exp in ds_models and exp in minimum_datetime:
            center = minimum_datetime[exp]
            if center['time'] is not None and not np.isnan(center['lat']):
                all_datasets.append((exp, ds_models[exp], center))
    
    n_exps = len(all_datasets)
    if n_exps == 0:
        print("No valid experiments found")
        return None, None
    
    # Subplot grid (2 rows, 3 columns for up to 6 subplots)
    n_cols = min(3, n_exps)
    n_rows = (n_exps + n_cols - 1) // n_cols
    
    plt.rcParams.update({
        "font.size": 12,
        "axes.labelsize": 14,
        "axes.titlesize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10
    })
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))
    if n_rows == 1 and n_cols == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # First pass: compute global min/max for colorbar
    vmin_all, vmax_all = float('inf'), float('-inf')
    theta_e_data = {}
    
    for idx, (exp_name, ds, center) in enumerate(all_datasets):
        mtime = center['time']
        mlat = center['lat']
        mlon = center['lon']
        
        if np.isnan(mlat) or np.isnan(mlon) or mtime is None:
            print(f"Warning: Invalid center for {exp_name}")
            continue
        
        print(f"Processing VTS for {exp_name}: time={mtime}, lat={mlat:.2f}, lon={mlon:.2f}")
        
        # Fix longitude if needed
        if np.any(ds.lon > 180):
            ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
            ds = ds.sortby(ds.lon)
        
        # Extract cross-section
        lons = ds.lon.values
        lon_mask = (lons > mlon - 3) & (lons < mlon + 3)
        
        try:
            cs = ds.sel(time=mtime, lat=mlat, lon=lon_mask, method="nearest")
            
            # Compute distance
            r = haversine_signed_1d_from_center(
                np.array([cs.lat.values.item()] * len(cs.lon.values)), 
                cs.lon.values
            )
            
            # Add distance coordinate
            if exp_name == 'ERA5':
                cs = cs.assign_coords(r=("lon", r))
                # Compute theta_e for ERA5
                cs["theta_e"] = compute_equiv_potential_temp_era5(cs)
                # Rename pressure_level to plev for consistency in plotting
                cs = cs.rename({'pressure_level': 'plev'})
            else:
                cs = cs.assign_coords(r=("lon", r))
                # Compute theta_e for RegCM5
                cs["theta_e"] = compute_equiv_potential_temp_regcm5(cs)
            
            theta_e_data[exp_name] = cs
            
            vmin_all = min(vmin_all, cs["theta_e"].min().values)
            vmax_all = max(vmax_all, cs["theta_e"].max().values)
            
        except Exception as e:
            print(f"Warning: Could not process {exp_name}: {e}")
            continue
    
    # Second pass: plot
    subplot_labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)']
    
    plot_idx = 0
    for idx, (exp_name, ds, center) in enumerate(all_datasets):
        if exp_name not in theta_e_data:
            continue
            
        ax = axes[plot_idx]
        cs = theta_e_data[exp_name]
        mtime = center['time']
        mlat = center['lat']
        mlon = center['lon']
        
        X, Y = np.meshgrid(cs.r.values, cs.plev.values)
        z = cs["theta_e"].values
        
        # Plot filled contour
        cf = ax.contourf(X, Y, z, levels=30, cmap="jet", vmin=vmin_all, vmax=vmax_all)
        
        # Add contour lines on top
        contour_levels = np.linspace(vmin_all, vmax_all, 15)
        ax.contour(X, Y, z, levels=contour_levels, colors='black', linewidths=0.5, alpha=0.5)
        
        # Add subplot label
        ax.text(0.02, 0.98, subplot_labels[plot_idx], transform=ax.transAxes, 
                fontsize=14, fontweight='bold', va='top', ha='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Title with experiment name and coordinates
        time_str = pd.to_datetime(mtime).strftime('%d/%m %H:%M')
        ax.set_title(f"{exp_name}\n({mlat:.2f}°, {mlon:.2f}°) {time_str}")
        ax.set_xlabel("R [km]")
        ax.set_ylabel("Pressure (hPa)")
        ax.yaxis.set_inverted(True)
        ax.grid(True, alpha=0.3, linestyle='--')
        
        plot_idx += 1
    
    # Hide unused subplots
    for idx in range(plot_idx, len(axes)):
        axes[idx].set_visible(False)
    
    # Colorbar
    if plot_idx > 0:
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(cf, cax=cbar_ax)
        cbar.set_label('Equivalent Potential Temperature (K)', labelpad=15, fontsize=12)
    
    plt.tight_layout(rect=[0, 0, 0.9, 1])
    
    if figname:
        os.makedirs(OUTPUT_PATH, exist_ok=True)
        plt.savefig(OUTPUT_PATH + figname, dpi=300, bbox_inches='tight')
        print(f"Saved: {OUTPUT_PATH + figname}")
    
    return fig, axes

# %%
# MAIN
if __name__ == "__main__":
    # Configuration
    experiments = ["ctrl", "holt_r2", "holt_r3", "uw_r2", "uw_r3"]
    map_extent = [-103, -92.5, 5, 17]  # [lon_min, lon_max, lat_min, lat_max]
    data_dir = f"/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/exps_v2/domain_{domain}_regridded/"
    era5_dir = "/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/era5"
    
    # Load ERA5 data (cropped to hurricane domain)
    print("Loading ERA5 data...")
    ds_era = load_era5_vts(era5_dir, map_extent)
    
    # Load RegCM5 models
    print("Loading RegCM5 models...")
    ds_models = {}
    ds_models['ERA5'] = ds_era  # Add ERA5 to models dictionary
    for exp in experiments:
        print(f"Loading {exp}...")
        ds_models[exp] = load_regcm5_multi_file(data_dir + exp)
    
    # Find cyclone centers for ALL datasets (including ERA5) using same function
    print("Finding cyclone centers for all datasets...")
    cyclone_centers, minimum_datetime = find_cyclone_centers(ds_models, map_extent)
    
    # Print centers for verification
    print("\nCyclone centers found:")
    for exp_name, center in minimum_datetime.items():
        if center['time'] is not None:
            print(f"  {exp_name}: time={center['time']}, lat={center['lat']:.2f}, lon={center['lon']:.2f}")
        else:
            print(f"  {exp_name}: No valid center found")
    
    # Plot VTS (ERA5 + 5 experiments = 6 subplots)
    print("\nPlotting VTS...")
    plot_vts_combined(ds_models, minimum_datetime, experiments, figname="pyplt_Hurricane_Otis_vts.png")
    plt.show()
    
    print("Done!")
