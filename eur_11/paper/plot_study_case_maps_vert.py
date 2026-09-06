# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Sep 06, 2026"
__description__ = "This script computes and plots vertical profiles"

import os, glob, pathlib, warnings
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')

# Configuration and Paths
PATH_BASE = pathlib.Path("/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11")
DATA_DIR_REGCM = PATH_BASE / "postproc/paper/cyc/RegCM5"
DATA_DIR_ERA5  = PATH_BASE / "postproc/paper/cyc/ERA5"
OUTPUT_DIR     = PATH_BASE / "figs/paper"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DT = "2006sep26"
EXPERIMENTS = ["ERA5", "NoTo-EUR", "WSM5-EUR", "WSM7-EUR", "WDM7-EUR"]
MAP_EXTENT = [-10, 30, 30, 60]  # [lon_min, lon_max, lat_min, lat_max]
BOX_SIZE = 1.5                  # 1.5 deg x 1.5 deg spatial region around center


def find_file(data_dir, var_candidates, exp_or_tag, dt):
    """Search for a single variable file matching candidate prefixes."""
    for var in var_candidates:
        pattern = f"{var}_{exp_or_tag}_6hr_{dt}.nc"
        file_path = data_dir / pattern
        if file_path.exists():
            return file_path
        
        # Fallback wildcard match
        matches = list(data_dir.glob(f"*{var}*6hr*{dt}.nc"))
        if matches:
            return matches[0]
            
    return None


def load_single_var_dataset(file_path):
    """Open a NetCDF file safely."""
    if file_path is None or not file_path.exists():
        return None
    return xr.open_dataset(file_path)


def load_era5_hydrometeors(data_dir, dt):
    """Load separated ERA5 variable files (ciwc, clwc, psl/msl) and merge them."""
    file_cli = find_file(data_dir, ['ciwc', 'cli'], 'ERA5', dt)
    file_clw = find_file(data_dir, ['clwc', 'clw'], 'ERA5', dt)
    file_psl = find_file(data_dir, ['psl', 'msl'], 'ERA5', dt)

    missing = []
    if not file_cli: missing.append('ciwc/cli')
    if not file_clw: missing.append('clwc/clw')
    if not file_psl: missing.append('psl/msl')

    if missing:
        raise FileNotFoundError(f"Missing required ERA5 files in {data_dir} for date {dt}: {missing}")

    found_files = [file_cli, file_clw, file_psl]
    print(f"Loading ERA5 files: {[f.name for f in found_files]}")

    ds_cli = xr.open_dataset(file_cli)
    ds_clw = xr.open_dataset(file_clw)
    ds_psl = xr.open_dataset(file_psl)

    ds_era5 = xr.merge([ds_cli, ds_clw, ds_psl], compat='override')

    if 'valid_time' in ds_era5.coords or 'valid_time' in ds_era5.variables:
        ds_era5 = ds_era5.rename({'valid_time': 'time'})

    rename_dict = {}
    for var_name, std_name in [('latitude', 'lat'), ('longitude', 'lon'), 
                               ('level', 'plev'), ('isobaricInhPa', 'plev'), 
                               ('pressure_level', 'plev'),
                               ('msl', 'psl'), ('ciwc', 'cli'), ('clwc', 'clw')]:
        if var_name in ds_era5.variables or var_name in ds_era5.coords:
            rename_dict[var_name] = std_name
            
    if rename_dict:
        ds_era5 = ds_era5.rename(rename_dict)

    if 'lon' in ds_era5.coords:
        ds_era5 = ds_era5.assign_coords(lon=(((ds_era5['lon'] + 180) % 360) - 180))
        ds_era5 = ds_era5.sortby('lon')

    if 'plev' in ds_era5.coords and np.max(ds_era5['plev'].values) > 2000:
        ds_era5['plev'] = ds_era5['plev'] / 100.0

    return ds_era5


def load_regcm5_hydrometeors(data_dir, exp, dt):
    """Load separated RegCM5 variable files (cli, clw, psl) and merge them."""
    file_cli = find_file(data_dir, ['cli'], exp, dt)
    file_clw = find_file(data_dir, ['clw'], exp, dt)
    file_psl = find_file(data_dir, ['psl', 'msl'], exp, dt)

    missing = []
    if not file_cli: missing.append('cli')
    if not file_clw: missing.append('clw')
    if not file_psl: missing.append('psl/msl')

    if missing:
        raise FileNotFoundError(f"Missing required RegCM5 files for '{exp}' in {data_dir}: {missing}")

    found_files = [file_cli, file_clw, file_psl]
    print(f"Loading {len(found_files)} files for {exp}: {[f.name for f in found_files]}")

    ds_cli = xr.open_dataset(file_cli)
    ds_clw = xr.open_dataset(file_clw)
    ds_psl = xr.open_dataset(file_psl)

    datasets = [ds_cli, ds_clw, ds_psl]
    common_times = sorted(set.intersection(*[set(ds.time.values) for ds in datasets]))
    
    if not common_times:
        raise ValueError(f"No matching timestamps found across separated files for experiment '{exp}'.")

    ds_merged = xr.merge([ds.sel(time=common_times) for ds in datasets], compat='override')

    rename_dict = {}
    for k, v in [('xlat', 'lat'), ('xlon', 'lon'), ('msl', 'psl')]:
        if k in ds_merged.variables or k in ds_merged.coords:
            rename_dict[k] = v

    if rename_dict:
        ds_merged = ds_merged.rename(rename_dict)

    if 'plev' in ds_merged.coords and np.max(ds_merged['plev'].values) > 2000:
        ds_merged['plev'] = ds_merged['plev'] / 100.0

    if 'lon' in ds_merged.coords and np.any(ds_merged['lon'].values > 180):
        ds_merged['lon'] = xr.where(ds_merged['lon'] > 180, ds_merged['lon'] - 360, ds_merged['lon'])

    return ds_merged


def find_cyclone_center(ds_time, map_extent):
    """Locate minimum PSL in 1D/2D lat/lon fields."""
    lon_min, lon_max, lat_min, lat_max = map_extent

    if 'psl' not in ds_time:
        raise KeyError(f"'psl' variable missing in dataset. Available: {list(ds_time.data_vars)}")

    psl = ds_time['psl'].squeeze().values
    lat_grid = ds_time['lat'].squeeze().values
    lon_grid = ds_time['lon'].squeeze().values

    if lat_grid.ndim == 1 and lon_grid.ndim == 1:
        lon_grid, lat_grid = np.meshgrid(lon_grid, lat_grid)

    if np.nanmean(psl) > 10000:
        psl /= 100.0

    if psl.shape != lat_grid.shape:
        psl = psl.reshape(lat_grid.shape)

    spatial_mask = (lon_grid >= lon_min) & (lon_grid <= lon_max) & (lat_grid >= lat_min) & (lat_grid <= lat_max)
    if not np.any(spatial_mask):
        return None

    psl_domain = np.where(spatial_mask, psl, np.nan)
    if np.all(np.isnan(psl_domain)):
        return None

    min_idx = np.unravel_index(np.nanargmin(psl_domain), psl_domain.shape)
    return {
        'lat': float(lat_grid[min_idx]),
        'lon': float(lon_grid[min_idx]),
        'min_pressure': float(psl_domain[min_idx])
    }


def compute_area_mean_profile(ds, center_lat, center_lon, delta=0.75):
    """Average hydrometeor variables across all timesteps over a 1.5x1.5 deg box."""
    lat_min, lat_max = center_lat - delta, center_lat + delta
    lon_min, lon_max = center_lon - delta, center_lon + delta

    if ds['lat'].ndim == 1 and ds['lon'].ndim == 1:
        mask = (ds['lat'] >= lat_min) & (ds['lat'] <= lat_max) & \
               (ds['lon'] >= lon_min) & (ds['lon'] <= lon_max)
        ds_box = ds.where(mask, drop=True)
        spatial_dims = ('lat', 'lon')
    else:
        mask = (ds['lat'] >= lat_min) & (ds['lat'] <= lat_max) & \
               (ds['lon'] >= lon_min) & (ds['lon'] <= lon_max)
        ds_box = ds.where(mask, drop=True)
        spatial_dims = ds['lat'].dims

    cli_profile = ds_box['cli'].mean(dim=['time'] + list(spatial_dims)).values
    clw_profile = ds_box['clw'].mean(dim=['time'] + list(spatial_dims)).values
    plev = ds_box['plev'].values

    # Convert kg/kg to g/kg (10^-3 kg/kg)
    cli_profile *= 1000000.0
    clw_profile *= 1000000.0

    return plev, cli_profile, clw_profile


def plot_vertical_profiles(profiles_dict, output_file=None):
    """Plot separate vertical profiles of hydrometeor mixing ratios in two subplots."""
    
    # Create 1 row, 2 columns with shared Y axis
    fig, axes = plt.subplots(1, 2, figsize=(12, 7), sharey=True)
    fig.patch.set_facecolor('#E0E0E0')
    fig.patch.set_alpha(0.75)

    ax1, ax2 = axes

    colors = {
        'ERA5': 'black',
        'NoTo-EUR': 'red',
        'WSM5-EUR': 'blue',
        'WSM7-EUR': 'green',
        'WDM7-EUR': 'orange'
    }

    for exp_name, data in profiles_dict.items():
        plev = data['plev']
        color = colors.get(exp_name, 'gray')
        label_title = exp_name if exp_name == "ERA5" else exp_name.split('-')[0]

        # Plot (a) Cloud Liquid Water
        ax1.plot(data['clw'], plev, color=color, linewidth=1, label=label_title)
        
        # Plot (b) Cloud Ice Water
        ax2.plot(data['cli'], plev, color=color, linewidth=1, label=label_title)

    # Formating for Subplot (a) - Liquid Water
    ax1.set_ylim(1000, 100)  # Inverted pressure axis (Surface to Upper Troposphere)
    ax1.set_xlim(0, 80)
    ax1.set_xlabel('Cloud liquid water (mg kg$^{-1}$)', fontsize=10)
    ax1.set_ylabel('Pressure (hPa)', fontsize=10)
    ax1.set_title('(a)', fontsize=10, fontweight='bold', loc='left')
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(loc='upper right', fontsize=10, frameon=True)

    # Formating for Subplot (b) - Ice Water
    ax2.set_xlabel('Cloud liquid ice (mg kg$^{-1}$)', fontsize=10)
    ax2.set_ylim(1000, 100)
    ax2.set_xlim(0, 120)
    ax2.set_title('(b)', fontsize=10, fontweight='bold', loc='left')
    ax2.grid(True, linestyle='--', alpha=0.6)

    plt.tight_layout(rect=[0, 0, 1, 0.95]) # Adjust layout to fit main title
    
    if output_file:
        plt.savefig(output_file, dpi=400, bbox_inches='tight')
    plt.show()

def main():
    print("=" * 80 + "\nVertical Profiles for ERA5 and RegCM5 Experiments (6-hr)\n" + "=" * 80)
    all_profiles = {}

    for exp in EXPERIMENTS:
        print(f"\nProcessing {exp}...")
        if exp == "ERA5":
            ds_model = load_era5_hydrometeors(DATA_DIR_ERA5, DT)
        else:
            ds_model = load_regcm5_hydrometeors(DATA_DIR_REGCM, exp, DT)

        # Compute mean on psl DataArray specifically
        ds_mean_psl = ds_model[['psl']].mean(dim='time')
        center = find_cyclone_center(ds_mean_psl, MAP_EXTENT)
        
        if not center:
            raise ValueError(f"Could not track cyclone center for experiment {exp}")

        print(f"{exp:10} | Center Lat: {center['lat']:.2f}°N | Center Lon: {center['lon']:.2f}°E")

        # Compute area-averaged profiles across 6-hr timesteps over 1.5x1.5 deg box
        plev, cli_prof, clw_prof = compute_area_mean_profile(ds_model, center['lat'], center['lon'], delta=BOX_SIZE / 2.0)

        all_profiles[exp] = {
            'plev': plev,
            'cli': cli_prof,
            'clw': clw_prof
        }

    output_file = OUTPUT_DIR / f'pyplt_maps_study_case_RegCM5_EUR-11_vert_{DT}.png'
    plot_vertical_profiles(all_profiles, output_file)

if __name__ == "__main__":
    main()
