# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Sep 05, 2026"
__description__ = "This script computes and plots vte cross-sections"

import os, glob, pathlib, warnings
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.units import units

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


def load_era5_data(data_dir, dt):
    """Load and process ERA5 NetCDF files, standardizing time, longitude, and variable names."""
    files = [
        data_dir / f"msl_ERA5_1hr_{dt}.nc",
        data_dir / f"t_ERA5_1hr_{dt}.nc",
        data_dir / f"q_ERA5_1hr_{dt}.nc",
        data_dir / f"r_ERA5_1hr_{dt}.nc"
    ]
    
    found_files = [f for f in files if f.exists()]
    if not found_files:
        found_files = sorted(glob.glob(str(data_dir / f"*ERA5*_{dt}.nc")))
        
    if not found_files:
        raise FileNotFoundError(f"No valid ERA5 files found in '{data_dir}' for date {dt}")

    print(f"Loading {len(found_files)} ERA5 files: {[os.path.basename(f) for f in found_files]}")
    datasets = [xr.open_dataset(f) for f in found_files]
    ds_era5 = xr.merge(datasets, compat='override')

    # Standardize time dimension name
    if 'valid_time' in ds_era5.coords or 'valid_time' in ds_era5.variables:
        ds_era5 = ds_era5.rename({'valid_time': 'time'})

    # Standardize ERA5 Variable and Coordinate Names (msl -> psl, q -> hus, t -> ta)
    rename_dict = {}
    for var_name, std_name in [('latitude', 'lat'), ('longitude', 'lon'), 
                               ('level', 'plev'), ('isobaricInhPa', 'plev'), 
                               ('pressure_level', 'plev'),
                               ('t', 'ta'), ('q', 'hus'), ('msl', 'psl')]:
        if var_name in ds_era5.variables or var_name in ds_era5.coords:
            rename_dict[var_name] = std_name
            
    if rename_dict:
        ds_era5 = ds_era5.rename(rename_dict)

    # Convert longitude from [0, 360) to [-180, 180) and sort array along lon dimension
    if 'lon' in ds_era5.coords:
        ds_era5 = ds_era5.assign_coords(lon=(((ds_era5['lon'] + 180) % 360) - 180))
        ds_era5 = ds_era5.sortby('lon')

    # Convert pressure levels to hPa if given in Pa
    if 'plev' in ds_era5.coords and np.max(ds_era5['plev'].values) > 2000:
        ds_era5['plev'] = ds_era5['plev'] / 100.0

    return ds_era5


def load_regcm5_data(data_dir, exp, dt):
    """Load and merge RegCM5 hourly variables for experiment `exp`."""

    patterns = [
        os.path.join(data_dir, f"*_{exp}_1hr_{dt}.nc"),
        os.path.join(data_dir, "**", f"*_{exp}_1hr_{dt}.nc"),
        os.path.join(data_dir, f"*_{exp}_{dt}.nc"),
        os.path.join(data_dir, "**", f"*_{exp}_{dt}.nc")
    ]
    
    nc_files = []
    for pattern in patterns:
        found = [f for f in sorted(glob.glob(pattern, recursive=True)) if "_day_" not in os.path.basename(f)]
        if found:
            nc_files = found
            break

    if not nc_files:
        raise FileNotFoundError(f"No valid NetCDF files found for experiment '{exp}' under '{data_dir}'")

    print(f"Found {len(nc_files)} files for {exp}: {[os.path.basename(f) for f in nc_files]}")
    datasets = [xr.open_dataset(f) for f in nc_files]

    # Intersect common timestamps across all variable files
    common_times = sorted(set.intersection(*[set(ds.time.values) for ds in datasets]))
    if not common_times:
        raise ValueError(f"No matching timestamps found across files for experiment '{exp}'.")

    ds_merged = xr.merge([ds.sel(time=common_times) for ds in datasets], compat='override')

    # Standardize RegCM coordinates
    rename_dict = {k: v for k, v in [('xlat', 'lat'), ('xlon', 'lon')] if k in ds_merged.variables or k in ds_merged.coords}
    if rename_dict:
        ds_merged = ds_merged.rename(rename_dict)

    # Convert pressure levels to hPa if given in Pa
    if 'plev' in ds_merged.coords and np.max(ds_merged['plev'].values) > 2000:
        ds_merged['plev'] = ds_merged['plev'] / 100.0

    # Standardize longitudes to [-180, 180]
    if 'lon' in ds_merged.coords and np.any(ds_merged['lon'].values > 180):
        ds_merged['lon'] = xr.where(ds_merged['lon'] > 180, ds_merged['lon'] - 360, ds_merged['lon'])

    return ds_merged


def find_cyclone_center(ds_time, map_extent):
    """Locate minimum MSLP in 1D/2D curvilinear or rectilinear lat/lon fields."""

    lon_min, lon_max, lat_min, lat_max = map_extent

    # Check for pressure variable name flexible fallback
    psl_var = 'psl' if 'psl' in ds_time else ('msl' if 'msl' in ds_time else None)
    if psl_var is None:
        raise KeyError(f"Neither 'psl' nor 'msl' variable found in dataset. Available: {list(ds_time.data_vars)}")

    psl = ds_time[psl_var].squeeze().values
    lat_grid = ds_time['lat'].squeeze().values
    lon_grid = ds_time['lon'].squeeze().values

    if 0 in (psl.size, lat_grid.size, lon_grid.size):
        return None

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


def track_cyclone(ds, map_extent, exp_name):
    """Track cyclone over all time steps to locate peak intensity."""

    track_data = []
    for t in ds.time.values:
        center = find_cyclone_center(ds.sel(time=t), map_extent)
        if center:
            track_data.append({'time': pd.to_datetime(t), **center})

    if not track_data:
        lat_sample, lon_sample = ds['lat'].squeeze().values, ds['lon'].squeeze().values
        print(f"\n[DIAGNOSTIC] Lat range: [{np.nanmin(lat_sample):.2f}, {np.nanmax(lat_sample):.2f}] | Lon range: [{np.nanmin(lon_sample):.2f}, {np.nanmax(lon_sample):.2f}]")
        raise ValueError(f"No cyclone center detected for experiment '{exp_name}'. Check MAP_EXTENT bounds.")

    df_track = pd.DataFrame(track_data)
    peak = df_track.loc[df_track['min_pressure'].idxmin()]
    print(f"{exp_name:10} | Peak Time: {peak['time']} | Lat: {peak['lat']:.2f}°N | Lon: {peak['lon']:.2f}°E | Min MSLP: {peak['min_pressure']:.1f} hPa")
    return peak


def extract_cross_section(ds, center_time, center_lat, center_lon, radius_deg=5):
    """Extract zonal cross-section through cyclone center for 1D or 2D grids."""

    ds_time = ds.sel(time=center_time, method='nearest')
    lon_min, lon_max = center_lon - radius_deg, center_lon + radius_deg
    lat_vals, lon_vals = ds_time['lat'].squeeze().values, ds_time['lon'].squeeze().values

    if lat_vals.ndim == 1:
        lat_idx = np.abs(lat_vals - center_lat).argmin()
        ds_cross = ds_time.isel(lat=lat_idx)
        lon_dim = 'lon' if 'lon' in ds_cross.dims else list(ds_cross.dims)[-1]
        lon_line = np.asarray(ds_cross[lon_dim].values).squeeze()
        lon_indices = np.where((lon_line >= lon_min) & (lon_line <= lon_max))[0]
        
        if len(lon_indices) == 0:
            raise ValueError(f"No longitude points found within ±{radius_deg}° of center.")

        lon_cross = lon_line[lon_indices]
        lat_cross = np.full(len(lon_cross), lat_vals[lat_idx])
        ds_cross = ds_cross.isel({lon_dim: lon_indices}).rename({lon_dim: 'distance_dim'})

    else:
        lat_dim, lon_dim = ds_time['lat'].dims[0], ds_time['lon'].dims[1]
        lat_idx = np.unravel_index(np.nanargmin(np.abs(lat_vals - center_lat)), lat_vals.shape)[0]

        lon_line = np.asarray(ds_time['lon'].isel({lat_dim: lat_idx}).values).squeeze()
        lat_line = np.asarray(ds_time['lat'].isel({lat_dim: lat_idx}).values).squeeze()
        lon_indices = np.where((lon_line >= lon_min) & (lon_line <= lon_max))[0]

        if len(lon_indices) == 0:
            raise ValueError(f"No longitude points found within ±{radius_deg}° of center.")

        ds_cross = ds_time.isel({lat_dim: lat_idx, lon_dim: lon_indices}).rename({lon_dim: 'distance_dim'})
        lon_cross, lat_cross = lon_line[lon_indices], lat_line[lon_indices]

    # Clean 1D Coordinates
    ds_cross = ds_cross.drop_vars([v for v in ['lat', 'lon'] if v in ds_cross.coords])
    ds_cross = ds_cross.assign_coords(lon=('distance_dim', lon_cross), lat=('distance_dim', lat_cross))

    # Calculate distance from cyclone center in km
    R = 6371.0
    distances = R * np.cos(np.radians(center_lat)) * (np.radians(lon_cross) - np.radians(center_lon))
    return ds_cross.assign_coords(r=('distance_dim', distances)).squeeze()


def compute_equivalent_potential_temperature(ds_slice):
    """Calculate theta_e using MetPy from temperature, specific humidity, and pressure."""

    T, p, q = ds_slice['ta'], ds_slice['plev'], ds_slice['hus']
    p_bcast = p.broadcast_like(T)

    T_u = T.values * units.kelvin
    p_u = p_bcast.values * units.hPa
    q_u = q.values * units('kg/kg')

    dewpoint = mpcalc.dewpoint_from_specific_humidity(p_u, T_u, q_u)
    theta_e = mpcalc.equivalent_potential_temperature(p_u, T_u, dewpoint)

    return xr.DataArray(
        theta_e.magnitude, coords=T.coords, dims=T.dims, name='theta_e'
    ).assign_coords(r=ds_slice.r)


def plot_vts_subplots(vts_data_dict, output_file=None):
    """Plot 2x3 grid cross-sections of theta_e including ERA5 and RegCM5 experiments."""

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes_flat = axes.flatten()
    levels = np.arange(310, 341, 1)
    labels = ['(a)', '(b)', '(c)', '(d)', '(e)']
    font_size, cf = 10, None

    for idx, (exp_name, data) in enumerate(vts_data_dict.items()):
        ax = axes_flat[idx]
        theta_e = data['theta_e']
        X, Y = np.meshgrid(theta_e.r.values, theta_e.plev.values)

        cf = ax.contourf(X, Y, np.squeeze(theta_e.values), levels=levels, cmap='jet', extend='both')
        cs = ax.contour(X, Y, np.squeeze(theta_e.values), levels=levels[::2], colors='black', linewidths=0.5, alpha=0.6)
        ax.clabel(cs, inline=True, fontsize=6, fmt='%d')

        ax.set_ylim(1000, 200)
        ax.set_xlim(-500, 500)
        ax.axvline(x=0, color='red', linestyle='--', linewidth=1.2, alpha=0.8)
        ax.grid(True, linestyle='--', alpha=0.5)

        title_text = exp_name if exp_name == "ERA5" else exp_name.split('-')[0]
        ax.set_title(f"{labels[idx]} {title_text}", loc='left', fontsize=font_size, fontweight='bold')

        if idx in [0, 3]:
            ax.set_ylabel('Pressure (hPa)', fontsize=font_size)
        if idx >= 3:
            ax.set_xlabel('Distance from Center (km)', fontsize=font_size)

    # Hide unused 6th subplot frame
    if len(vts_data_dict) < len(axes_flat):
        fig.delaxes(axes_flat[-1])

    cbar_ax = fig.add_axes([0.15, 0.04, 0.70, 0.02])
    cbar = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('Equivalent potential temperature (K)', fontsize=font_size)

    plt.tight_layout(rect=[0, 0.08, 1, 1])
    if output_file:
        plt.savefig(output_file, dpi=400, bbox_inches='tight')
    plt.show()


def main():
    print("=" * 80 + "\nComputing VTE for ERA5 and RegCM5 Experiments\n" + "=" * 80)
    all_vts_data = {}

    for exp in EXPERIMENTS:
        print(f"\nProcessing {exp}...")
        if exp == "ERA5":
            ds_model = load_era5_data(DATA_DIR_ERA5, DT)
        else:
            ds_model = load_regcm5_data(DATA_DIR_REGCM, exp, DT)

        peak_model = track_cyclone(ds_model, MAP_EXTENT, exp)
        ds_cross = extract_cross_section(ds_model, peak_model['time'], peak_model['lat'], peak_model['lon'], radius_deg=5)
        theta_e = compute_equivalent_potential_temperature(ds_cross)

        all_vts_data[exp] = {
            'theta_e': theta_e,
            'lat': peak_model['lat'],
            'lon': peak_model['lon'],
            'time': peak_model['time']
        }

    output_file = OUTPUT_DIR / f'pyplt_maps_study_case_ERA5_RegCM5_EUR-11_vte_{DT}.png'
    plot_vts_subplots(all_vts_data, output_file)

if __name__ == "__main__":
    main()

