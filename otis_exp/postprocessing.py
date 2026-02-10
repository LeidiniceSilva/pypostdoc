# %%
import os
os.chdir("/mnt/shared/teamx/otis_exp")
import glob
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.ndimage import maximum_filter, minimum_filter
from datetime import datetime
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from metpy.units import units
import metpy.calc as mpcalc
from PIL import Image
import io
from tqdm import tqdm

COLORS = [
    "#33a02c", "#005B72","#905dc7", "#692510", 
    "#fdbf6f", "b", "#fb9a99", "#b2df8a",  "#a6cee3", "#ff7f00",  "#cab2d6", 
]
OUTPUT_PATH = "output/"

# %%

##########################################################################################
##########################################################################################
##############################    CYCLONE CENTER DETECTION    ##########################
##########################################################################################
##########################################################################################
def haversine_distance(lat1, lon1, lat2, lon2):
    """Compute great-circle distance (km) between two lat/lon points."""
    R = 6371.0
    phi1, phi2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlambda = np.radians(lon2 - lon1)
    a = np.sin(dphi/2)**2 + np.cos(phi1)*np.cos(phi2)*np.sin(dlambda/2)**2
    return 2 * R * np.arcsin(np.sqrt(a))


class CycloneTracker:
    """
    Class to track cyclone centers and compute associated meteorological variables.
    """
    
    def __init__(self, filepath, experiments, map_extent):
        """
        Initialize the cyclone tracker.
        
        Parameters:
        -----------
        filepath : str
            Base path to data files
        experiments : list
            List of experiment names
        map_extent : list
            [lon_min, lon_max, lat_min, lat_max] for analysis domain
        """
        self.data_dir = data_dir
        self.experiments = experiments
        self.map_extent = map_extent
        self.cyclone_centers = {}
        self.minimum_datetime = {}
        self.variable_data = {}
        
        # Load model datasets
        self.ds_models = {}
        for exp in experiments:
            try:
                self.ds_models[exp] = load_regcm5_multi_file(data_dir+exp)
            except FileNotFoundError:
                print(f"Warning: Could not load data for experiment {exp}")
        
        # Load ERA5 data
        try:
            era_dir = "data/data/ERA5"
            self.ds_era = load_era5(era_dir)
            self.ds_era = self.ds_era.rename({
            "valid_time": "time",
            "latitude": "lat",
            "longitude": "lon",
            "pressure_level": "plev",
            "msl": "psl",
            "t": "ta",
            "tp": "pr",
            "u10": "uas",
            "u": "ua",
            "v10": "vas",
            "v": "va",
            "w": "wa"
            })
            self.ds_era = self.fix_longitude(self.ds_era, 'lon')
            self.ds_models["ERA5"] = self.ds_era
            
        except FileNotFoundError:
            print("Warning: Could not load ERA5 data")
            self.ds_era = None

        
    def fix_longitude(self, ds, lon_var):
        """Fix longitude range to [-180, 180] if necessary."""
        if np.any(ds[lon_var] > 180):
            ds.coords[lon_var] = (ds.coords[lon_var] + 180) % 360 - 180
            ds = ds.sortby(ds[lon_var])
        return ds
    
    def find_cyclone_centers(self, psl_threshold=1010.0, max_jump_km=1000, smooth_window=None):
        """
        Find cyclone centers based on minimum sea level pressure,
        with robustness checks for dissipation and spurious jumps.
        
        Parameters
        ----------
        psl_threshold : float
            Minimum central pressure (hPa) to be considered a cyclone.
        max_jump_km : float
            Maximum allowed displacement per timestep (km).
        smooth_window : int or None
            Rolling window size for smoothing (timesteps). None disables smoothing.
        """
        print("Finding cyclone centers...")

        def process_dataset(ds, lat_name, lon_name, time_name, exp_name, max_jump_km):
            centers = []
            last_lat, last_lon = None, None

            for time_val in ds[time_name].values:
                ds_sel = ds.sel(
                    **{lat_name: slice(self.map_extent[3], self.map_extent[2]) if exp_name == "ERA5" else slice(self.map_extent[2], self.map_extent[3]),
                    lon_name: slice(self.map_extent[0], self.map_extent[1]),
                    time_name: time_val}
                )

                # Pressure in hPa
                pressure = ds_sel['psl']
                if pressure.attrs.get('units', '').lower() == 'pa':
                    pressure = pressure / 100.0
                
                if np.all(np.isnan(pressure.values)):
                    # raise ValueError("All data is NaN in the search region")
                    continue

                # Find min pressure
                min_pressure = pressure.min().item()
                min_point = pressure.where(pressure == min_pressure, drop=True)
                lat_center = min_point[lat_name].values.flatten()[0]
                lon_center = min_point[lon_name].values.flatten()[0]

                # --- Filters ---
                valid = True
                if min_pressure > psl_threshold:
                    valid = False
                if valid and last_lat is not None:
                    dist = haversine_distance(last_lat, last_lon, lat_center, lon_center)
                    if dist > max_jump_km:
                        valid = False

                if not valid:
                    lat_center, lon_center, min_pressure = np.nan, np.nan, np.nan

                centers.append({
                    'time': pd.to_datetime(time_val),
                    'lat': lat_center,
                    'lon': lon_center,
                    'min_pressure': min_pressure,
                    'experiment': exp_name
                })

                if valid:
                    last_lat, last_lon = lat_center, lon_center

            df = pd.DataFrame(centers)

            # --- Optional smoothing ---
            if smooth_window is not None and smooth_window > 1:
                df['lat'] = df['lat'].rolling(window=smooth_window, min_periods=1, center=True).mean()
                df['lon'] = df['lon'].rolling(window=smooth_window, min_periods=1, center=True).mean()
                df['min_pressure'] = df['min_pressure'].rolling(window=smooth_window, min_periods=1, center=True).mean()

            return df

        # Process models
        for exp_name, ds in self.ds_models.items():
            print(f"Processing experiment: {exp_name}")
            ds = self.fix_longitude(ds, 'lon')
            centers = process_dataset(ds, 'lat', 'lon', 'time', exp_name, max_jump_km=300)
            self.cyclone_centers[exp_name] = centers
            self.minimum_datetime[exp_name] = centers.iloc[centers["min_pressure"].idxmin()]
        # # Process ERA5
        # if self.ds_era is not None:
        #     print("Processing ERA5 data...")
        #     ds_era = self.fix_longitude(self.ds_era, 'longitude')
        #     self.cyclone_centers['ERA5'] = process_dataset(ds_era, 'latitude', 'longitude', 'valid_time', 'ERA5', max_jump_km=1000)

    
    def compute_area_statistics(self, var_name, stat_type='max', radius_deg=2.0):
        """
        Compute statistics for a variable in an area around the cyclone center.
        
        Parameters:
        -----------
        var_name : str
            Variable name ('sfcWind', 'pr', 'tas', 'psl')
        stat_type : str
            Type of statistic ('max', 'mean', 'min')
        radius_deg : float
            Radius around center in degrees
        """
        print(f"Computing {stat_type} {var_name} around cyclone centers...")
        models = self.ds_models
        # models["ERA5"] = self.ds_era
        for exp_name, ds in models.items():
            if exp_name not in self.cyclone_centers:
                continue
                
            centers_df = self.cyclone_centers[exp_name]
            results = []
            
            for _, center in centers_df.iterrows():
                if center.isna().any():
                    continue
                # Define area around center
                lat_min = center['lat'] - radius_deg
                lat_max = center['lat'] + radius_deg
                lon_min = center['lon'] - radius_deg
                lon_max = center['lon'] + radius_deg
                
                # Select data around center
                try:
                    ds_area = ds.sel(
                        lat=slice(lat_max, lat_min) if exp_name == "ERA5" else slice(lat_min, lat_max),
                        lon=slice(lon_min, lon_max),
                        time=center['time']
                    )
                    
                    var_data = ds_area[var_name]
                    
                
                    # Compute statistic
                    if stat_type == 'max':
                        stat_value = var_data.max().item()
                    elif stat_type == 'mean':
                        stat_value = var_data.mean().item()
                    elif stat_type == 'min':
                        stat_value = var_data.min().item()
                    else:
                        raise ValueError(f"Unknown stat_type: {stat_type}")
                    
                    results.append({
                        'time': center['time'],
                        'lat': center['lat'],
                        'lon': center['lon'],
                        f'{stat_type}_{var_name}': stat_value,
                        'experiment': exp_name
                    })
                    
                except Exception as e:
                    print(f"Warning: Could not compute {var_name} for {exp_name} at {center['time']}: {e}")
                    continue
                    
        
            # Store results
            key = f"{exp_name}-{stat_type}-{var_name}"
            self.variable_data[key] = pd.DataFrame(results)



    def get_combined_dataframe(self):
        """
        Get a combined dataframe with all cyclone centers and computed variables.
        
        Parameters:
        -----------
        include_vars : list or None
            List of variable keys to include. If None, includes all.
        """
        all_data = []
        
        for exp_name in self.cyclone_centers:
            base_df = self.cyclone_centers[exp_name].copy().set_index('time')
            idx = base_df.index
            # Add computed variables
        
            include_vars = [k for k in tracker.variable_data.keys() if k.split("-")[0] == exp_name]
            
            for var_key in include_vars:
                var_df = self.variable_data[var_key].copy()
                
                # Extract variable column (skip time, lat, lon, experiment)
                var_cols = [col for col in var_df.columns 
                        if col not in ['time', 'lat', 'lon', 'experiment']]
                
                for col in var_cols:
                    base_df[col] = var_df.set_index('time')[col].reindex(idx).copy()
            
            all_data.append(base_df.reset_index())
        
        return pd.concat(all_data, ignore_index=True)
    
    
##########################################################################################
##########################################################################################
##############################    PLOTTING FUNCTIONS    ################################
##########################################################################################
##########################################################################################

def plot_time_series_comparison(tracker, var_name, experiments=None, nhc_data=None, 
                               era5_included=True, ax=None, **kwargs):
    """
    Plot time series comparison of a variable across experiments.
    
    Parameters:
    -----------
    tracker : CycloneTracker
        The cyclone tracker object
    var_name : str
        Name of variable to plot
    experiments : list or None
        List of experiments to include
    nhc_data : DataFrame or None
        NHC best track data
    era5_included : bool
        Whether to include ERA5 data
    ax : matplotlib axis or None
        Axis to plot on
    **kwargs : dict
        Additional plotting arguments
    """
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
    
    if experiments is None:
        experiments = list(tracker.cyclone_centers.keys())
        if not era5_included and 'ERA5' in experiments:
            experiments.remove('ERA5')
    
    # Plot model experiments
    for i, exp in enumerate(experiments):
        if exp == 'ERA5':
            continue
            
        if exp in tracker.cyclone_centers:
            df = tracker.cyclone_centers[exp]
            
            if var_name in df.columns:
                color = COLORS[i % len(COLORS)]
                ax.plot(df['time'], df[var_name], 'o-', 
                       label=exp, color=color, alpha=0.8, **kwargs)
    
    # Plot ERA5 if requested
    if era5_included and 'ERA5' in tracker.cyclone_centers:
        df_era5 = tracker.cyclone_centers['ERA5']
        df_era5 = df_era5.drop(index=range(23, 28), errors="ignore")
        if var_name in df_era5.columns:
            ax.plot(df_era5['time'], df_era5[var_name], 'o-', 
                   label='ERA5', color='blue', alpha=0.8)
    
    # Plot NHC data if provided
    if nhc_data is not None and var_name in nhc_data.columns:
        ax.plot(nhc_data['time'], nhc_data[var_name], 'o-', 
               label='NHC', color="r", linewidth=2)
    
    #ax.set_xlabel('Time')
    ax.set_ylabel(var_name)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='x', rotation=45)
    
    return ax

def plot_cyclone_tracks(tracker, experiments=None, nhc_data=None, era5_included=True):
    """
    Plot cyclone tracks on a map.
    
    Parameters:
    -----------
    tracker : CycloneTracker
        The cyclone tracker object
    experiments : list or None
        List of experiments to include
    nhc_data : DataFrame or None
        NHC best track data
    era5_included : bool
        Whether to include ERA5 data
    """
    
    fig = plt.figure(figsize=(12, 10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    # Base map
    ax.add_feature(cfeature.COASTLINE, linewidth=0.7)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle=":")
    ax.add_feature(cfeature.LAND, facecolor="lightgray")
    ax.add_feature(cfeature.OCEAN, facecolor="lightblue")
    ax.set_extent(tracker.map_extent, crs=ccrs.PlateCarree())
    
    if experiments is None:
        experiments = [exp for exp in tracker.cyclone_centers.keys() if exp != 'ERA5']
    
    # Plot model tracks
    for i, exp in enumerate(experiments):
        if exp in tracker.cyclone_centers:
            df = tracker.cyclone_centers[exp]
            color = COLORS[i % len(COLORS)]
            
            # Plot track
            ax.plot(df['lon'], df['lat'], 'o-', color=color, 
                   linewidth=2, markersize=4, label=exp, alpha=0.8)
            
            # Mark start point
            ax.scatter(df['lon'].iloc[0], df['lat'].iloc[0], 
                      marker='*', s=150, color=color, 
                      edgecolor='black', zorder=5)
            
            # Add date labels for 00 UTC
            df_00utc = df[df['time'].dt.hour == 0]
            for _, row in df_00utc.iterrows():
                ax.text(row['lon'] + 0.3, row['lat'] + 0.3,
                       row['time'].strftime('%d/%m'),
                       fontsize=8, color=color, weight='bold')
    
    # Plot ERA5 track
    if era5_included and 'ERA5' in tracker.cyclone_centers:
        df_era5 = tracker.cyclone_centers['ERA5']
        df_era5 = df_era5.drop(index=range(23, 28), errors="ignore")
        
        ax.plot(df_era5['lon'], df_era5['lat'], 's-', color='blue',
               linewidth=2, label='ERA5', alpha=0.8)
        ax.scatter(df_era5['lon'].iloc[0], df_era5['lat'].iloc[0],
                  marker='*', s=150, color='blue', edgecolor='black', zorder=5)
        
        # Add date labels
        df_00utc = df_era5[df_era5['time'].dt.hour == 0]
        for _, row in df_00utc.iterrows():
            ax.text(row['lon'] + 0.3, row['lat'] + 0.3,
                   row['time'].strftime('%d/%m'),
                   fontsize=8, color='blue', weight='bold')
    
    # Plot NHC track
    if nhc_data is not None:
        ax.plot(nhc_data['lon'], nhc_data['lat'], 's-', color="r",
               linewidth=3, label='NHC', alpha=0.9)
        ax.scatter(nhc_data['lon'].iloc[0], nhc_data['lat'].iloc[0],
                  marker='*', s=150, color="r", edgecolor='black', zorder=5)
        
        # # Add date labels
        # nhc_00utc = nhc_data[nhc_data['time'].dt.hour == 0]
        # for _, row in nhc_00utc.iterrows():
        #     ax.text(row['lon'] + 0.3, row['lat'] + 0.3,
        #            row['time'].strftime('%d/%m'),
        #            fontsize=8, color="r", weight='bold')
    
    # Add reference point (Acapulco)
    acapulco_lon, acapulco_lat = -99.91, 16.85
    ax.plot(acapulco_lon, acapulco_lat, marker="*", color="red",
           markersize=12, zorder=10, label="Acapulco")
    
    # Formatting
    ax.gridlines(draw_labels=True, alpha=0.3)
    ax.legend(loc='best', fontsize=10)
    ax.set_title('Cyclone Tracks Comparison', fontsize=14, fontweight='bold')
    
    return fig, ax

# %%
def plot_max_sfc_wind(tracker, nhc_data=None):
    fig, ax = plt.subplots(figsize=(7, 4))
    for i, var in enumerate(tracker.variable_data.keys()):
        expname = var.split("-")[0]
        moist_data = tracker.variable_data[var]
        ax.plot(moist_data['time'], moist_data['max_sfcWind'], 
                            'o-', label=expname, alpha=0.8, color=COLORS[i])
    
    
    # Add NHC wind data
    if nhc_data is not None:
        ax.plot(nhc_data['time'], nhc_data['wind']* 0.514444, 
                    'o-', label='NHC', color="r", linewidth=2)

    # Define lines and labels (values in m/s)
    levels = {
        "TS": 18,
        "Cat1": 33,
        "Cat2": 43,
        "Cat3": 50,
        "Cat4": 58,
        "Cat5": 70
    }

    # Fill colors for intervals
    fill_colors = {
        "TS": "lightgreen",        # under TS
    "Cat1": "#FFFF00",    # Yellow
    "Cat2": "#FFD700",    # Gold / light orange
    "Cat3": "#FFA500",    # Orange
    "Cat4": "#FF4500",    # Orange-Red
    "Cat5": "#9B0000",    # Red
    "Above_Cat5": "#FF00FF"  # Bright Magenta
    }

    # Line colors for horizontal lines
    line_colors = {
        "TS": "#FFFF00",     # Yellow
    "Cat1": "#FFD700",   # Light Orange / Gold
    "Cat2": "#FFA500",   # Orange
    "Cat3": "#FF4500",   # Orange-Red
    "Cat4": "#FF0000",   # Red
    "Cat5": "#FF00FF"    # Magenta
    }
    
    x_min, x_max = ax.get_xlim()
    y_top = 85

    # Fill intervals
    intervals = [
        (0, levels["TS"], fill_colors["TS"], line_colors["TS"], "TS"),
        (levels["TS"], levels["Cat1"], fill_colors["Cat1"], line_colors["Cat1"], "Cat1"),
        (levels["Cat1"], levels["Cat2"], fill_colors["Cat2"], line_colors["Cat2"], "Cat2"),
        (levels["Cat2"], levels["Cat3"], fill_colors["Cat3"], line_colors["Cat3"], "Cat3"),
        (levels["Cat3"], levels["Cat4"], fill_colors["Cat4"], line_colors["Cat4"], "Cat4"),
        (levels["Cat4"], levels["Cat5"], fill_colors["Cat5"], line_colors["Cat5"], "Cat5"),
        (levels["Cat5"], y_top, fill_colors["Above_Cat5"], fill_colors["Above_Cat5"], "Above Cat5")
    ]

    for i, (y0, y1, fill_color, line_color, label) in enumerate(intervals):
        ax.fill_between([x_min, x_max], y0, y1, color=fill_color, alpha=0.4)
        if label != "Above Cat5":
            # Determine the interval above
            if i < len(intervals) - 1:
                next_y0, next_y1 = intervals[i + 1][0], intervals[i + 1][1]
                y_text = next_y0 + (next_y1 - next_y0)/2  # middle of the next interval
            else:
                y_text = y1 + (y_top - y1)/2  # for the last interval, place in the middle of top
            x_text = pd.to_datetime(tracker.variable_data['ERA5-max-sfcWind']["time"].max()) + pd.Timedelta(days=-0.5)
            ax.text(x_text, y_text, label, va='center', ha='left', fontsize=15, color='k')

    # Set axes limits and labels
    ax.set_ylim(0, y_top)
    #ax.set_xlabel("Date")
    ax.set_ylabel("Max wind speed (m/s)")
    ax.set_title("Wind Speed Comparison")
    #ax.grid(True)
    ax.margins(x=0)
    ax.legend()

    ax.tick_params(axis='x', rotation=45)
        
    return fig, ax

# %%

def load_nhc_data():
    """Load NHC best track data."""
    raw_data = [
        ("21/0000",  9.3,  95.8, 1007,  25),
        ("21/0600",  9.6,  96.1, 1007,  25),
        ("21/1200",  9.7,  96.6, 1007,  25),
        ("21/1800",  9.5,  96.9, 1007,  25),
        ("22/0000",  9.5,  96.7, 1007,  25),
        ("22/0600",  9.6,  96.6, 1006,  30),
        ("22/1200",  9.8,  96.7, 1006,  30),
        ("22/1800", 10.2,  96.9, 1004,  35),
        ("23/0000", 10.8,  97.1, 1004,  35),
        ("23/0600", 11.4,  97.2, 1004,  35),
        ("23/1200", 11.9,  97.3, 1003,  40),
        ("23/1800", 12.6,  97.6, 1001,  45),
        ("24/0000", 13.3,  98.0,  998,  50),
        ("24/0600", 13.8,  98.5,  994,  55),
        ("24/1200", 14.2,  98.9,  990,  65),
        ("24/1800", 14.9,  99.4,  971, 100),
        ("25/0000", 15.7,  99.6,  938, 130),
        ("25/0300", 16.1,  99.7,  922, 145),
        ("25/0600", 16.7,  99.8,  924, 145),
        ("25/0645", 16.8,  99.9,  929, 140),
        ("25/1200", 17.7, 100.3,  980,  85),
        ("25/1800", 18.6, 100.7, 1004,  40)
    ]
    
    def parse_date(datestr):
        day, hour = datestr.split("/")
        return datetime(2023, 10, int(day), int(hour[:2]), int(hour[2:]))
    
    df_nhc = pd.DataFrame(
        [(parse_date(d), lat, -lon, pres, wind) for d, lat, lon, pres, wind in raw_data],
        columns=["time", "lat", "lon", "min_pressure", "wind"]
    )
    
    return df_nhc

######### Vertical Thermal Structure ########
def compute_equiv_potential_temp(da, relative_hum=False):
    T = da['ta']       # K, shape (time, plev, lat, lon)
    p = da['plev']
    p_broadcasted = p.broadcast_like(T)
    if relative_hum:
        # Convert RH from % to unitless fraction and clip to [0, 1] (discard supersaturation in ERA5)
        rh = da['r'].clip(min=0.0, max=100)

        dpt = mpcalc.dewpoint_from_relative_humidity(
            T * units.K,
            rh * units.percent
        )
    else:
        q = da['hus'] 
        r = q / (1 - q)
        dpt = mpcalc.dewpoint_from_specific_humidity(p_broadcasted * units.hPa, q * units("kg/kg"))
    theta_e = mpcalc.equivalent_potential_temperature(p_broadcasted * units.hPa, T * units.K, dpt)
    return theta_e


def haversine_signed_1d_from_center(lat, lon):
    """
    Compute 1D signed distance (in km) from the central point (lat_c, lon_c)
    along a lat/lon cross-section (no meshgrid used).

    Parameters:
        lat: 1D array of latitude values (degrees)
        lon: 1D array of longitude values (degrees)

    Returns:
        radius_km: 1D array of signed distances (km), same length as lat/lon
    """
    R = 6371.0  # Earth radius in km

    # Center point
    idx_c = len(lat) // 2
    lat_c = lat[idx_c]
    lon_c = lon[idx_c]

    # Convert to radians
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)
    lat_c_rad = np.radians(lat_c)
    lon_c_rad = np.radians(lon_c)

    # Differences
    delta_lat = lat_rad - lat_c_rad
    delta_lon = lon_rad - lon_c_rad

    # Compute signed distances
    dy = R * delta_lat  # N/S component
    dx = R * np.cos((lat_rad + lat_c_rad) / 2) * delta_lon  # E/W component

    # Compute total distance and preserve sign based on direction
    # Sign is based on projection along the section (dot product with center-to-point vector)
    signs = np.sign(dy + dx)  # approximate directionality
    radius_km = signs * np.sqrt(dx**2 + dy**2)

    return radius_km


def plot_var(da, extent):
    fig, ax = plt.subplots(
        figsize=(10, 6),
        subplot_kw={'projection': ccrs.PlateCarree()}
    )
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')

    cs = da.plot.contour(
        ax=ax,
        transform=ccrs.PlateCarree(),
        colors='black',
        linewidths=0.5,
        levels=np.arange(int(da.min().values), int(da.max().values)+1, 5)
    )

    to_label = cs.levels[[-3,-2,-1]]
    ax.clabel(cs, levels=to_label, inline=True, fontsize=8, fmt='%1.0f', inline_spacing=10, use_clabeltext=None)
    ax.set_title(None)

    gl = ax.gridlines(draw_labels=True, linestyle='--', linewidth=0.5, color='gray')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10}
    gl.ylabel_style = {'size': 10}

    return fig, ax

# %%
############## Plot Vertical Thermal Structure #############
def plot_vts(tracker, figname=None):
    minimum_datetime = tracker.minimum_datetime
    for exp_name, ds in tracker.ds_models.items():
        mtime, mlat, mlon = minimum_datetime[exp_name]["time"], minimum_datetime[exp_name]["lat"],minimum_datetime[exp_name]["lon"]
        lons = ds.lon.values
        cs = ds.sel(time=mtime, lat=mlat, lon=lons[(lons > mlon - 3) & (lons < mlon + 3)], method="nearest")
        r = haversine_signed_1d_from_center(np.array([cs.coords['lat'].values.item()] * len(cs.coords['lon'].values)), cs.coords['lon'].values)
        cs = cs.assign_coords(r=("lon", r))
        if exp_name == "ERA5":
            continue
           # cs["theta_e"] = compute_equiv_potential_temp(cs, relative_hum=True)
        else: 
            cs["theta_e"] = compute_equiv_potential_temp(cs)
        fig, ax = plot_cross_section(cs)
        ax.set_title(f"{exp_name}: ({mlat:.4f}, {mlon:.4f}) @ {str(mtime)} {'(' + figname.split('_')[1]+ ')' if figname else ''}")
        if figname:
            fname = figname + '_' + exp_name
            plt.savefig(OUTPUT_PATH + fname, dpi=300)
        plt.show()

def plot_cross_section(cs):
    var_name = "theta_e"
    plt.rcParams.update({
        "font.size": 14,
        "axes.labelsize": 16,
        "axes.titlesize": 18,
        "legend.fontsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12
    })

    fig, ax = plt.subplots(figsize=(9,6))
    X, Y = np.meshgrid(cs.r.values, cs.plev.values)
    z = cs[var_name].values

    # filled contour
    cf = ax.contourf(X, Y, z, levels=30, cmap="inferno")

    # colorbar with label
    cbar = fig.colorbar(cf, ax=ax, orientation='vertical',)
    cbar.set_label('Equivalent Potential Temperature (K)', labelpad=20)  # customize label

    # axes labels
    ax.set_xlabel("R [km]")
    ax.set_ylabel("Pressure (hPa)")
    ax.yaxis.set_inverted(True)

    plt.tight_layout()
    return fig, ax
# %%


############## Dataset loader #####################
# %%
def load_era5(era_dir, verbose=False):
    nc_files = sorted(glob.glob(os.path.join(era_dir, '*6hr*.nc')))
    
    if not nc_files:
        raise ValueError(f"No NetCDF files found in {era_dir}")
    
    # Load all datasets
    datasets = []
    for file in nc_files:
        ds_temp = xr.open_dataset(file)
        datasets.append(ds_temp)
    
    # Find the intersection of time coordinates across all datasets
    time_coords = [set(ds.valid_time.values) for ds in datasets]
    common_times = sorted(time_coords[0].intersection(*time_coords[1:]))
    
    if verbose:
        print(f"\nFound {len(nc_files)} files")
        print(f"Time steps per file: {[len(t) for t in time_coords]}")
        print(f"Common time steps: {len(common_times)}")
    
    # Select only common time steps from each dataset
    datasets_aligned = [ds.sel(valid_time=common_times) for ds in datasets]
    
    # Merge the aligned datasets
    ds_merged = xr.merge(datasets_aligned, compat='override')
    
    if verbose:
        print(f"\nMerged dataset:")
        print(f"  Variables: {list(ds_merged.data_vars)}")
        print(f"  Dimensions: {dict(ds_merged.dims)}")
        print(f"  Coordinates: {list(ds_merged.coords)}")
    
    return ds_merged

def load_regcm5_multi_file(data_dir, verbose=False):
    """
    Load RegCM5 data from separate variable files and merge into single dataset.
    Only uses the intersection of time steps across all files.
    
    Parameters
    ----------
    data_dir : str
        Directory containing separate NetCDF files for each variable
    verbose : bool, optional
        Print diagnostic information
    
    Returns
    -------
    xr.Dataset
        Merged dataset with all variables at intersecting time steps
    """
    nc_files = sorted(glob.glob(os.path.join(data_dir, '*.nc')))
    
    if not nc_files:
        raise ValueError(f"No NetCDF files found in {data_dir}")
    
    # Load all datasets
    datasets = []
    for file in nc_files:
        ds_temp = xr.open_dataset(file)
        datasets.append(ds_temp)
    
    # Find the intersection of time coordinates across all datasets
    time_coords = [set(ds.time.values) for ds in datasets]
    common_times = sorted(time_coords[0].intersection(*time_coords[1:]))
    
    if verbose:
        print(f"\nFound {len(nc_files)} files")
        print(f"Time steps per file: {[len(t) for t in time_coords]}")
        print(f"Common time steps: {len(common_times)}")
    
    # Select only common time steps from each dataset
    datasets_aligned = [ds.sel(time=common_times) for ds in datasets]
    
    # Merge the aligned datasets
    ds_merged = xr.merge(datasets_aligned, compat='override')
    
    if verbose:
        print(f"\nMerged dataset:")
        print(f"  Variables: {list(ds_merged.data_vars)}")
        print(f"  Dimensions: {dict(ds_merged.dims)}")
        print(f"  Coordinates: {list(ds_merged.coords)}")
    
    return ds_merged

# %%

if __name__ == "__main__":
    # %%
    # Configuration
    experiments = ["ctrl", "holt_r2", "holt_r3", "uw_r2", "uw_r3"]
    map_extent = [-103, -92.5, 5, 17]
    data_dir = "data/domain_large_regridded/"
    
    # Initialize tracker
    tracker = CycloneTracker(data_dir, experiments, map_extent)
    
    # Find cyclone centers
    tracker.find_cyclone_centers()
    
    # Compute additional variables
    tracker.compute_area_statistics('sfcWind', 'max', radius_deg=2)  # Max wind in 1.5° radius
    
    # Load NHC data
    nhc_data = load_nhc_data()
    # %%
    # Plot tracks
    fig, ax = plot_cyclone_tracks(tracker, experiments, nhc_data, era5_included=True)
    plt.savefig(OUTPUT_PATH + "track_large.png", dpi=300)
    plt.show()

    # %%
    plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 16,
    "axes.titlesize": 18,
    "legend.fontsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12
    })
    fig, ax = plt.subplots(figsize=(7, 6))
    
    # # Minimum pressure
    plot_time_series_comparison(tracker, 'min_pressure', experiments, 
                              nhc_data, ax=ax)
    ax.set_title('Minimum Sea Level Pressure')
    ax.set_ylabel('Pressure (hPa)')
    plt.savefig(OUTPUT_PATH + "msl_large.png", dpi=300)
    plt.show()

    # %%
    plot_max_sfc_wind(tracker, nhc_data)
    plt.savefig(OUTPUT_PATH + "wind_large.png", dpi=300)
    plt.show()
    
    # %%
    # Vertical Thermal structure
    plot_vts(tracker, figname="vts_large")






# %%
