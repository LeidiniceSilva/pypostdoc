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

def plot_precipitation_timeseries(
    ds, lat, lon, precip_var='pr',
    time_dim='time', lat_dim='lat', lon_dim='lon',
    ax=None, axis_settings=None, dataframe=False
):
    if dataframe:
        x = ds[precip_var]
        t = ds.index
    else:
        precip_point = ds[precip_var].sel({lat_dim: lat, lon_dim: lon}, method='nearest')
        t = precip_point[time_dim].values
        x = precip_point.values

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))

    markerline, stemlines, baseline = ax.stem(t, x, markerfmt='none')
    stemlines.set_linewidth(0.5)
    ax.set_ylim(bottom=0)
    ax.xaxis.set_tick_params(which="minor", bottom=False)
    ax.xaxis.set_tick_params(bottom=False)

    if axis_settings:
        if 'xlabel' in axis_settings:
            ax.set_xlabel(axis_settings['xlabel'])
        if 'ylabel' in axis_settings:
            ax.set_ylabel(axis_settings['ylabel'])
        if 'title' in axis_settings:
            ax.set_title(axis_settings['title'])
        if 'xlim' in axis_settings:
            ax.set_xlim(axis_settings['xlim'])
        if 'ylim' in axis_settings:
            ax.set_ylim(axis_settings['ylim'])

    if not dataframe:
        actual_lat = float(precip_point[lat_dim])
        actual_lon = float(precip_point[lon_dim])
        if not axis_settings or 'title' not in axis_settings:
            ax.set_title(f'Precipitation at ({actual_lat:.2f}°, {actual_lon:.2f}°)')

    return ax


def plot_precipitation_timeseries_comparison(
    ds_dict, station, lat, lon, precip_var='pr',
    time_dim='time', lat_dim='lat', lon_dim='lon',
    figsize=(12, 8), model_name="Model"
):
    """
    ds_dict keys: 'model', 'ERA5', 'CMORPH'
    """
    fig, axes = plt.subplots(4, 1, figsize=figsize, sharex=True)

    # Model
    plot_precipitation_timeseries(
        ds_dict['model'], lat, lon, precip_var=precip_var,
        time_dim=time_dim, lat_dim=lat_dim, lon_dim=lon_dim,
        ax=axes[0],
        axis_settings={'ylabel': 'Precip [mm/hr]', 'title': model_name}
    )

    plot_precipitation_timeseries(
        station, lat, lon, precip_var=precip_var,
        time_dim=time_dim, lat_dim=lat_dim, lon_dim=lon_dim,
        ax=axes[0],
        axis_settings={'ylabel': 'Precip [mm/hr]', 'title': 'Observations'},
        dataframe=True
    )

    # ERA5
    plot_precipitation_timeseries(
        ds_dict['ERA5'], lat, lon, precip_var=precip_var,
        time_dim=time_dim, lat_dim=lat_dim, lon_dim=lon_dim,
        ax=axes[1],
        axis_settings={'ylabel': 'Precip [mm/hr]', 'title': 'ERA5'}
    )

    # CMORPH
    plot_precipitation_timeseries(
        ds_dict['CMORPH'], lat, lon, precip_var=precip_var,
        time_dim=time_dim, lat_dim=lat_dim, lon_dim=lon_dim,
        ax=axes[2],
        axis_settings={'ylabel': 'Precip [mm/hr]', 'title': 'CMORPH', 'xlabel': 'Time'}
    )

    plt.tight_layout()
    return fig, axes


def create_precipitation_gif_comparison(
    ds_dict, output_path, precip_var='pr',
    time_dim='time', lat_dim='lat', lon_dim='lon',
    extent=None, vmin=None, vmax=None,
    cmap='viridis', duration=100,
    figsize=(18, 6), dpi=100,
    add_colorbar=True, coastlines=True,
    time_slice=None, skip_frames=1,
    model_name="Model"
):
    """
    ds_dict keys: 'model', 'ERA5', 'CMORPH'
    """
    # Extract data
    data_model = ds_dict['model'][precip_var]
    data_era5 = ds_dict['ERA5'][precip_var]
    data_cmorph = ds_dict['CMORPH'][precip_var]

    # **NEW: Align to common time period (intersection)**
    common_times = (
        set(data_model[time_dim].values) &
        set(data_era5[time_dim].values) &
        set(data_cmorph[time_dim].values)
    )
    common_times = sorted(common_times)
    
    if len(common_times) == 0:
        print(f"Warning: No overlapping times found!")
        return
    
    print(f"Using {len(common_times)} common time steps")
    
    # Select only common times
    data_model = data_model.sel({time_dim: common_times})
    data_era5 = data_era5.sel({time_dim: common_times})
    data_cmorph = data_cmorph.sel({time_dim: common_times})

    # Spatial subset
    if extent is not None:
        lon_min, lon_max, lat_min, lat_max = extent
        sel_kwargs = {lat_dim: slice(lat_min, lat_max), lon_dim: slice(lon_min, lon_max)}
        data_model = data_model.sel(sel_kwargs)
        data_era5 = data_era5.sel(sel_kwargs)
        data_cmorph = data_cmorph.sel(sel_kwargs)

    # Time subset
    if time_slice is not None:
        data_model = data_model.isel({time_dim: time_slice})
        data_era5 = data_era5.isel({time_dim: time_slice})
        data_cmorph = data_cmorph.isel({time_dim: time_slice})

    # Skip frames
    if skip_frames > 1:
        slicer = {time_dim: slice(None, None, skip_frames)}
        data_model = data_model.isel(slicer)
        data_era5 = data_era5.isel(slicer)
        data_cmorph = data_cmorph.isel(slicer)

    # Shared vmin/vmax
    if vmin is None:
        vmin = min(float(data_model.min()), float(data_era5.min()), float(data_cmorph.min()))
    if vmax is None:
        vmax = max(float(data_model.max()), float(data_era5.max()), float(data_cmorph.max()))

    # Cartopy
    projection = None
    if coastlines:
        try:
            projection = ccrs.PlateCarree()
        except Exception:
            coastlines = False

    frames = []
    n_frames = data_model.sizes[time_dim]
    print(f"Generating {n_frames} frames for comparison GIF...")

    acapulco_lon, acapulco_lat = -99.91, 16.85

    for i in tqdm(range(n_frames)):
        if coastlines and projection:
            fig, axes = plt.subplots(
                1, 3, figsize=figsize, dpi=dpi,
                subplot_kw={'projection': projection}
            )
        else:
            fig, axes = plt.subplots(1, 3, figsize=figsize, dpi=dpi)

        da_m = data_model.isel({time_dim: i})
        da_e = data_era5.isel({time_dim: i})
        da_c = data_cmorph.isel({time_dim: i})
        time_val = da_m[time_dim].values

        ims = []
        for ax, da, title in zip(
            axes,
            [da_m, da_e, da_c],
            [model_name, "ERA5", "CMORPH"]
        ):
            im = da.plot(
                ax=ax, vmin=vmin, vmax=vmax, cmap=cmap,
                add_colorbar=False, add_labels=False
            )
            ims.append(im)
            ax.set_title(f"{title}\n{pd.to_datetime(time_val)}", fontsize=10)
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")
            if coastlines and projection:
                ax.coastlines()
                ax.add_feature(cfeature.BORDERS, linestyle=':')
            ax.plot(acapulco_lon, acapulco_lat, marker="*", color="red",
                    markersize=10, zorder=10, alpha=0.7)

        if add_colorbar:
            fig.subplots_adjust(right=0.9)
            cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
            cbar = fig.colorbar(ims[-1], cax=cbar_ax)
            cbar.set_label("Precipitation [mm/hr]", fontsize=10)

        plt.tight_layout(rect=[0, 0, 0.9, 1])

        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=dpi, bbox_inches='tight')
        buf.seek(0)
        frames.append(Image.open(buf).copy())
        buf.close()
        plt.close(fig)

    print(f"Saving comparison GIF to {output_path}...")
    frames[0].save(
        output_path,
        save_all=True,
        append_images=frames[1:],
        duration=duration,
        loop=0,
        optimize=False
    )
    print(f"Comparison GIF saved successfully! ({len(frames)} frames)")


def fix_longitude(ds, lon_var):
        """Fix longitude range to [-180, 180] if necessary."""
        if np.any(ds[lon_var] > 180):
            ds.coords[lon_var] = (ds.coords[lon_var] + 180) % 360 - 180
            ds = ds.sortby(ds[lon_var])
        return ds

# %% load reference datasets (ctrl, ERA5, CMORPH as in your code)
def plot_precipitation_timeseries_comparison(
    ds_dict, lat, lon, precip_var='pr',
    time_dim='time', lat_dim='lat', lon_dim='lon',
    figsize=(12, 8), model_name="Model"
):
    """
    ds_dict keys: 'model', 'ERA5', 'CMORPH'
    """
    # **NEW: Align to common time period**
    data_model = ds_dict['model'][precip_var]
    data_era5 = ds_dict['ERA5'][precip_var]
    data_cmorph = ds_dict['CMORPH'][precip_var]
    station_df = ds_dict['observations']
    
    common_times = (
        set(data_model[time_dim].values) &
        set(data_era5[time_dim].values) &
        set(data_cmorph[time_dim].values)
    )
    common_times = sorted(common_times)
    
    # Create aligned datasets
    ds_aligned = {
        'model': ds_dict['model'].sel({time_dim: common_times}),
        'ERA5': ds_dict['ERA5'].sel({time_dim: common_times}),
        'CMORPH': ds_dict['CMORPH'].sel({time_dim: common_times})
    }

    station_aligned = station.reindex(common_times)
    
    fig, axes = plt.subplots(4, 1, figsize=figsize, sharex=True)

    # Model
    plot_precipitation_timeseries(
        ds_aligned['model'], lat, lon, precip_var=precip_var,
        time_dim=time_dim, lat_dim=lat_dim, lon_dim=lon_dim,
        ax=axes[0],
        axis_settings={'ylabel': 'Precip [mm/hr]', 'title': model_name}
    )
    plot_precipitation_timeseries(
        station_aligned, lat, lon, precip_var=precip_var,
        time_dim=time_dim, lat_dim=lat_dim, lon_dim=lon_dim,
        ax=axes[1],
        axis_settings={'ylabel': 'Precip [mm/hr]', 'title': 'Observation'},
        dataframe=True
    )

    # ERA5
    plot_precipitation_timeseries(
        ds_aligned['ERA5'], lat, lon, precip_var=precip_var,
        time_dim=time_dim, lat_dim=lat_dim, lon_dim=lon_dim,
        ax=axes[2],
        axis_settings={'ylabel': 'Precip [mm/hr]', 'title': 'ERA5'}
    )

    # CMORPH
    plot_precipitation_timeseries(
        ds_aligned['CMORPH'], lat, lon, precip_var=precip_var,
        time_dim=time_dim, lat_dim=lat_dim, lon_dim=lon_dim,
        ax=axes[3],
        axis_settings={'ylabel': 'Precip [mm/hr]', 'title': 'CMORPH', 'xlabel': 'Time'}
    )

    plt.tight_layout()
    return fig, axes



# %%
station = pd.read_csv("data/data/obs_pr_gauge/acapulco.csv")
station = station[["Fecha UTC", "Precipitacion (mm)"]].rename(columns={"Fecha UTC": "time", "Precipitacion (mm)": "pr"})
station["time"] = pd.to_datetime(station["time"])
station = station.set_index("time")
station = station.resample("h").sum()

# %%
experiments = ["ctrl", "holt_r2", "holt_r3", "uw_r2", "uw_r3",  "ERA5", "CMORPH"]
data_dir = "data/domain_large_regridded/"
ds_models = {}

for exp in experiments:
    if exp == "ERA5":
        ds_models[exp] = xr.open_dataset("data/data/ERA5/tp_ERA5_Otis_1hr_6hr_Oct2023.nc")
        ds_models[exp] = ds_models[exp].rename({
            "valid_time": "time",
            "latitude": "lat",
            "longitude": "lon",
            "tp": "pr",
        })
        ds_models[exp] = fix_longitude(ds_models[exp], "lon")
        ds_models[exp] = ds_models[exp].sortby(ds_models[exp]['lat'])
        ds_models[exp]["pr"] = ds_models[exp]["pr"] * 1000  # m -> mm
    elif exp == "CMORPH":
        ds_models[exp] = xr.open_dataset(
            "data/data/satellite_pr/CMORPH_cut_domain_20251022_20251025.nc"
        )
        ds_models[exp] = ds_models[exp].rename({"cmorph": "pr"})
        ds_models[exp] = fix_longitude(ds_models[exp], "lon")
    else:  # ctrl
        ds_models[exp] = xr.open_dataset(data_dir + exp + f"/pr_{exp}_2023101900.nc")
        ds_models[exp]["pr"] = ds_models[exp]["pr"] * 3600  # kg m-2 s-1 -> mm/hr

    ds_models[exp] = ds_models[exp].sel(
        time=slice("2023-10-19T01:00:00", "2023-10-27T00:00:00")
    )

map_extent = (-103, -92.5, 5, 17)
acapulco_lon, acapulco_lat = -99.91, 16.85

# %% process model experiments with comparison
model_experiments = ["ctrl", "holt_r2", "holt_r3", "uw_r2", "uw_r3"]

for model_exp in model_experiments:
    print(f"\nProcessing {model_exp}...")

    if model_exp not in ds_models:
        ds_models[model_exp] = xr.open_dataset(
            data_dir + model_exp + f"/pr_{model_exp}_2023101900.nc"
        )
        ds_models[model_exp]["pr"] = ds_models[model_exp]["pr"] * 3600
        ds_models[model_exp] = ds_models[model_exp].sel(
            time=slice("2023-10-19T01:00:00", "2023-10-27T00:00:00")
        )

    ds_comparison = {
        "model": ds_models[model_exp],
        "ERA5": ds_models["ERA5"],
        "CMORPH": ds_models["CMORPH"],
        'observations': station
    }

    # 3-row timeseries
    fig, axes = plot_precipitation_timeseries_comparison(
        ds_comparison, acapulco_lat, acapulco_lon,
        figsize=(12, 8), model_name=model_exp.upper()
    )
    plt.savefig(OUTPUT_PATH + f"timeseries_comparison_{model_exp}.png",
                dpi=150, bbox_inches='tight')
    plt.close(fig)

    # 3-column GIF
    create_precipitation_gif_comparison(
        ds_comparison,
        OUTPUT_PATH + f"gif_comparison_{model_exp}.gif",
        precip_var='pr',
        time_dim='time', lat_dim='lat', lon_dim='lon',
        extent=map_extent,
        cmap='viridis',
        duration=100,
        figsize=(18, 6),
        dpi=100,
        add_colorbar=True,
        coastlines=True,
        time_slice=None,
        skip_frames=1,
        model_name=model_exp.upper()
    )

print("\nAll comparisons completed!")

# %%
