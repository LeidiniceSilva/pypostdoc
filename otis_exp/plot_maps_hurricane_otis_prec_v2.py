# %% ONE FIGURE - All hourly precipitation comparisons (BAR PLOT)
import os
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

os.chdir("/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp")

OUTPUT_PATH = "/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/figs/"

# %% Helper function
def fix_longitude(ds, lon_var='lon'):
    if np.any(ds[lon_var] > 180):
        ds.coords[lon_var] = (ds.coords[lon_var] + 180) % 360 - 180
        ds = ds.sortby(ds[lon_var])
    return ds

# %% Load station data
station = pd.read_csv(
    "ws/Estacion_ACAPULCO_SEMAR_24_horas.csv",
    encoding="latin1",
    skiprows=8
)
station.columns = station.columns.str.strip()
station = station.loc[:, ~station.columns.str.contains("^Unnamed")]
station = station[["Fecha UTC", "Precipitación (mm)"]].rename(
    columns={"Fecha UTC": "time", "Precipitación (mm)": "pr"}
)
station["time"] = pd.to_datetime(station["time"])
station["pr"] = pd.to_numeric(station["pr"], errors="coerce")
station = station.set_index("time").resample("1h").sum()

# %% Load all datasets
experiments = ["noto", "wsm5", "wsm7", "wdm7", "ERA5", "CMORPH"]
ds_dict = {}
data_dir = "exps_mp_regrided/"

for exp in experiments:
    if exp == "ERA5":
        ds = xr.open_dataset("era5/1hr/tp_ERA5_Otis_1hr_Oct2023.nc")
        ds = ds.rename({"valid_time": "time", "latitude": "lat", "longitude": "lon", "tp": "pr"})
        ds = fix_longitude(ds, "lon")
        ds = ds.sortby(ds['lat'])
        ds["pr"] = ds["pr"] * 1000  # m -> mm
    elif exp == "CMORPH":
        ds = xr.open_dataset("cmorph/CMORPH_cut_domain_20251022_20251025.nc")
        ds = ds.rename({"cmorph": "pr"})
        ds = fix_longitude(ds, "lon")
    else:
        ds = xr.open_dataset(data_dir + exp + f"/pr_{exp}_2023101900.nc")
        ds["pr"] = ds["pr"] * 3600  # kg m-2 s-1 -> mm/hr
    
    ds = ds.sel(time=slice("2023-10-19T01:00:00", "2023-10-27T00:00:00"))
    ds_dict[exp] = ds

# %% Location
acapulco_lat, acapulco_lon = 16.85, -99.91

# %% Extract time series at Acapulco point
times_series = {}
for exp in experiments:
    point = ds_dict[exp]['pr'].sel(lat=acapulco_lat, lon=acapulco_lon, method='nearest')
    times_series[exp] = (point.time.values, point.values)

# Find common time range (intersection of all)
common_times = station.index
for exp in experiments:
    common_times = common_times.intersection(times_series[exp][0])

print(f"Plotting {len(common_times)} common time steps")

# %% PLOT - ONE FIGURE with BARS (grouped bar chart)
fig, ax = plt.subplots(figsize=(16, 6))

# Define bar width and positions
n_datasets = len(experiments) + 1  # +1 for station
bar_width = 0.12
x_positions = np.arange(len(common_times))

# Colors for each dataset
colors = {
    'noto': 'magenta',
    'wsm5': 'brown', 
    'wsm7': 'gray',
    'wdm7': 'green',
    'ERA5': 'red',
    'CMORPH': 'blue',
    'station': 'black'
}

# Plot each model as bars
for i, exp in enumerate(experiments):
    # Align to common times
    mask = np.isin(times_series[exp][0], common_times)
    values = times_series[exp][1][mask]
    
    # Offset position for this dataset
    offset = (i - n_datasets/2) * bar_width
    ax.bar(x_positions + offset, values, width=bar_width, 
           label=exp.upper(), color=colors[exp], alpha=0.8, edgecolor='black', linewidth=0.5)

# Plot station as bars
station_aligned = station.loc[common_times]
offset_station = (len(experiments) - n_datasets/2) * bar_width
ax.bar(x_positions + offset_station, station_aligned['pr'], 
       width=bar_width, label='Station', color=colors['station'], 
       alpha=0.8, edgecolor='black', linewidth=0.5)

# Formatting
ax.set_xlabel('Time', fontsize=10)
ax.set_ylabel('Precipitation (mm/hr)', fontsize=10)
ax.set_title(f'(a) Hourly Precipitation at Acapulco ({acapulco_lat:.2f}°, {acapulco_lon:.2}°)', fontsize=0)
ax.set_xticks(x_positions[::12])  # Show every 12th hour to avoid crowding
ax.set_xticklabels(pd.to_datetime(common_times[::12]).strftime('%m/%d %H:00'), rotation=45, ha='right')
ax.legend(loc='upper right', fontsize=10, ncol=7)
ax.grid(True, linestyle='--', axis='y')
ax.set_ylim(bottom=0)

plt.tight_layout()

# Save
plt.savefig(OUTPUT_PATH + "precipitation_all_comparison_bars.png", dpi=400, bbox_inches='tight')
plt.show()
print(f"Saved to {OUTPUT_PATH}precipitation_all_comparison_bars.png")
