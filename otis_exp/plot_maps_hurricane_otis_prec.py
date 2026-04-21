# %% ONE FIGURE - All hourly precipitation comparisons
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

# %% PLOT - ONE FIGURE with all datasets
fig, ax = plt.subplots(figsize=(14, 6))

colors = {
    'noto': '#33a02c',
    'wsm5': '#005B72', 
    'wsm7': '#905dc7',
    'wdm7': '#692510',
    'ERA5': '#fdbf6f',
    'CMORPH': '#fb9a99'
}

# Plot each model
for exp in experiments:
    # Align to common times
    mask = np.isin(times_series[exp][0], common_times)
    times = times_series[exp][0][mask]
    values = times_series[exp][1][mask]
    ax.plot(times, values, label=exp.upper(), color=colors[exp], linewidth=1.5, alpha=0.7)

# Plot station
station_aligned = station.loc[common_times]
ax.plot(station_aligned.index, station_aligned['pr'], 
        label='Station', color='black', linewidth=2, marker='o', markersize=2, alpha=0.8)

# Formatting
ax.set_xlabel('Time', fontsize=12)
ax.set_ylabel('Precipitation [mm/hr]', fontsize=12)
ax.set_title(f'Hourly Precipitation at Acapulco ({acapulco_lat:.2f}°, {acapulco_lon:.2f}°)\nHurricane Otis (Oct 19-27, 2023)', fontsize=14)
ax.legend(loc='upper right', fontsize=10, ncol=2)
ax.grid(True, alpha=0.3, linestyle='--')
ax.set_ylim(bottom=0)
plt.xticks(rotation=45)
plt.tight_layout()

# Save
plt.savefig(OUTPUT_PATH + "pyplt_Hurricane_Otis_otis_prec.png", dpi=400, bbox_inches='tight')
plt.show()
print(f"Saved to {OUTPUT_PATH}precipitation_all_comparison.png")
