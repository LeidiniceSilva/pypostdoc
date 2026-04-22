# %% EXTRACT TIME SERIES - Process one by one and save to txt files
import os
import xarray as xr
import pandas as pd
import numpy as np

os.chdir("/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp")

# Helper function
def fix_longitude(ds, lon_var='lon'):
    if np.any(ds[lon_var] > 180):
        ds.coords[lon_var] = (ds.coords[lon_var] + 180) % 360 - 180
        ds = ds.sortby(ds[lon_var])
    return ds

# Location
acapulco_lat, acapulco_lon = 16.85, -99.91

# Create output directory for extracted data
os.makedirs("extracted_timeseries", exist_ok=True)

# Process each dataset individually
experiments = ["noto"] # "wsm5", "wsm7", "wdm7", "ERA5", "CMORPH"
data_dir = "exps_mp_regrided/"

for exp in experiments:
    print(f"\n{'='*50}")
    print(f"Processing {exp}...")
    print(f"{'='*50}")
    
    try:
        if exp == "ERA5":
            print("  Opening ERA5 file...")
            ds = xr.open_dataset("era5/1hr/tp_ERA5_Otis_1hr_Oct2023.nc")
            ds = ds.rename({"valid_time": "time", "latitude": "lat", "longitude": "lon", "tp": "pr"})
            ds = fix_longitude(ds, "lon")
            ds = ds.sortby(ds['lat'])
            ds["pr"] = ds["pr"] * 1000  # m -> mm
            
        elif exp == "CMORPH":
            print("  Opening CMORPH file...")
            ds = xr.open_dataset("cmorph/CMORPH_cut_domain_20251022_20251025.nc")
            ds = ds.rename({"cmorph": "pr"})
            ds = fix_longitude(ds, "lon")
            
        else:
            print(f"  Opening {exp} file...")
            ds = xr.open_dataset(data_dir + exp + f"/pr_{exp}_2023101900.nc")
            ds["pr"] = ds["pr"] * 3600  # kg m-2 s-1 -> mm/hr
        
        print("  Selecting time period...")
        ds = ds.sel(time=slice("2023-10-19T01:00:00", "2023-10-27T00:00:00"))
        
        print(f"  Extracting point at ({acapulco_lat}, {acapulco_lon})...")
        point = ds['pr'].sel(lat=acapulco_lat, lon=acapulco_lon, method='nearest')
        
        actual_lat = float(point.lat.values)
        actual_lon = float(point.lon.values)
        
        print(f"  Nearest grid point: ({actual_lat:.2f}, {actual_lon:.2f})")
        print(f"  Extracted {len(point.time.values)} hours")
        
        # Create DataFrame
        df = pd.DataFrame({
            'time': point.time.values,
            'precipitation_mm_per_hr': point.values
        })
        
        # Save to CSV
        output_file = f"extracted_timeseries/{exp}_acapulco_timeseries.csv"
        df.to_csv(output_file, index=False)
        print(f"  ✓ Saved to {output_file}")
        
        # Also save a simpler txt file
        txt_file = f"extracted_timeseries/{exp}_acapulco_timeseries.txt"
        with open(txt_file, 'w') as f:
            f.write(f"# Precipitation timeseries for {exp.upper()} at Acapulco\n")
            f.write(f"# Location: {actual_lat:.4f}°N, {actual_lon:.4f}°W\n")
            f.write(f"# Time period: 2023-10-19 01:00 to 2023-10-27 00:00\n")
            f.write(f"# Units: mm/hr\n")
            f.write(f"# Format: time, precipitation_mm_per_hr\n")
            for idx, row in df.iterrows():
                f.write(f"{row['time']}, {row['precipitation_mm_per_hr']:.4f}\n")
        
        print(f"  ✓ Saved to {txt_file}")
        
        # Clean up to free memory
        del ds
        del point
        del df
        
    except Exception as e:
        print(f"  ✗ Error processing {exp}: {str(e)}")
        continue

print(f"\n{'='*50}")
print("Processing station data...")
print(f"{'='*50}")

# Load and save station data separately
try:
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
    
    # Save station data
    station_df = station.reset_index()
    station_df.to_csv("extracted_timeseries/station_acapulco_timeseries.csv", index=False)
    
    with open("extracted_timeseries/station_acapulco_timeseries.txt", 'w') as f:
        f.write("# Station precipitation timeseries at Acapulco\n")
        f.write("# Source: SEMAR station\n")
        f.write("# Units: mm/hr\n")
        f.write("# Format: time, precipitation_mm_per_hr\n")
        for idx, row in station_df.iterrows():
            f.write(f"{row['time']}, {row['pr']:.4f}\n")
    
    print(f"  ✓ Station data saved to extracted_timeseries/station_acapulco_timeseries.csv")
    print(f"  ✓ Station data saved to extracted_timeseries/station_acapulco_timeseries.txt")
    
except Exception as e:
    print(f"  ✗ Error processing station: {str(e)}")

print(f"\n{'='*50}")
print("✓ EXTRACTION COMPLETE!")
print(f"{'='*50}")
print("Files saved in: extracted_timeseries/")
print("\nFiles created:")
for exp in experiments:
    print(f"  - {exp}_acapulco_timeseries.csv")
    print(f"  - {exp}_acapulco_timeseries.txt")
print("  - station_acapulco_timeseries.csv")
print("  - station_acapulco_timeseries.txt")


