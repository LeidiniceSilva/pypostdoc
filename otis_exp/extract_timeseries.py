# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 17, 2026"
__description__ = "This script extract timeseries"

import os
import xarray as xr
import pandas as pd
import numpy as np

os.chdir("/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp")

# Load and save station 
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
    station_df.to_csv("ws/station_acapulco_timeseries.csv", index=False)
    
    with open("ws/station_acapulco_timeseries.txt", 'w') as f:
        f.write("# Station precipitation timeseries at Acapulco\n")
        f.write("# Source: SEMAR station\n")
        f.write("# Units: mm/hr\n")
        f.write("# Format: time, precipitation_mm_per_hr\n")
        for idx, row in station_df.iterrows():
            f.write(f"{row['time']}, {row['pr']:.4f}\n")
    
    print(f"Station data saved to ws/station_acapulco_timeseries.csv")
    print(f"Station data saved to ws/station_acapulco_timeseries.txt")
    
except Exception as e:
    print(f"  ✗ Error processing station: {str(e)}")

