# %% PLOT BAR CHART from extracted CSV files
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

os.chdir("/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp")

OUTPUT_PATH = "/leonardo/home/userexternal/mdasilva/leonardo_work/Otis_exp/figs/"

# Load all CSV files
print("Loading data from CSV files...")

# Experiments
experiments = ["noto", "wsm5", "wsm7", "wdm7", "ERA5", "CMORPH"]
data_dict = {}

for exp in experiments:
    file_path = f"extracted_timeseries/{exp}_acapulco_timeseries.csv"
    df = pd.read_csv(file_path)
    df['time'] = pd.to_datetime(df['time'])
    df = df.set_index('time')
    data_dict[exp] = df['precipitation_mm_per_hr']
    print(f"  Loaded {exp}: {len(df)} hours")

# Load station data
station = pd.read_csv("extracted_timeseries/station_acapulco_timeseries.csv")
station['time'] = pd.to_datetime(station['time'])
station = station.set_index('time')
station = station['pr']
print(f"  Loaded station: {len(station)} hours")

# Find common time steps (intersection of all)
print("\nFinding common time steps...")
common_times = station.index
for exp in experiments:
    common_times = common_times.intersection(data_dict[exp].index)

common_times = pd.DatetimeIndex(sorted(common_times))
print(f"Common time steps: {len(common_times)} hours")

# Create the plot
fig, ax = plt.subplots(figsize=(14, 6))
ax.set_facecolor('silver')

# Bar settings
n_datasets = len(experiments) + 1  # +1 for station
bar_width = 0.12
x_positions = np.arange(len(common_times))

# Colors
colors = {
    'noto': 'cyan',
    'wsm5': 'magenta', 
    'wsm7': 'brown',
    'wdm7': 'green',
    'ERA5': 'red',
    'CMORPH': 'blue',
    'station': 'black'
}

ax.axvspan(14, 19, color='white')

# Plot each experiment
for i, exp in enumerate(experiments):
    values = data_dict[exp].loc[common_times].values
    offset = (i - n_datasets/2) * bar_width
    ax.bar(x_positions + offset, values, width=bar_width,
           label=exp.upper(), color=colors[exp], alpha=0.75,
           edgecolor='black', linewidth=0.5)

# Plot station
font_size = 10
station_values = station.loc[common_times].values
offset_station = (len(experiments) - n_datasets/2) * bar_width
ax.bar(x_positions + offset_station, station_values, width=bar_width,
       label='Station', color=colors['station'], alpha=0.75,
       edgecolor='black', linewidth=0.5)

# Formatting
ax.set_xlabel('Time (24-25 Oct 2023)', fontsize=font_size, fontweight='bold')
ax.set_ylabel('Precipitation (mm/hr)', fontsize=font_size, fontweight='bold')
ax.set_title('(a)', loc='left', fontsize=font_size, fontweight='bold')

# Set x-ticks every 12 hours
tick_positions = x_positions
tick_labels = common_times.strftime('%m/%d %H:00')
ax.set_xticks(tick_positions)
ax.set_xticklabels(tick_labels, rotation=45, ha='right')

ax.legend(loc='upper right', fontsize=font_size, ncol=7)
ax.grid(True, linestyle='--', axis='y')
ax.set_ylim(bottom=0)

plt.tight_layout()

# Save figure
plt.savefig(OUTPUT_PATH + "pyplt_Hurricane_Otis_otis_prec.png", dpi=400, bbox_inches='tight')
plt.show()

print(f"\n✓ Plot saved to {OUTPUT_PATH}pyplt_Hurricane_Otis_otis_prec.png")
