import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# -------- Simulated cyclone data --------
NUM_CYCLONES = 500
DAYS_PER_CYCLONE = 7

np.random.seed(42)

def simulate_cyclone_data(num_cyclones, days):
    mslp_all = []
    gust_all = []
    for _ in range(num_cyclones):
        mslp = np.random.normal(1005, 3, days) - np.array([3, 5, 6, 5, 3, 2, 1])[:days]
        gust = np.random.normal(15, 3, days) + np.random.choice([-1, 0, 1, 2, 3], days)
        mslp_all.append(mslp)
        gust_all.append(gust)
    return np.array(mslp_all), np.array(gust_all)

# Generate data
mslp_series, gust_series = simulate_cyclone_data(NUM_CYCLONES, DAYS_PER_CYCLONE)

# -------- Analysis --------
x_vals = []  # delta_day = day_mslp - day_gust
y_vals = []

for mslp, gust in zip(mslp_series, gust_series):
    day_mslp = np.argmin(mslp)
    day_gust = np.argmax(gust)
    delta_day = day_mslp - day_gust
    max_wind = np.max(gust)

    x_vals.append(delta_day)
    y_vals.append(max_wind)

x_vals = np.array(x_vals)
y_vals = np.array(y_vals)

# -------- Frequency Normalization --------
counter = Counter(x_vals)
max_count = max(counter.values())
norm_freq = np.array([counter[x] / max_count for x in x_vals])

# -------- Plotting --------
fig, ax = plt.subplots(figsize=(10, 6))

# Create colormap and plot scatter
cmap = plt.get_cmap('viridis')
colors = cmap(norm_freq)

sc = ax.scatter(x_vals, y_vals, c=colors, edgecolor='k', s=70, alpha=0.9)

# Add reference line and labels
ax.axvline(0, color='gray', linestyle='--', linewidth=1.2, label='Min MSLP = Max Wind')
ax.set_xlabel('Days (MSLP - Max Wind)', fontsize=12)
ax.set_ylabel('Max Wind Gust (m/s)', fontsize=12)
ax.set_title('Max Wind vs. Timing of Min MSLP\nColor by Normalized Frequency', fontsize=14)
ax.grid(True, linestyle='--', alpha=0.5)
ax.legend()

# Add colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(vmin=0, vmax=1))
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, label='Normalized Frequency')

plt.tight_layout()
plt.show()

