import os
import numpy as np
import xarray as xr
import scipy as sc
import scipy.ndimage
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeat

from scipy import signal, misc
from matplotlib.patches import Rectangle
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


# Create a figure
fig = plt.figure(figsize=(8, 6))

# Define the grid layout: 2 rows, 1 column for the main subplots
gs = gridspec.GridSpec(2, 1, figure=fig)

# Create the first two subplots (2 rows, 1 column)
ax1 = fig.add_subplot(gs[0, 0])  # First subplot
ax2 = fig.add_subplot(gs[1, 0])  # Second subplot

# Add a third subplot centered on the second line
# Use gridspec to create a new subplot in the center of the second row
gs_center = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1, 0])
ax3 = fig.add_subplot(gs_center[0, 1])  # Centered subplot

# Plot something in the first two subplots
ax1.plot([0, 1], [0, 1], label="Line 1", color="blue")
ax1.set_title("Subplot 1")
ax2.plot([0, 1], [1, 0], label="Line 2", color="red")
ax2.set_title("Subplot 2")

# Plot something in the third subplot (centered)
ax3.plot([0, 1], [0.5, 0.5], label="Line 3", color="green")
ax3.set_title("Subplot 3 (Centered)")

# Adjust layout for better spacing
plt.tight_layout()

# Show the plot
plt.show()
