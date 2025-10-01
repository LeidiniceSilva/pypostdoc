# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Sep 22, 2025"
__description__ = "This script plots composite around cyclone center"

import os
import netCDF4
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeat

from matplotlib.patches import Circle
from scipy.interpolate import griddata
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

path = '/afs/ictp.it/home/m/mda_silv/Downloads'
font_size = 10
radius = 4.0  # graus
res_r = 0.05  # radial resolution
res_theta = 1 # angular resolution in degrees


def import_data(param):

	arq = '{0}/{1}_Mar26.nc'.format(path, param)
	data = netCDF4.Dataset(arq)
	var = data.variables[param][:]    
	lat = data.variables['lat'][:]
	lon = data.variables['lon'][:]

	if param == 'pr':
		var = var * 3600  # Convert to mm/h (assuming kg/m²/s input)
	if param == 'psl':
		var = var / 100.0  # Convert Pa to hPa
	data.close()
    
	return lat, lon, var


def create_polar_composite(data_var, lat2d, lon2d, psl, radius=4.0, res_r=0.05, res_theta=1):
    
	# Create polar grid
	r = np.arange(0, radius + res_r, res_r)
	theta = np.deg2rad(np.arange(0, 360, res_theta))
	R, TH = np.meshgrid(r, theta)
    
	# Convert polar to Cartesian coordinates
	X = R * np.cos(TH)
	Y = R * np.sin(TH)
    
	composite_list = []
	valid_centers = 0
    
	for t in range(psl.shape[0]):
		# Find cyclone center (MSLP minimum)
		iy, ix = np.unravel_index(np.nanargmin(psl[t, :, :]), psl[t, :, :].shape)
		lat_c = lat2d[iy, ix]
		lon_c = lon2d[iy, ix]
        
		# Skip if center is at edge
		if (iy == 0 or iy == lat2d.shape[0]-1 or ix == 0 or ix == lat2d.shape[1]-1):
			continue

		# Relative coordinates with latitude correction
		x_rel = (lon2d - lon_c) * np.cos(np.deg2rad(lat_c))
		y_rel = lat2d - lat_c
        
		# Create mask for points within radius
		mask = x_rel**2 + y_rel**2 <= radius**2
		points = np.column_stack((x_rel[mask], y_rel[mask]))
		values = data_var[t, :, :][mask]
        
		# Interpolate to polar grid
		if len(points) > 10:  # Minimum points for interpolation
			try:
				grid_z = griddata(points, values, (X, Y), method='linear')
				if np.any(~np.isnan(grid_z)):
					composite_list.append(grid_z)
					valid_centers += 1
			except:
				continue
    
	print(f"Valid cyclone centers found: {valid_centers}")
    
	if len(composite_list) == 0:
		print("Warning: No valid composites created")
		return None, R, TH
    
	# Create final composite
	composite = np.nanmean(np.array(composite_list), axis=0)
    
	return composite, R, TH


def create_wind_composite(u_var, v_var, lat2d, lon2d, psl, radius=4.0, res_r=0.2, res_theta=20):
    
    # Create coarser polar grid for wind arrows
    r = np.arange(0, radius + res_r, res_r)
    theta = np.deg2rad(np.arange(0, 360, res_theta))
    R_wind, TH_wind = np.meshgrid(r, theta)
    
    # Convert polar to Cartesian coordinates
    X_wind = R_wind * np.cos(TH_wind)
    Y_wind = R_wind * np.sin(TH_wind)
    
    u_composite_list = []
    v_composite_list = []
    valid_centers = 0
    
    for t in range(psl.shape[0]):
        # Find cyclone center (MSLP minimum)
        iy, ix = np.unravel_index(np.nanargmin(psl[t, :, :]), psl[t, :, :].shape)
        lat_c = lat2d[iy, ix]
        lon_c = lon2d[iy, ix]
        
        # Skip if center is at edge
        if (iy == 0 or iy == lat2d.shape[0]-1 or ix == 0 or ix == lat2d.shape[1]-1):
            continue

        # Relative coordinates with latitude correction
        x_rel = (lon2d - lon_c) * np.cos(np.deg2rad(lat_c))
        y_rel = lat2d - lat_c
        
        # Create mask for points within radius
        mask = x_rel**2 + y_rel**2 <= radius**2
        points = np.column_stack((x_rel[mask], y_rel[mask]))
        u_values = u_var[t, :, :][mask]
        v_values = v_var[t, :, :][mask]
        
        # Interpolate to polar grid
        if len(points) > 10:
            try:
                u_grid = griddata(points, u_values, (X_wind, Y_wind), method='linear')
                v_grid = griddata(points, v_values, (X_wind, Y_wind), method='linear')
                
                if np.any(~np.isnan(u_grid)) and np.any(~np.isnan(v_grid)):
                    u_composite_list.append(u_grid)
                    v_composite_list.append(v_grid)
                    valid_centers += 1
            except:
                continue
    
    print(f"Valid wind centers found: {valid_centers}")
    
    if len(u_composite_list) == 0:
        print("Warning: No valid wind composites created")
        return None, None, R_wind, TH_wind
    
    # Create final composites
    u_composite = np.nanmean(np.array(u_composite_list), axis=0)
    v_composite = np.nanmean(np.array(v_composite_list), axis=0)
    
    return u_composite, v_composite, R_wind, TH_wind


def calculate_precipitation_sum(pr_composite, psl_composite, R, TH):
    """Calculate precipitation sum patterns relative to MSLP composite"""
    
    if pr_composite is None or psl_composite is None:
        return None, None
    
    # Use absolute precipitation values (no normalization)
    pr_absolute = pr_composite
    
    # Calculate precipitation sum in different quadrants
    quadrants = {
        'NE': (0, np.pi/2),
        'SE': (np.pi/2, np.pi),
        'SW': (np.pi, 3*np.pi/2),
        'NW': (3*np.pi/2, 2*np.pi)
    }
    
    quadrant_sums = {}
    for quadrant, (theta_min, theta_max) in quadrants.items():
        mask = (TH >= theta_min) & (TH < theta_max)
        if np.any(mask):
            # Sum of precipitation in the quadrant
            quadrant_sums[quadrant] = np.nansum(pr_composite[mask])
        else:
            quadrant_sums[quadrant] = 0
    
    return pr_absolute, quadrant_sums


def calculate_total_precipitation(pr_composite, R, TH):
    """Calculate total precipitation statistics"""
    
    if pr_composite is None:
        return 0, 0, 0
    
    # Total precipitation sum over entire area
    total_sum = np.nansum(pr_composite)
    
    # Maximum precipitation value
    max_precip = np.nanmax(pr_composite)
    
    # Area with precipitation > 1 mm
    high_precip_area = np.sum(pr_composite > 1.0)
    
    return total_sum, max_precip, high_precip_area


# Import data
print("Importing data...")
lat, lon, psl = import_data('psl')
lat, lon, pr = import_data('pr')  # Precipitation
lat, lon, ta925 = import_data('ta925')
lat, lon, ua925 = import_data('ua925')
lat, lon, va925 = import_data('va925')
       
# Create 2D grid if needed
if lat.ndim == 1 and lon.ndim == 1:
	lon2d, lat2d = np.meshgrid(lon, lat)
else:
	lon2d, lat2d = lon, lat
        
# Create composites
print("Creating composites...")
psl_composite, R, TH = create_polar_composite(psl, lat2d, lon2d, psl, radius, res_r, res_theta)
pr_composite, R, TH = create_polar_composite(pr, lat2d, lon2d, psl, radius, res_r, res_theta)
ta_composite, R, TH = create_polar_composite(ta925, lat2d, lon2d, psl, radius, res_r, res_theta)
u_wind_composite, v_wind_composite, R_wind, TH_wind = create_wind_composite(ua925, va925, lat2d, lon2d, psl)

# Calculate precipitation sums (not normalized)
pr_absolute, quadrant_sums = calculate_precipitation_sum(pr_composite, psl_composite, R, TH)
total_precip, max_precip, high_precip_area = calculate_total_precipitation(pr_composite, R, TH)

# Calculate wind speed for statistics
if u_wind_composite is not None and v_wind_composite is not None:
    wind_speed = np.sqrt(u_wind_composite**2 + v_wind_composite**2)

# Print composite statistics
print('COMPOSITE STATISTICS:')
print(f'MSLP range: {np.nanmin(psl_composite):.1f} - {np.nanmax(psl_composite):.1f} hPa')
print(f'Precipitation sum range: {np.nanmin(pr_composite):.2f} - {np.nanmax(pr_composite):.2f} mm/h')
print(f'Temperature range: {np.nanmin(ta_composite):.1f} - {np.nanmax(ta_composite):.1f} K')
print(f'Wind speed range: {np.nanmin(wind_speed):.1f} - {np.nanmax(wind_speed):.1f} m/s')

# Create plot with precipitation sum (not normalized)
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(12, 10))

# Precipitation sum (absolute values, not normalized)
if np.any(~np.isnan(pr_composite)):
    # Use linear levels for absolute precipitation values
    pr_max = np.nanmax(pr_composite)
    if pr_max > 0:
        # Create levels based on actual precipitation values
        if pr_max <= 10:
            pr_levels = np.arange(0, 10.5, 0.5)  # Small range
        elif pr_max <= 30:
            pr_levels = np.arange(0, 20, 1)      # Medium range  
        else:
            pr_levels = np.arange(0, pr_max + 5, 5)  # Large range
        
        cf = ax.contourf(TH, R, pr_composite, levels=pr_levels, cmap='YlGnBu', extend='max')
        cbar = plt.colorbar(cf, ax=ax, shrink=0.8, pad=0.02)
        cbar.set_label('Precipitation (mm)', fontsize=font_size, weight='bold')
        
        # Add precipitation contour lines for key values
        key_levels = [1, 5, 10, 20, 30]  # Key precipitation thresholds
        valid_levels = [level for level in key_levels if level <= pr_max]
        if valid_levels:
            ct_pr = ax.contour(TH, R, pr_composite, levels=valid_levels, colors='blue', linewidths=1.0, linestyles='-', alpha=0.8)
            ax.clabel(ct_pr, inline=True, fontsize=font_size-2, fmt='%.0f mm')

# Temperature contours  
if np.any(~np.isnan(ta_composite)):
    temp_levels = np.linspace(np.nanmin(ta_composite), np.nanmax(ta_composite), 10)
    ct1 = ax.contour(TH, R, ta_composite, levels=temp_levels, colors='red', linewidths=1.2, linestyles='--')
    ax.clabel(ct1, inline=True, fontsize=font_size-2, fmt='%.0f K')
            
# Pressure contours 
if np.any(~np.isnan(psl_composite)):
    psl_min = np.nanmin(psl_composite)
    psl_max = np.nanmax(psl_composite)
    psl_step = max(1, round((psl_max - psl_min) / 6))
    psl_levels = np.arange(round(psl_min), round(psl_max) + psl_step, psl_step)
    ct2 = ax.contour(TH, R, psl_composite, levels=psl_levels, colors='black', linewidths=2.0, linestyles='-')
    ax.clabel(ct2, inline=True, fontsize=font_size-2, fmt='%.0f hPa')

# Wind vectors 
if u_wind_composite is not None and v_wind_composite is not None:
    u_polar = u_wind_composite * np.cos(TH_wind) - v_wind_composite * np.sin(TH_wind)
    v_polar = u_wind_composite * np.sin(TH_wind) + v_wind_composite * np.cos(TH_wind)
    scale = 25 / np.nanmax(wind_speed) if np.nanmax(wind_speed) > 0 else 25
    q = ax.quiver(TH_wind, R_wind, u_polar, v_polar, color='darkgreen', width=0.003, headwidth=4, headlength=5, alpha=0.9)
    ax.quiverkey(q, 0.92, 0.98, 5, '5 m/s', labelpos='E', fontproperties={'size': font_size-1, 'weight': 'bold'}, color='darkgreen')

# Configure polar plot
ax.set_theta_zero_location('N')  # 0° at top (North)
ax.set_theta_direction(-1)       # clockwise direction
ax.set_rlim(0, radius)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.grid(True, alpha=0.5)

# Add distance circles and labels
for dist in [1, 2, 3, 4]:
    ax.plot(np.linspace(0, 2*np.pi, 100), [dist]*100, 'k--', alpha=0.4, linewidth=0.8)

# Save figura
path_out = '{0}'.format(path)
name_out = 'pyplt_maps_composite_Hurricane_Catarina_26Mar2004.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
