# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot composites"

import netCDF4
import numpy as np
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature

from datetime import datetime, timedelta

path='/marconi/home/userexternal/mdasilva'


def generate_hourly_dates(start_date, end_date):
    
	dates_hr = []
	current_date = start_date
	while current_date <= end_date:
		dates_hr.append(current_date.strftime('%Y%m%d%H'))
		current_date += timedelta(hours=24)
	
	return dates_hr
	

def find_indices_in_date_list(date_list, target_dates):
    
	indices = []
	for target_date in target_dates:
		try:
			index = date_list.index(target_date)
			indices.append(index)
		except ValueError:
			pass  # Date not found in date_list
	
	return indices
	
	
def read_dat_file(filename):

	data = []
	with open(filename, 'r') as file:
		lines = file.readlines()
		header = []
		rows = []
	
	# Iterate over lines in the file
	for line in lines:
		line = line.strip().split()
		
		# If the line contains 6 elements, it's considered a header
		if len(line) == 6:
			
			# If we have rows, append them to data
			if rows:
				data.append((header, rows))
				rows = []
			header = line
		else:  
			rows.append(line)
	
	# Append the last header and rows to data
	if header and rows:
		data.append((header, rows))
	
	return data


def open_file_dt(dataset):

	dt_hr = []
	for yr in range(2018, 2021+1):
	
		data = read_dat_file('{0}/user/mdasilva/SAM-3km/post_cyclone/ECv2/{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
	
		rows_list = []
		for i, (header, rows) in enumerate(data):
			if (rows[0][2] < str(-56)):
				rows_list.append(rows)
						
		for j  in rows_list:
			for k in j:
				dt_hr.append(str(k[0][:]))

	return dt_hr



def open_file_coords(dataset):

	latitudes, longitudes = [], []
	
	# Open tracking file 
	for yr in range(2018, 2018+1):
		data = read_dat_file('{0}/user/mdasilva/SAM-3km/post_cyclone/ECv2/{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
	
		# Select the first row
		rows_list = []
		for i, (header, rows) in enumerate(data):
			if (rows[0][2] < str(-56)):
				rows_list.append(rows[0])
						
		# Select the coords
		for j  in rows_list:
			latitudes.append(float(j[1]))
			longitudes.append(float(j[2]))

	coordinates = [(lat, lon) for lat, lon in zip(latitudes, longitudes)]
	
	return coordinates
	
	
def import_data(param, dataset, indices):

	if dataset == 'RegCM5':
		arq = '{0}/user/mdasilva/SAM-3km/post_cyclone/regcm5/regcm5/{1}_SAM-3km_{2}_6hr_2018-2021_lonlat.nc'.format(path, param, dataset)
	elif dataset == 'WRF415':
		arq = '{0}/user/mdasilva/SAM-3km/post_cyclone/wrf/wrf/{1}_SAM-3km_{2}_6hr_2018-2021_lonlat.nc'.format(path, param, dataset)
	else:
		arq   = '{0}/user/mdasilva/SAM-3km/post_cyclone/obs/era5/era5/{1}_SAM-25km_{2}_6hr_2018-2021_lonlat.nc'.format(path, param, dataset)		
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]	
	lon   = data.variables['lon'][:]
	mean = var[:][0,:,:]
		
	return lat, lon, mean
	
	
# Generate list of dates from 2018 to 2021
hourly_dates = generate_hourly_dates(datetime(2018, 1, 1, 0), datetime(2021, 12, 31, 23))

# Import cyclone tracking date 
dt_era5 = open_file_dt('ERA5')
dt_regcm5 = open_file_dt('RegCM5')
dt_wrf415 = open_file_dt('WRF415')

# Import cyclone tracking coords 
coords_era5 = open_file_coords('ERA5')
coords_regcm5 = open_file_coords('RegCM5')
coords_wrf415 = open_file_coords('WRF415')

era5_idx = find_indices_in_date_list(hourly_dates, dt_era5)
regcm5_idx = find_indices_in_date_list(hourly_dates, dt_regcm5)
wrf415_idx = find_indices_in_date_list(hourly_dates, dt_wrf415)

# Import model and obs dataset 
lat, lon, u10_era5 = import_data('u10', 'ERA5', era5_idx_i)
lat, lon, v10_era5 = import_data('v10', 'ERA5', era5_idx_i)

ws10_era5 = np.sqrt(u10_era5**2 + v10_era5**2)

# Calculate composites
radius = 10  # Radius in degrees

for center_era5 in coords_era5:
    center_lat_era5, center_lon_era5 = center_era5
    
    # Create a mask for points within the radius
    distance = np.sqrt((lat_grid - center_lat)**2 + (lon_grid - center_lon)**2)
    mask = distance <= radius
    
    # Simulate some wind data within the radius for demonstration
    wind_data = np.random.rand(*lon_grid.shape) * mask  # Replace with actual data extraction

    # Accumulate wind data in the composite and count occurrences
    composite_wind += wind_data * mask
    count += mask

# Finalize the composite by averaging
composite_wind /= count
composite_wind[count == 0] = np.nan  # Set no data points to NaN for plotting

# Plotting the composite
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

# Add geographical features
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)

# Create filled contours for the composite wind
contour = ax.contourf(lon_grid, lat_grid, composite_wind, cmap='viridis', 
                       levels=np.linspace(0, 1, 21), transform=ccrs.PlateCarree())
plt.colorbar(contour, ax=ax, label='Composite Wind Speed (units)')

# Set title and labels
plt.title('Centralized Cyclone Composite with 20° Radius')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
plt.show()


exit
var = 'msl'
dataset = 'ERA5'
freq = 'day'
path='/marconi/home/userexternal/mdasilva'


def read_dat_file(filename):

	data = []
	with open(filename, 'r') as file:
		lines = file.readlines()
		header = []
		rows = []
	
	# Iterate over lines in the file
	for line in lines:
		line = line.strip().split()
		
		# If the line contains 6 elements, it's considered a header
		if len(line) == 6:
			
			# If we have rows, append them to data
			if rows:
				data.append((header, rows))
				rows = []
			header = line
		else:  
			rows.append(line)
	
	# Append the last header and rows to data
	if header and rows:
		data.append((header, rows))
	
	return data


def open_dat_file():

	latitudes, longitudes = [], []
	
	# Open tracking file 
	for yr in range(2018, 2018+1):
		data = read_dat_file('{0}/user/mdasilva/SAM-3km/post_cyclone/ECyclone/ECyclone_{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
	
		# Select the first row
		rows_list = []
		for i, (header, rows) in enumerate(data):
			rows_list.append(rows[0])
		
		# Select the coords
		for j  in rows_list:
			latitudes.append(float(j[1]))
			longitudes.append(float(j[2]))

	coordinates = [(lat, lon) for lat, lon in zip(latitudes, longitudes)]
	
	return coordinates


def extract_cyclone_data(dataset_cyclone, cyclone_center, radius):

    lat_center, lon_center = cyclone_center
    lat_slice = slice(lat_center - radius*0.25, lat_center + radius*0.25)
    lon_slice = slice(lon_center - radius*0.25, lon_center + radius*0.25)
    subset = dataset_cyclone.sel(latitude=lat_slice, longitude=lon_slice)

    return subset


# Import dataset
dataset_i = xr.open_dataset('{0}/user/mdasilva/SAM-3km/post_cyclone/obs/era5/era5/'.format(path) + '{0}_SAM-25km_{1}_{2}_2018-2021.nc'.format(var, dataset, freq), engine='netcdf4', decode_times=False)

# Import cyclone tracking coords 
coords_i = open_dat_file()

# Define the radius
radius = 79

# Extract data around cyclone
cyclone_data = []
for center in coords_i:
	data = extract_cyclone_data(dataset_i, center, radius)
	cyclone_data.append(data)

# Compute the composite
composite = sum(cyclone_data) / len(cyclone_data)

# Save file
composite.to_netcdf('{0}/user/mdasilva/SAM-3km/post_cyclone/ECyclone/ECyclone_{1}/composite_cyclones_{1}_{2}_2018-2021.nc'.format(path, dataset, var))

