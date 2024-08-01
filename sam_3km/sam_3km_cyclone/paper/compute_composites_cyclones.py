# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot composites"

import numpy as np
import xarray as xr
import pandas as pd

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

