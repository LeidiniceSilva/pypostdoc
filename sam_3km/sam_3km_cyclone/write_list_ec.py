# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot map of precipitation"

import os
import netCDF4
import datetime
import numpy as np

from datetime import datetime, timedelta

dataset = 'ERA5'
iyr, fyr = 2018, 2021
path='/marconi/home/userexternal/mdasilva'


def add_second_list(data_list, second_date_list):

	# The length of date list and second date list must be the same
	if len(data_list) != len(second_date_list):
		raise ValueError("The length must be the same")
    
	for i in range(len(data_list)):
		data_list[i].append(second_date_list[i])
	
	return data_list
    
    
def gather_and_repeat(date_string, suffix, n):
        
	concatenated_string = str(date_string) + str(suffix)
	
	return [concatenated_string] * int(n)

    
def read_dat_file(filename):

	data = []
	with open(filename, 'r') as file:
		lines = file.readlines()
		header = []
		rows = []
	
	# Iterate over lines in the file
	for line in lines:
		line = line.strip().split()
		print(line)
		exit()
		
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


def save_dat_file(data_list, file_name):

    with open(file_name, 'w') as file:
        
	# Write headers
        file.write('\t'.join(headers) + '\n')
        
        # Write data
        for sublist in data_list:
            for data in sublist:
                line = '\t'.join(data)
                file.write(line + '\n')
		
		
# Open .dat file
for yr in range(iyr, fyr+1):
	print(yr)

	data = read_dat_file('{0}/user/mdasilva/SAM-3km/post_cyclone/ECyclone/ECyclone_{1}/track/resultado_{2}.dat'.format(path, dataset, yr))

	idx_list = []
	for i, (header, rows) in enumerate(data):
		index = header[0]
		datelist = header[1]
		lifetime = header[4]
			
		repeated_dt = gather_and_repeat(datelist, index, lifetime)
					
		idx_list.append(add_second_list(rows, repeated_dt))
		
	# Save .dat file 
	headers = ['date', 'lat', 'lon', 'vo', 'pres', 'idx']
	save_dat_file(idx_list, 'IDX_list_cyclone_{0}_{1}.dat'.format(dataset, yr))
	
exit()

