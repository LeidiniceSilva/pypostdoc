# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script write .txt with cyclone tracking"

idx = 1
dataset = 'ERA5'

if idx == 1:
	lat1 = str(-27.0)
	lat2 = str(-22.0)
	lon1 = str(-47.0)
	lon2 = str(-41.0)

	legend = 'region_i'
else:
	lat1 = str(-34.0)
	lat2 = str(-29.0)
	lon1 = str(-54.0)
	lon2 = str(-48.0)
	legend = 'region_ii'

path = '/marconi/home/userexternal/mdasilva/user/mdasilva/SAM-3km/post_cyclone/ECyclone'


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


def open_dat_file(dataset):

	dt_r, lat_r, lon_r, vo_r, pcen_r = [], [], [], [], []
	
	for yr in range(2018, 2021+1):
	
		data = read_dat_file('{0}/ECyclone_{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
	
		rows_list = []
		rows_list_i = []
		for i, (header, rows) in enumerate(data):
			if ((lat2 < rows[0][1] < lat1) and (lon2 < rows[0][2] < lon1)):
				rows_list.append(rows)
		
		print(rows_list)
		exit()
		
		for j  in rows_list:
			for k in j:
				dt_r.append(int(k[0]))
				lat_r.append(float(k[1]))
				lon_r.append(float(k[2]))
				vo_r.append(float(k[3]))
				pcen_r.append(float(k[4]))
	
	return dt_r, lat_r, lon_r, vo_r, pcen_r
    

def transpose_list(data):
    
	return list(map(list, zip(*data)))


def write_list_to_dat(data):
	
	with open('{0}/ECyclone_{1}/tracking_density_{1}_{2}.txt'.format(path, dataset, legend), 'w') as f:
		
		for row in data:
			# Convert elements to strings and join with tabs
			row_string = '\t'.join(map(str, row)) + '\n'  
			f.write(row_string)

# Import model and obs dataset
list_dt, list_lat, list_lon, list_vo, list_pcen = open_dat_file(dataset)
data = [list_dt, list_lat, list_lon, list_vo, list_pcen]

# Transpose the data to convert it into columns
transposed_data = transpose_list(data)

# Write transposed data to .dat file
write_list_to_dat(transposed_data)
exit()
			
	
