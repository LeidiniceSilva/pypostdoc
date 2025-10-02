# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script .txt with cyclone genesis"

import argparse

parser = argparse.ArgumentParser(description="Script to save cyclone genesis")
parser.add_argument("--dataset", type=str, help="ERA5, CNRM-ESM2-1, EC-Earth3-Veg, GFDL-ESM4, HadGEM3-GC31-MM, MPI-ESM1-2-HR, MPI-ESM1-2-LR, NorESM-2MM, UKESM1-0-LL")
parser.add_argument("--domain", type=str, help="AFR, AUS, CAM, EAS, EUR, NAM, SAM, WAS")
args = parser.parse_args()

dataset = args.dataset
domain = args.domain

path = '/home/mda_silv/users/TRACK-CYCLONE/CORDEX-TF/S2R-Vortrack/Data/{0}/{1}/track'.format(dataset, domain)


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

	dt_list, lat_list, lon_list, vo_list, pcen_list = [], [], [], [], []
	
	for yr in range(2000, 2009+1):
	
		data = read_dat_file('{0}/resultado_{1}.dat'.format(path, yr))
	
		rows_list = []
		for i, (header, rows) in enumerate(data):
			rows_list.append(rows[0])
	
		for j  in rows_list:
			dt_list.append(int(j[0]))
			lat_list.append(float(j[1]))
			lon_list.append(float(j[2]))
			vo_list.append(float(j[3]))
			pcen_list.append(float(j[4]))

	
	return dt_list, lat_list, lon_list, vo_list, pcen_list

	
def write_list_to_dat(data):
	
	with open('{0}/genesis_density_{1}_{2}.txt'.format(path, dataset, domain), 'w') as f:
		f.write("data\tlat\tlon\tvo\tpres\n")
		for row in data:
			row_string = '\t'.join(map(str, row)) + '\n'
			f.write(row_string)

# Import model and obs dataset
list_dt, list_lat, list_lon, list_vo, list_pcen = open_dat_file()
data = [list_dt, list_lat, list_lon, list_vo, list_pcen]
data_ = list(map(list, zip(*data)))

# Write transposed data to .dat file
write_list_to_dat(data_)
exit()
