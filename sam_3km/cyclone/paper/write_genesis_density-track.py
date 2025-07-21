# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script .txt with cyclone genesis"

dataset = 'ERA5'
path = '/leonardo/home/userexternal/mdasilva/leonardo_work/SAM-3km'


def read_dat_file(filename):

	with open(filename, 'r') as file:
		lines = file.readlines()
	
	return lines


def open_dat_file(dataset):

	dt_list, lat_list, lon_list, vo_list, pcen_list = [], [], [], [], []
	
	for yr in range(2018, 2021+1):
		data = read_dat_file('{0}/postproc/cyclone/{1}/track/resul_era10{2}.dat'.format(path, dataset, yr))

		for line in data:
			parts = line.strip().split()
			if len(parts) == 5:
				dt_list.append(int(parts[0]))
				lat_list.append(float(parts[1]))
				lon_list.append(float(parts[2]))
				vo_list.append(float(parts[3]))
				pcen_list.append(float(parts[4]))

	return dt_list, lat_list, lon_list, vo_list, pcen_list

	
def write_list_to_dat(data):
	
	with open('{0}/postproc/cyclone/{1}/genesis_density-track_{1}.txt'.format(path, dataset), 'w') as f:
		f.write("data\tlat\tlon\tvo\tpres\n")
		for row in data:
			row_string = '\t'.join(map(str, row)) + '\n'
			f.write(row_string)

# Import model and obs dataset
list_dt, list_lat, list_lon, list_vo, list_pcen = open_dat_file(dataset)
data = [list_dt, list_lat, list_lon, list_vo, list_pcen]
data_ = list(map(list, zip(*data)))

# Write transposed data to .dat file
write_list_to_dat(data_)
exit()
