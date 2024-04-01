# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot annual cycle of EC over SAM"

import os
import numpy as np
import pandas as pd

from netCDF4 import Dataset

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
                if rows:  # If we have rows, append them to data
                    data.append((header, rows))
                    rows = []  # Reset rows
                header = line
            else:  # Otherwise, it's a data row
                rows.append(line)

        # Append the last header and rows to data
        if header and rows:
            data.append((header, rows))

    return data


yr_1 = 2018
filename_1 = '{0}/ECyclone_ERA5/track/resultado_{1}.dat'.format(path, yr_1)
data_1 = read_dat_file(filename_1)
header_list_1 = []
for i, (header, rows) in enumerate(data_1):
	header_list_1.append(header[1])
	
yr_2 = 2019
filename_2 = '{0}/ECyclone_ERA5/track/resultado_{1}.dat'.format(path, yr_2)
data_2 = read_dat_file(filename_2)
header_list_2 = []
for i, (header, rows) in enumerate(data_2):
	header_list_2.append(header[1])
	
yr_3 = 2020
filename_3 = '{0}/ECyclone_ERA5/track/resultado_{1}.dat'.format(path, yr_3)
data_3 = read_dat_file(filename_3)
header_list_3 = []
for i, (header, rows) in enumerate(data_3):
	header_list_3.append(header[1])
	
yr_4 = 2021
filename_4 = '{0}/ECyclone_ERA5/track/resultado_{1}.dat'.format(path, yr_4)
data_4 = read_dat_file(filename_4)
header_list_4 = []
for i, (header, rows) in enumerate(data_4):
	header_list_4.append(header[1])

header_list = header_list_1 + header_list_2 + header_list_3 + header_list_4

print(len(header_list))
exit()

	
