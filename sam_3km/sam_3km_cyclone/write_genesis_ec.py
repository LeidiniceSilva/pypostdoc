# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Apr 01, 2024"
__description__ = "This script plot annual cycle of EC over SAM"

import os
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from netCDF4 import Dataset
from datetime import datetime

dataset = 'RegCM5'
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


def open_dat_file(yr):

	data = read_dat_file('{0}/ECyclone_{1}/track/resultado_{2}.dat'.format(path, dataset, yr))
	
	rows_list = []
	for i, (header, rows) in enumerate(data):
		rows_list.append(rows[0])
	
	return rows_list
    

# Import model and obs dataset
ec_2018 = open_dat_file('2018')
ec_2019 = open_dat_file('2019')
ec_2020 = open_dat_file('2020')
ec_2021 = open_dat_file('2021')

genesis_list = ec_2018 + ec_2019 + ec_2020 + ec_2021

with open('{0}/ECyclone_{1}/genesis_density_{1}.txt'.format(path, dataset), 'w') as f:
	f.write(str(genesis_list))
	f.close()	
