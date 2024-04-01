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

# create date list
dt = pd.date_range('2018-01-01','2021-12-31', freq='D')

# Reading era5 tracking
for yr in range(2018, 2021+1):
	df = pd.read_csv(os.path.join('{0}/ECyclone_ERA5/track/'.format(path), 'resultado_{0}.dat'.format(yr)))
print(df)
exit()
	
# Reading RegCM5 tracking
for yr in range(2018, 2021+1):
	df = pd.read_csv(os.path.join('{0}/ECyclone_RegCM5/track/'.format(path), 'resultado_{0}.dat'.format(yr)))
	print(df)
	exit()
			
