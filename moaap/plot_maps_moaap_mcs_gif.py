# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 03, 2026"
__description__ = "This script plot MCSs"

import os
import glob
import pickle
import netCDF4
import cartopy
import warnings
import numpy as np
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from tqdm import tqdm
from PIL import Image
from IPython import display
from IPython.display import Image

warnings.filterwarnings("ignore")


def make_gif(frame_folder):
    frames = [Image.open(image) for image in np.sort(glob.glob(f"{frame_folder}/*.jpg"))]
    frame_one = frames[0]
    frame_one.save("phenomenon_mcs.gif", format="GIF", append_images=frames,
               save_all=True, duration=100, loop=0)


# Import data
data_vars = xr.open_dataset('/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/ERA5/CAR-4/input/CAR-4_ERA5_reanalysis_1hr_2000010100.nc')
time_datetime = pd.to_datetime(np.array(data_vars['time'].values, dtype='datetime64'))

# load MCS characteristics
data_moaap = xr.open_dataset('/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/ERA5/CAR-4/output/20000101_ERA5_evaluation_ObjectMasks__dt-1h_MOAAP-masks.nc')

object_names = [['MCS', '#33a02c', '-', 2]]

for tt in tqdm(range(len(time_datetime))):

    # Create a figure and axis with a PlateCarree projection (latitude and longitude)
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()},
                          figsize=(14,6))

    # Add coastlines and gridlines for reference
    ax.coastlines(color='#969696')
    ax.gridlines()

    # MCSs
    sc = plt.contour(data_moaap['lon'],
                        data_moaap['lat'],
                        np.array(data_moaap['MCS_Tb_Objects'][tt,:,:]),
                        colors = '#33a02c')   

    ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

    # Set the title of the plot
    plt.title('Objects identied by MOAAP at '+str(time_datetime[tt])[:16])

    # create legend
    for ob in range(len(object_names)):
      plt.plot([],[], color = object_names[ob][1], \
               linestyle = object_names[ob][2],\
               lw = object_names[ob][3],\
               label = object_names[ob][0])

    # plt.legend()
    ax.legend(bbox_to_anchor=(1, 0.00), ncol=4)

    # Show the plot
    fig.savefig('/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/figs/'+str(tt).zfill(3)+'_CausesOfExtreme_100_daily-extremes.jpg', bbox_inches='tight', dpi=100)
    plt.close()

# Make gif
make_gif("/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/figs/")

with open("phenomenon_mcs.gif",'rb') as f:
    display.Image(data=f.read(), format='jpg')

fname = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/figs/phenomenon_mcs.gif'
Image(open(fname, 'rb').read())

