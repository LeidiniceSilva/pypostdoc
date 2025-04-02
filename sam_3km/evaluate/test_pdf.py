#!/usr/bin/env python3

import glob
import xarray as xr
import numpy as np
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt

rain_minimum = 1.0
rain_maximum = 500.0
rain_step = 1.0
spd = 86400.0
sph = 3600.0

nrainbin = int((rain_maximum-rain_minimum)/rain_step) + 1

rain_pdf = np.linspace(rain_minimum,rain_maximum, num=nrainbin, endpoint=True)
rain_bin = np.zeros((nrainbin-1))

for ncfile in glob.glob('/leonardo/home/userexternal/mdasilva/leonardo_work/SAM-3km/output/SAM-3km_STS.2020080100.nc'):
    print(ncfile)
    print()
    pr = xr.load_dataset(ncfile)
    size = pr.pr.shape
    print(size)
    print()

    for time in range(size[0]):
        data = pr.pr[time,Ellipsis] * spd
        print(data)
        print()
        nzd = data.values[tuple(np.where(data >= rain_minimum))]
        print(nzd)
        print()
        hist = np.histogram(nzd, bins = nrainbin-1, range = (rain_minimum, rain_maximum))
        print(nzd)
        print()
        rain_bin = rain_bin + hist[0]
        print(nzd)
        print()

total = sum(rain_bin)
print(total)
print()
rain_bin[rain_bin < 1.0] = np.nan
rain_bin = rain_bin / (total * rain_step)
print(rain_bin)
print()

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.plot(rain_bin, 'ro')
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('Precipitation [mm/d]')
plt.ylabel('Event Normalized probability [1]')
plt.title('Precipitation daily PDF')
plt.show( )
exit()
