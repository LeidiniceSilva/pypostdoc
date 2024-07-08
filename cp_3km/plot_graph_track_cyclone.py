# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jan 02, 2024"
__description__ = "This script plot cyclone tracking"

import os
import numpy as np
import matplotlib.pyplot as plt

# Import cyclone information
cyclone = 'cyclone_i'

if cyclone == 'cyclone_i':
	dt = 'Jun2023'
	c = 'red'
	lat_cyclone = (-24,-24,-23,-25,-25,-25,-25,-27,-29,-30,-31,-32,-33,-34,-36,-37,-38,-39,-40,-40,-41,-41,-39,-39,-39,-40,-40,-42,-43,-43,-46,-48,-52,-55.5,-58)
	lon_cyclone = (-39,-41,-45,-45,-46,-46,-46,-47,-49,-49,-48,-47,-45,-44,-43,-40,-37,-36,-34,-33,-33,-32,-30,-27,-25,-24,-22,-21,-19,-16,-12,-6,-3,3,-1)
	ps_cyclone = (1019,1019,1019,1019,1016,1018,1016,1012,1010,1008,1008,1008,1008,1008,1005,1002,1002,1000,1000,995,995,995,1000,1000,1000,1000,1000,1000,995,990,990,985,975,970,965)
	dt_cyclone = ('00Z14','06Z14','12Z14','18Z14','00Z15','06Z15','12Z15','18Z15','00Z16','06Z16','12Z16','18Z16','00Z17','06Z17','12Z17','18Z17','00Z18','06Z18','12Z18','18Z18','00Z19','06Z19','12Z10','18Z19','00Z20','06Z20','12Z20','18Z20','00Z21','06Z21','12Z21','18Z21','00Z22','06Z22','12Z22')
else:
	dt= 'Jul2023'
	c = 'black'
	lat_cyclone = (-27,-29,-30,-31,-32,-33,-33,-35,-38,-39,-40,-41,-42,-43,-44,-46)
	lon_cyclone = (-58,-54,-53,-48,-46,-43,-39,-34,-31,-26,-23,-21,-17,-13,-9,-5)
	ps_cyclone = (1005,999,1000,995,995,990,995,990,995,995,993,990,989,985,985,980)
	dt_cyclone = ('12Z12','18Z12','00Z13','06Z13','12Z13','18Z13','00Z14','06Z14','12Z14','18Z14','00Z15','06Z15','12Z15','18Z15','00Z16','06Z16')

# Plot study area
fig = plt.figure() 
font_size = 10
  
plt.plot(dt_cyclone, ps_cyclone, linewidth=1, color=c, marker='o', markerfacecolor=c, markersize=3, label='Cyclone {}'.format(dt))
plt.xlabel(u'Hour/Day', fontsize=font_size, fontweight='bold')
plt.ylabel(u'Mean sea level pressure (hPa)', fontsize=font_size, fontweight='bold')
plt.xticks(rotation=90, fontsize=font_size)
plt.legend(loc=1, ncol=1, fontsize=font_size)
plt.grid(linestyle=':', color='gray', linewidth=1)

# Path out to save figure
path_out = '/marconi/home/userexternal/mdasilva/user/mdasilva/SAM-3km-cyclone/figs'
name_out = 'pyplt_maps_tracking_cyclone_era5_{0}.png'.format(dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

