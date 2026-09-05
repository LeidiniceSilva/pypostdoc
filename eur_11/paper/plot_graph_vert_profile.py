# -*- coding: utf-8 -*-

# author       = "Leidinice Silva"
# email        = "leidinicesilva@gmail.com"
# date         = "Jul 28, 2026"
# description  = "This script plots vertical profile"

import os
import glob
import netCDF4
import numpy as np
import matplotlib.pyplot as plt


def import_obs(path, obs_var, obs_name, dt, tag):
    # Fixed to use obs_var (e.g., 'clwc') instead of param
    arq = f"{path}/{obs_var}_{obs_name}_{dt}_{tag}.nc"
    files = glob.glob(arq)
    if not files:
        raise FileNotFoundError(f"Missing OBS file: {arq}")
    
    data = netCDF4.Dataset(files[0])
    var = data.variables[obs_var][:]
    
    # Check dimensions before slicing
    if var.ndim == 4:
        value = var[:, ::-1, :, :]
        mean = np.nanmean(np.nanmean(value, axis=2), axis=2)[0]
    else:
        mean = np.nanmean(var, axis=0) if var.ndim > 1 else var

    return mean


def import_rcm(path, param, exp_name, dt, tag):
    arq = f"{path}/{param}_RegCM5_{exp_name}_{dt}_{tag}.nc"
    files = glob.glob(arq)
    if not files:
        raise FileNotFoundError(f"Missing RCM file: {arq}")
        
    data = netCDF4.Dataset(files[0])
    var = data.variables[param][:]
    
    if var.ndim == 4:
        value = var[:, :, :, :]
        mean = np.nanmean(np.nanmean(value, axis=2), axis=2)[0]
    else:
        mean = np.nanmean(var, axis=0) if var.ndim > 1 else var

    return mean


# Year
yr = "2000-2009"

# Paths and setup
path_data = "/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11/postproc/paper/vert"
path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/EUR-11/figs/paper'
os.makedirs(path_out, exist_ok=True)

seasons = ['DJF', 'MAM', 'JJA', 'SON', 'ANN']
exps = [('NoTo-EUR', 'red', 'NoTo'), ('WSM5-EUR', 'blue', 'WSM5'), 
        ('WSM7-EUR', 'green', 'WSM7'), ('WDM7-EUR', 'orange', 'WDM7')]

vars_info = [
    ('cl',  'cc',   100,  '(a)', '(b)', '(c)', '(d)', '(e)'),
    ('cli', 'ciwc', 1e6,  '(f)', '(g)', '(h)', '(i)', '(j)'),
    ('clw', 'clwc', 1e6,  '(k)', '(l)', '(m)', '(n)', '(o)')
]

dict_plot = {
    'cl':  ['Cloud fraction (%)', 0, 30, np.arange(0, 33, 3)],
    'cli': ['Cloud liquid ice (mg kg$^-$$^1$)', 0, 10, np.arange(0, 11, 1)],
    'clw': ['Cloud liquid water (mg kg$^-$$^1$)', 0, 75, np.arange(0, 80, 5)]
}

levels_i = (1000,975,950,925,900,875,850,825,800,775,750,700,650,600,550,500,450,400,350,300,250,225,200,175,150,125,100,70,50,30,20,10,7,5,3,2,1)
levels_ii = (1000,925,850,700,600,500,400,300,250,200,150,100)

font_size = 8
fig = plt.figure(figsize=(12, 10))
fig.patch.set_facecolor('#E0E0E0')
fig.patch.set_alpha(0.75)

# Plotting with loops
idx = 1
for var_name, obs_var, scale, *tags in vars_info:
    for s_idx, season in enumerate(seasons):
        ax = fig.add_subplot(3, 5, idx)
        
        # Plot ERA5
        era5 = import_obs(path_data, obs_var, "ERA5", f"{season}_{yr}", "land_box") * scale
        plt.plot(era5, levels_i, color='black', label='ERA5', linewidth=1)
        
        # Plot RCM Exps
        for exp_code, color, label in exps:
            rcm = import_rcm(path_data, var_name, exp_code, f"{season}_{yr}", "land_box") * scale
            plt.plot(rcm, levels_ii, color=color, label=label, linewidth=1)
            
        plt.title(f"{tags[s_idx]} {season}", loc='left', fontsize=font_size, fontweight='bold')
        plt.xlabel(dict_plot[var_name][0], fontsize=font_size, fontweight='bold')
        plt.ylabel('Level pressure (hPa)', fontsize=font_size, fontweight='bold')
        plt.xlim(dict_plot[var_name][1], dict_plot[var_name][2])
        plt.ylim(0, 1000)
        plt.yticks(fontsize=font_size)
        plt.xticks(dict_plot[var_name][3], fontsize=font_size)
        plt.grid(linestyle='--')
        plt.gca().invert_yaxis()
        
        if idx == 1:
            plt.legend(loc=1, ncol=1, fontsize=font_size)
            
        idx += 1

# Save
name_out = f'pyplt_vert_profile_RegCM5_EUR-11_{yr}.png'
plt.tight_layout()
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
