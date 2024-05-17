# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 20, 2023"
__description__ = "This script plot cluster analysis"

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap

path = '/afs/ictp.it/home/m/mda_silv/Documents/BDMET'

skip_list = [1,2,415,19,21,23,28,35,41,44,47,54,56,59,64,68,7793,100,105,106,107,112,117,124,135,137,139,
149,152,155,158,168,174,177,183,186,199,204,210,212,224,226,239,240,248,249,253,254,276,277,280,293,298,
303,305,306,308,319,334,335,341,343,359,362,364,384,393,396,398,399,400,402,413,416,417,422,423,426,427,
443,444,446,451,453,457,458,467,474,479,483,488,489,490,495,505,509,513,514,516,529,534,544,559,566]

def import_inmet():
	
	ix = []		  
	iy = []
	clim_i = []
	clim_ii = []

	# Select lat and lon 
	for station in range(1, 567):
		if station in skip_list:
			continue
		if inmet[station][2] >= -11.25235:
			continue
			
		iy.append(inmet[station][2])
		ix.append(inmet[station][3])

		print('Reading weather station:', station, inmet[station][0])
		# Reading inmet
		d_i = xr.open_dataset('{0}/database/nc/hourly/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[station][0]))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		values_i = d_i.groupby('time.month').mean('time')
		values_i = values_i.values
		clim_i.append(values_i*24)
		
		values_ii = d_i.groupby('time.year').sum('time')
		values_ii = values_ii.valueszzz
		clim_ii.append(values_ii)		
			
	return iy, ix, clim_i, clim_ii
	

# Import latitude, longitude and database
lat_yy, lon_xx, clim_tot, clim_tot_i = import_inmet()			

longitude = [-41.106944, -52.211667, -52.875833, -47.625801, -41.457787, -41.090833, -41.488889, -55.525556, -40.741944, -47.212222, -40.687778, -53.224444, -47.523333, -39.616944, -44.303333, -42.035278, -52.245278, -49.528611, -49.49792, -46.949722, -42.021389, -54.013333, -46.008889, -43.767778, -45.026944, -48.545, -52.471389, -49.028877, -48.471389, -43.969444, -43.958611, -51.534202, -49.479722, -43.405278, -46.530556, -47.925833, -58.231389, -48.131111, -41.672222, -46.435556, -50.98552, -53.466944, -45.005556, -61.434167, -51.8175, -46.382996, -51.83424, -50.149722, -41.958333, -49.518133, -51.064042, -54.722615, -45.603889, -41.343889, -41.051667, -50.827222, -52.700833, -42.389167, -42.137222, -39.258611, -47.075278, -51.720833, -50.046389, -47.927614, -45.6175, -52.931944, -49.157733, -59.7625, -48.151574, -37.683889, -43.261111, -44.6175, -57.6375, -53.171389, -47.613056, -53.673597, -39.089444, -56.062951, -49.230556, -50.604167, -44.453785, -43.648269, -46.847222, -53.633114, -44.875, -53.095278, -54.618056, -45.593889, -51.552222, -43.282222, -40.736389, -49.914722, -52.305833, -38.9675, -44.416944, -48.62, -45.453889, -49.608333, -47.3825, -53.429444, -48.1375, -51.353611, -48.99, -49.220278, -50.141389, -41.977007, -42.749722, -42.943056, -47.199167, -53.766111, -49.051667, -44.011111, -53.111944, -43.213333, -47.545833, -39.181389, -51.077946, -49.268056, -39.6925, -51.148889, -41.864504, -37.795, -48.761944, -39.558056, -41.485556, -49.539444, -40.229444, -48.885833, -48.641667, -54.181944, -40.119722, -49.525, -49.191944, -49.646944, -47.775278, -50.853889, -53.822778, -53.375833, -50.595, -44.366389, -50.180556, -56.137778, -51.7175, -51.558889, -46.119444, -49.946389, -58.774722, -43.364329, -50.335556, -51.512845, -52.400833, -41.388889, -40.068611, -49.734722, -45.829722, -47.966944, -48.813333, -41.811944, -45.944389, -54.019167, -42.182778, -40.986505, -38.972222, -45.373056, -40.539722, -49.965111, -51.932778, -52.601111, -44.016111, -44.404167, -46.890278, -46.043333, -43.843889, -48.808611, -49.101698, -50.906389, -42.375833, -50.577778, -42.676944, -57.092222, -43.296944, -54.752222, -52.850278, -44.864444, -49.843333, -43.756111, -53.318056, -47.871944, -52.134444, -47.557417, -44.726944, -50.425556, -44.961944, -52.403582, -46.633916, -46.440556, -46.985935, -48.544444, -41.774125, -43.291389, -44.835556, -48.284167, -53.748056, -55.716389, -59.346111, -51.174722, -39.1825, -46.366389, -48.113889, -41.039444, -51.408611, -54.381111, -56.437115, -50.974722, -49.041944, -44.445, -54.528056, -44.040916, -43.190556, -43.402778, -43.411389, -43.595556, -50.145556, -52.167778, -49.574167, -52.381944, -42.435833, -50.965, -47.434167, -42.310278, -38.50576, -53.720556, -54.48, -48.185, -40.579444, -55.401389, -54.885653, -53.791111, -56.016389, -47.883927, -50.727778, -54.310833, -44.250928, -49.934722, -50.057778, -54.9625, -39.864167, -50.368889, -47.999167, -46.62, -45.123889, -47.114167, -50.633449, -54.694444, -58.763333, -42.608889, -51.870833, -43.684722, -40.25, -39.023056, -44.173333, -54.971944, -42.415556, -48.618056, -52.542387, -47.585556, -55.722778, -53.372222, -57.431667, -45.520833, -41.515423, -42.986944, -51.823333, -42.6225, -49.733333, -50.135833, -45.459836, -43.208333, -50.490251, -53.826667, -47.961944, -48.255556, -46.881944, -57.081944, -49.315278, -50.882778, -43.695556, -39.126667, -50.930278, -45.404167, -41.19, -50.210278, -42.864013, -40.403889, -60.157778, -40.305833, -40.801389, -52.39809]
latitude = [-20.104167, -14.016389, -20.444444, -15.596491, -15.751536, -19.532778, -20.750556, -29.709167, -20.636389, -11.284167, -16.166667, -17.339444, -14.133056, -13.009444, -22.975556, -16.848889, -15.902778, -12.592222, -28.931353, -19.605833, -22.975278, -31.347778, -20.031111, -21.228333, -12.124722, -20.559167, -21.75, -22.358052, -20.948056, -19.883889, -19.98, -29.164581, -28.126944, -13.251111, -22.951944, -15.789444, -12.521944, -15.599722, -14.181944, -15.524167, -26.819156, -30.545278, -22.688889, -11.445833, -16.966944, -21.918066, -30.807953, -29.049167, -21.5875, -19.53921, -29.674293, -20.447195, -22.750278, -21.714722, -22.041667, -29.368889, -31.403333, -17.705556, -19.735833, -17.739444, -21.780556, -19.1225, -24.78, -18.154779, -15.300278, -23.359167, -25.322464, -13.708056, -19.98586, -12.035833, -21.546667, -13.3325, -18.996667, -18.492778, -16.785, -28.60344, -12.675556, -15.559295, -25.448611, -27.288611, -18.747711, -18.231052, -11.594444, -26.286562, -20.173333, -25.699167, -31.0025, -19.481944, -21.457778, -22.589722, -18.291389, -17.336944, -27.657778, -12.196111, -19.885278, -27.6025, -20.455, -11.8875, -20.584444, -27.395556, -15.935278, -26.398611, -15.220278, -16.642778, -15.939722, -18.830354, -14.208056, -18.786944, -17.561389, -16.341667, -11.745, -20.031389, -28.653333, -12.193056, -24.671667, -14.658889, -25.567879, -26.913611, -14.171389, -16.423056, -11.328998, -11.2725, -26.950833, -17.006944, -16.575556, -14.979722, -15.244722, -23.981944, -26.081389, -23.449444, -13.527778, -18.952778, -18.409722, -27.418333, -20.359722, -25.010833, -22.300556, -32.534722, -20.165, -15.448056, -23.773333, -21.478611, -17.923611, -27.169167, -17.784444, -23.505278, -11.375, -21.769965, -27.802222, -28.222381, -25.371389, -12.557778, -19.356944, -21.666111, -12.1525, -16.260556, -28.604444, -22.376111, -21.680722, -24.533333, -20.263333, -18.78062, -13.906944, -22.314444, -19.407222, -22.235222, -23.405278, -17.454722, -15.085833, -14.408333, -13.253611, -22.861667, -16.686389, -25.508889, -17.745066, -31.248333, -21.105, -23.415278, -22.334722, -13.038611, -15.802778, -13.411111, -26.406389, -20.715, -23.000556, -20.556667, -27.920278, -12.615, -22.658333, -16.012222, -23.223611, -16.9625, -22.395833, -28.226805, -20.745237, -18.520556, -18.996684, -12.015278, -13.155681, -22.464722, -17.258064, -17.304167, -25.721944, -22.5525, -15.234444, -30.053611, -16.388889, -14.089167, -21.338333, -21.100833, -22.12, -15.58, -30.368578, -22.372778, -27.678611, -22.451389, -21.775, -22.653579, -22.988333, -22.94, -22.861389, -23.050278, -26.9375, -32.078889, -26.248611, -29.872222, -15.723056, -17.785278, -19.875278, -16.160278, -13.005515, -29.725, -27.890556, -11.428889, -19.988333, -30.750556, -29.191599, -27.854444, -28.65, -21.980353, -11.618889, -30.341389, -21.106502, -28.275556, -28.748611, -28.417222, -18.676111, -25.835556, -23.890833, -23.496389, -16.362778, -20.91, -18.969142, -29.702222, -13.303889, -22.871111, -28.704722, -22.757778, -17.798889, -11.664722, -19.455278, -20.981667, -22.645833, -16.679722, -28.859211, -23.426111, -12.555, -33.742222, -14.65, -23.041667, -17.89284, -22.448611, -29.449167, -19.573889, -29.350278, -30.010278, -18.200855, -22.098333, -21.927251, -29.089444, -19.71, -18.916944, -16.554167, -29.84, -28.5325, -28.513611, -22.358056, -13.343611, -21.319167, -21.566389, -20.385556, -24.280278, -20.762607, -20.466944, -12.735, -20.270833, -14.886389, -26.938666]
list_hc = [2, 3, 2, 3, 2, 2, 3, 0, 1, 3, 2, 3, 3, 4, 1, 2, 3, 1, 0, 3, 0, 0, 3, 3, 2, 2, 2, 2, 3, 1, 1, 0, 1, 2, 3, 3, 1, 1, 2, 2, 0, 0, 3, 1, 3, 3, 0, 0, 2, 3, 0, 2, 3, 2, 2, 0, 0, 2, 3, 4, 2, 3, 2, 1, 3, 0, 0, 1, 3, 4, 3, 2, 2, 3, 1, 0, 4, 3, 0, 0, 3, 3, 1, 0, 3, 0, 0, 3, 2, 1, 2, 3, 0, 4, 1, 0, 3, 3, 3, 0, 3, 0, 3, 1, 1, 2, 2, 2, 3, 3, 3, 1, 0, 2, 1, 4, 0, 0, 4, 1, 2, 4, 0, 4, 2, 3, 2, 2, 1, 0, 4, 3, 3, 0, 3, 2, 2, 0, 2, 2, 2, 2, 3, 0, 3, 2, 1, 1, 0, 0, 0, 4, 2, 2, 1, 3, 0, 2, 3, 0, 3, 2, 4, 3, 2, 3, 0, 3, 2, 2, 3, 3, 2, 0, 3, 0, 3, 2, 1, 1, 2, 1, 0, 3, 2, 3, 0, 3, 2, 1, 1, 1, 2, 0, 3, 3, 1, 3, 2, 1, 3, 3, 0, 0, 3, 0, 4, 3, 2, 2, 0, 1, 0, 0, 0, 3, 2, 3, 2, 2, 2, 2, 0, 0, 0, 0, 2, 3, 3, 2, 4, 0, 0, 3, 2, 0, 0, 0, 0, 3, 1, 0, 1, 0, 0, 0, 2, 2, 2, 3, 2, 3, 3, 0, 1, 2, 0, 2, 2, 4, 1, 2, 1, 1, 0, 2, 1, 0, 2, 3, 2, 1, 0, 3, 0, 0, 3, 3, 3, 0, 1, 3, 3, 0, 0, 0, 2, 0, 2, 3, 3, 2, 3, 2, 1, 2, 2, 0]
 
print(len(longitude))
print(len(latitude))
print(len(list_hc))

count_i = []
count_ii = []
count_iii = []
count_iv = []
count_v = []

for count, idx in enumerate(list_hc):
	
	if idx == 0:
		count_i.append(count)

	if idx == 1:
		count_ii.append(count)
		
	if idx == 2:
		count_iii.append(count)
	
	if idx == 3:
		count_iv.append(count)
	
	if idx == 4:
		count_v.append(count)

lon_c_i = []
lon_c_ii = []
lon_c_iii = []
lon_c_iv = []
lon_c_v = []

lat_c_i = []
lat_c_ii = []
lat_c_iii = []
lat_c_iv = []
lat_c_v = []

pre_c_i = []
pre_c_ii = []
pre_c_iii = []
pre_c_iv = []
pre_c_v = []

pre_c_i_i = []
pre_c_ii_i = []
pre_c_iii_i = []
pre_c_iv_i = []
pre_c_v_i = []

for c_i in count_i:
	lon_c_i.append(longitude[c_i])
	lat_c_i.append(latitude[c_i])
	pre_c_i.append(clim_tot[c_i])
	pre_c_i_i.append(clim_tot_i[c_i])	

for c_ii in count_ii:
	lon_c_ii.append(longitude[c_ii])
	lat_c_ii.append(latitude[c_ii])
	pre_c_ii.append(clim_tot[c_ii])
	pre_c_ii_i.append(clim_tot_i[c_ii])

for c_iii in count_iii:
	lon_c_iii.append(longitude[c_iii])
	lat_c_iii.append(latitude[c_iii])
	pre_c_iii.append(clim_tot[c_iii])
	pre_c_iii_i.append(clim_tot_i[c_iii])
	
for c_iv in count_iv:
	lon_c_iv.append(longitude[c_iv])
	lat_c_iv.append(latitude[c_iv])
	pre_c_iv.append(clim_tot[c_iv])
	pre_c_iv_i.append(clim_tot_i[c_iv])
		
for c_v in count_v:
	lon_c_v.append(longitude[c_v])
	lat_c_v.append(latitude[c_v])
	pre_c_v.append(clim_tot[c_v])
	pre_c_v_i.append(clim_tot_i[c_v])

cluster_i = np.nanmean(pre_c_i, axis=0)
cluster_ii = np.nanmean(pre_c_ii, axis=0)
cluster_iii = np.nanmean(pre_c_iii, axis=0)
cluster_iv = np.nanmean(pre_c_iv, axis=0)
cluster_v = np.nanmean(pre_c_v, axis=0)

cluster_i_i = np.nanmean(pre_c_i_i, axis=0)
cluster_ii_i = np.nanmean(pre_c_ii_i, axis=0)
cluster_iii_i = np.nanmean(pre_c_iii_i, axis=0)
cluster_iv_i = np.nanmean(pre_c_iv_i, axis=0)
cluster_v_i = np.nanmean(pre_c_v_i, axis=0)

# Plot figure   
fig = plt.figure(figsize=(10, 7))

ax = fig.add_subplot(1, 2, 1) 
my_map = Basemap(projection='cyl', llcrnrlon=-85., llcrnrlat=-60., urcrnrlon=-30.,urcrnrlat=15., resolution='c')
my_map.drawmeridians(np.arange(-85.,-30.,10.), labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-60.,15.,10.), labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.readshapefile('/afs/ictp.it/home/m/mda_silv/Documents/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='black', linewidth=.5)

my_map.plot(lon_c_i, lat_c_i, 'o', color='blue', label='Cluster I', markersize=2)
my_map.plot(lon_c_ii, lat_c_ii, 'o', color='gray', label='Cluster II', markersize=2)
my_map.plot(lon_c_iii, lat_c_iii, 'o', color='green', label='Cluster III', markersize=2)
my_map.plot(lon_c_iv, lat_c_iv, 'o', color='red', label='Cluster IV', markersize=2)
my_map.plot(lon_c_v, lat_c_v, 'o', color='yellow', label='Cluster V', markersize=2)
plt.title('(a)', loc='left', fontsize=10, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=20, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=30, fontweight='bold')
plt.text(-35, 10, u'\u25B2 \nN', fontsize=10, fontweight='bold')
plt.legend(loc=3, ncol=3, fontsize=8, frameon=False)

# SESA
a1,b1 = (-78,-35)
a2,b2 = (-78,-11)
a3,b3 = (-35,-11)
a4,b4 = (-35,-35)
poly1 = Polygon([(a1,b1),(a2,b2),(a3,b3),(a4,b4)], facecolor='none', edgecolor='red', linewidth=1.)
plt.gca().add_patch(poly1)

ax = fig.add_subplot(2, 2, 2) 
time = np.arange(0.5, 12 + 0.5)
plt.plot(time, cluster_i, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='blue', label = 'Cluster I')
plt.plot(time, cluster_ii, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='gray', label = 'Cluster II')
plt.plot(time, cluster_iii, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='green', label = 'Cluster III')
plt.plot(time, cluster_iv, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='red', label = 'Cluster IV')
plt.plot(time, cluster_v, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='yellow', label = 'Cluster V')
plt.title('(b)', loc='left', fontsize=10, fontweight='bold')
plt.ylabel('Precipitation (mm d⁻¹)', fontsize=10, fontweight='bold')
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))
plt.yticks(np.arange(0, 15, 1))
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
ax.tick_params(axis='both', which='major', labelsize=10)
plt.legend(loc=9, ncol=3, fontsize=10, handletextpad=0.2, frameon=False)
	
ax = fig.add_subplot(2, 2, 4) 
box_plot_data=[cluster_i_i, cluster_ii_i, cluster_iii_i, cluster_iv_i, cluster_v_i]
box = plt.boxplot(box_plot_data, patch_artist=True, labels=['Cluster I','Cluster II','Cluster III','Cluster IV', 'Cluster V'])
colors = ['blue', 'gray', 'green', 'red', 'yellow']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
for median in box['medians']:
    median.set(color='k', linewidth=1.)
plt.title('(c)', loc='left', fontsize=10, fontweight='bold')
plt.ylabel('Annual mean precipitation (mm)', fontsize=10, fontweight='bold')
plt.yticks(np.arange(200, 2100, 100))
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
ax.tick_params(axis='both', which='major', labelsize=10)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_cluster_analysis_sesa.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

