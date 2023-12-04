# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 16, 2023"
__description__ = "This script plot cluster analysis"

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap

path = '/afs/ictp.it/home/m/mda_silv/Documents/FPS_SESA'


def import_inmet():
	
	ix = []		  
	iy = []
	clim_i = []
	clim_ii = []

	# Select lat and lon 
	for i in range(1, 100):
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])

		print('Reading weather station:', i, inmet[i][0], inmet[i][1])
		# Reading inmet
		d_i = xr.open_dataset('{0}/database/obs/inmet/inmet_nc_sesa/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_i = d_i.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		clim_i.append(values_i*24)
		
		d_ii = xr.open_dataset('{0}/database/obs/inmet/inmet_nc_sesa/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_ii = d_ii.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.groupby('time.year').sum('time')
		values_ii = d_ii.values
		clim_ii.append(values_ii)		
			
	return iy, ix, clim_i, clim_ii


def import_smn_i():
	
	iy = []
	ix = []
	clim_i = []
	clim_ii = []
	
	# Select lat and lon 
	for i in range(1, 73):
		iy.append(smn_i[i][1])
		ix.append(smn_i[i][2])
		
		print('Reading weather station:', i, smn_i[i][0])	
		# Reading smn
		d_i = xr.open_dataset('{0}/database/obs/smn_i/smn_nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[i][0]))
		d_i = d_i.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		clim_i.append(values_i*24)
				
		d_ii = xr.open_dataset('{0}/database/obs/smn_i/smn_nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[i][0]))
		d_ii = d_ii.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.groupby('time.year').sum('time')
		values_ii = d_ii.values
		clim_ii.append(values_ii)	
		
	return iy, ix, clim_i, clim_ii
	

def import_smn_ii():
	
	iy = []
	ix = []
	clim_i = []
	clim_ii = []
	
	# Select lat and lon 
	for i in range(1, 86):
		iy.append(smn_ii[i][1])
		ix.append(smn_ii[i][2])
		
		print('Reading weather station:', i, smn_ii[i][0])	
		# Reading smn
		d_i = xr.open_dataset('{0}/database/obs/smn_ii/smn_nc/'.format(path) + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[i][0]))
		d_i = d_i.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		clim_i.append(values_i)
				
		d_ii = xr.open_dataset('{0}/database/obs/smn_ii/smn_nc/'.format(path) + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[i][0]))
		d_ii = d_ii.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.groupby('time.year').sum('time')
		values_ii = d_ii.values
		clim_ii.append(values_ii)	
		
	return iy, ix, clim_i, clim_ii
	
	
var = 'pr'

# Import latitude, longitude and database
lat_x, lon_x, clim_inmet_x, clim_inmet_xx = import_inmet()			
lat_y, lon_y, clim_smn_y, clim_smn_yy = import_smn_i()
lat_z, lon_z, clim_smn_z, clim_smn_zz = import_smn_ii()

lon_xx = lon_x + lon_y + lon_z
lat_yy = lat_x + lat_y + lat_z

clim_tot = clim_inmet_x + clim_smn_y + clim_smn_z
clim_tot_i = clim_inmet_xx + clim_smn_yy + clim_smn_zz

longitude = [-53.37138888, -53.37638888, -52.1, -52.70111111, -54.01333333, -54.61805554, -55.61277777, -51.83472221, -53.46749999, -56.43722221, -54.31083333, -51.16666666, -53.82666666, -52.38249999, -57.08249999, -55.5261111, -53.70361111, -54.69416666, -51.06416666, -50.82749999, -49.73305554, -54.88555555, -51.53472221, -50.14944444, -49.49805555, -52.55833333, -50.05833333, -51.87055554, -53.11194444, -56.01555555, -48.81333333, -53.6736111, -49.31555555, -50.88277777, -54.96222222, -49.93444444, -52.40388888, -51.51222222, -49.48333333, -53.31694443, -54.48055554, -53.7911111, -50.33555555, -49.04194444, -52.30638888, -48.61999999, -49.64666666, -53.41277777, -50.60416666, -51.57555554, -48.76194444, -52.415, -49.26749999, -50.98555555, -52.34888888, -51.35444444, -53.63277777, -49.58055554, -50.3686111, -51.08944444, -48.80861111, -49.26666666, -52.39194444, -49.1575, -50.8711111, -54.01972221, -48.88527777, -48.16444444, -50.18055554, -49.94638888, -54.18166666, -51.91666666, -52.93194443, -52.13444444, -50.97416666, -49.02888888, -49.965, -51.4, -50.49027777, -54.52805554, -49.73416666, -51.55222222, -48.11388888, -50.93027777, -48.54472221, -54.6, -52.87555554, -49.96583333, -48.15138888, -49.51805554, -57.63749999, -48.25555555, -49.19194444, -51.71777777, -49.1, -52.60111111, -53.22416666, -49.91472222, -48.28416666, -56.542211, -56.224017, -55.866233, -56.555522, -57.4559, -57.485838, -57.097101, -56.618974, -56.456119, -58.165178, -57.325281, -57.78309, -56.972142, -57.429814, -56.214658, -56.326675, -58.295964, -56.657495, -56.629005, -58.141889, -57.021887, -57.85945, -57.680318, -56.230152, -56.905948, -57.421509, -58.070281, -57.903, -57.055016, -54.193437, -58.247317, -57.880716, -57.913847, -55.642269, -56.635879, -57.703535, -56.787211, -57.08769, -58.468999, -56.411658, -57.743836, -57.674446, -57.820722, -57.966907, -57.982836, -57.622519, -56.786308, -57.291694, -57.082856, -57.074737, -56.66754, -57.124992, -56.84021, -53.933283, -58.002561, -57.880984, -57.120524, -57.939, -56.392531, -57.977472, -57.916594, -58.317145, -55.129216, -57.565952, -57.92773, -56.028995, -56.215761, -56.879045, -58.193534, -57.440628, -57.480465, -56.81167, -58.42, -58.48, -65.77, -61.95, -58.02, -64.2, -64.18, -58.77, -58.6, -58.53, -58.23, -58.62, -54.47, -68.75, -65.08, -60.92, -57.9, -65.6, -66.82, -63.37, -60.58, -62.15, -68.78, -68.85, -57.65, -64.32, -60.48, -57.15, -63.88, -55.97, -59.7, -59.05, -64.23, -60.78, -65.48, -68.42, -66.35, -68.42, -68.4, -64.3, -60.82, -63.81, -67.57, -65.13, -65.38, -58.92, -58.67, -55.43, -58.93, -65.42, -69.12, -61.1, -62.12, -58.02, -60.54, -60.55, -61.55, -59.7, -60.45, -59.8, -60.88, -56.51, -56.01, -57.77, -56.5, -54.19, -58.07, -58.04, -56.21, -55.54, -54.31, -57.97, -54.39, -56.19, -54.95, -56.76, -56.92, -57.65, -65.1, -62.73, -64.23, -56.46, -54.35, -57.51, -56.35]
latitude = [-33.74166666, -32.53749999, -32.03333333, -31.40583333, -31.34777777, -31.0025, -30.84249999, -30.81055555, -30.54777777, -30.3686111, -30.34138888, -30.05, -29.89388888, -29.87333332, -29.84249999, -29.71166666, -29.70833333, -29.70194444, -29.67444443, -29.36888888, -29.35027777, -29.19138888, -29.16722221, -29.04888888, -28.93138888, -28.85361111, -28.75138888, -28.70472222, -28.65333333, -28.64944444, -28.60416666, -28.60361111, -28.53249999, -28.51361111, -28.41694443, -28.27555554, -28.22944443, -28.22194443, -28.13333333, -27.92166666, -27.89305555, -27.85416666, -27.80222222, -27.6786111, -27.66027777, -27.6025, -27.41833332, -27.39555555, -27.2886111, -27.16944443, -26.95083333, -26.9386111, -26.91361111, -26.81944443, -26.41722221, -26.3986111, -26.28638888, -26.2486111, -25.83499999, -25.56583333, -25.51277777, -25.43333333, -25.36888888, -25.32222221, -25.01333333, -24.53583333, -23.98138888, -23.85138888, -23.77305554, -23.50527777, -23.44944444, -23.40777777, -23.37583332, -22.49166666, -22.37249999, -22.35805555, -22.23527777, -22.11999999, -21.92722221, -21.77472221, -21.66527777, -21.45777777, -21.33833333, -21.31916666, -20.55888888, -20.45, -20.44444444, -20.40305555, -19.98583333, -19.5391666, -18.99666666, -18.91722221, -18.40972222, -17.92388888, -17.71666667, -17.45472222, -17.33944444, -17.33694444, -17.30416666, -28.416682, -27.416682, -26.416682, -25.416682, -24.416682, -23.416682, -22.416682, -21.416682, -20.416682, -19.416682, -18.416682, -30.787252, -30.299302, -29.819517, -30.751196, -30.637259, -29.98349, -30.62427, -31.337833, -30.776241, -31.09242, -30.014261, -30.580333, -30.58639, -30.613837, -31.444407, -29.787069, -31.459, -30.755734, -27.298582, -30.971026, -30.971187, -30.992189, -28.177515, -31.037601, -31.286127, -30.432802, -29.05663, -30.590296, -30.865546, -29.660671, -29.845637, -30.690307, -30.214445, -30.626772, -30.249296, -30.787214, -30.289, -29.721745, -30.205281, -30.223953, -30.476397, -31.444635, -27.153304, -31.404212, -30.390258, -31.493837, -32.135, -31.350331, -31.38978, -31.272706, -30.337372, -27.869343, -29.296744, -30.91468, -28.545073, -30.988557, -31.004877, -29.377358, -30.431886, -28.691829, -29.470444, -34.57, -34.58, -28.6, -29.88, -31.3, -31.3, -31.4, -27.45, -34.6, -34.82, -26.2, -33.0, -25.73, -30.23, -24.38, -34.55, -34.97, -22.1, -29.38, -34.13, -24.7, -32.7, -32.83, -32.88, -30.27, -23.15, -31.78, -29.68, -31.67, -27.37, -29.18, -27.45, -33.12, -32.92, -24.85, -31.57, -33.27, -33.08, -34.58, -27.77, -31.7, -22.65, -28.07, -31.95, -33.73, -28.43, -34.61, -27.65, -27.42, -27.05, -33.73, -27.1, -32.68, -29.17, -31.85, -33.93, -31.18, -29.18, -26.87, -33.74, -33.02, -30.38, -34.83, -34.45, -33.35, -32.37, -33.25, -32.35, -34.86, -30.9, -34.49, -31.4, -33.71, -34.09, -34.97, -34.35, -33.54, -32.69, -26.85, -34.92, -28.02, -24.67, -24.03, -25.24, -26.18]
list_hc = [4, 4, 4, 4, 4, 0, 0, 4, 4, 0, 0, 4, 0, 4, 0, 0, 0, 0, 4, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 3, 0, 1, 1, 0, 1, 0, 1, 3, 3, 3, 3, 1, 1, 1, 3, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 2, 4, 0, 4, 4, 4, 3, 4, 4, 4, 1, 4, 4, 2, 4, 4, 3, 1, 3, 3, 3, 3, 0, 3, 3, 3, 3, 1, 3, 3, 3, 3, 3, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 0, 2, 2, 1, 3, 3, 1, 3, 0, 2, 1, 3, 3, 2, 2, 2, 2, 2, 2, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 3, 0, 1, 1, 2, 1, 2, 0, 2, 2, 2, 1, 2, 2, 2, 0, 4, 3, 4, 4, 3, 4, 4, 0, 4, 0, 4, 4, 4, 3, 4, 4, 1, 2, 2, 0, 1, 1, 0]
 
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
my_map = Basemap(ax=ax, projection='cyl', llcrnrlon=-70, llcrnrlat=-40, urcrnrlon=-45,urcrnrlat=-15, resolution='c')
my_map.drawmeridians(np.arange(-70,-45,5), labels=[0,0,0,1], size=10, linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-40,-15,5), labels=[1,0,0,0], size=10, linewidth=0.5, color='black') 
my_map.readshapefile('/afs/ictp.it/home/m/mda_silv/Documents/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='black', linewidth=.5)

my_map.plot(lon_c_i, lat_c_i, 'o', color='blue', label='Cluster I', markersize=5)
my_map.plot(lon_c_ii, lat_c_ii, 'o', color='gray', label='Cluster II', markersize=5)
my_map.plot(lon_c_iii, lat_c_iii, 'o', color='green', label='Cluster III', markersize=5)
my_map.plot(lon_c_iv, lat_c_iv, 'o', color='red', label='Cluster IV', markersize=5)
my_map.plot(lon_c_v, lat_c_v, 'o', color='yellow', label='Cluster V', markersize=5)
plt.title('(a)', loc='left', fontsize=10, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=20, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=30, fontweight='bold')
plt.text(-68, -19, u'\u25B2 \nN', fontsize=10, fontweight='bold')
plt.legend(loc=3, ncol=3, fontsize=8, frameon=False)

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
plt.yticks(np.arange(0, 11, 1))
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
plt.yticks(np.arange(200, 1800, 100))
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
ax.tick_params(axis='both', which='major', labelsize=10)

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_cluster_analysis_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

