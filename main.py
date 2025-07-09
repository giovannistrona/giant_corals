from netCDF4 import Dataset
import csv,os
import wget
from numpy.ma import is_masked
import numpy as np
import random
from rasterio.transform import from_origin
import rasterio
from datetime import datetime, timedelta
import geopandas as gpd
from shapely.geometry import Point
from pyproj import Transformer
from random import randrange


os.makedirs('dhw_mon_layers',exist_ok=True)

##download historical DHW data
for year in range(1986,2022):
	url = 'https://www.star.nesdis.noaa.gov/pub/sod/mecb/crw/data/5km/v3.1_op/nc/v1.0/annual/ct5km_dhw-max_v3.1_'+str(year)+'.nc'
	filename = wget.download(url,out='./dhw_mon_layers')
	print (year)


def cell_to_coord(col, row, cellx, celly, xmin, ymax):
	lon = cellx*col + xmin
	lat = ymax-celly*row
	return lat,lon


def pre_months_fun(m, y):
    date_start = datetime(y, m, 1)
    return [
        ((date_start - timedelta(days=30 * i)).month, (date_start - timedelta(days=30 * i)).year)
        for i in range(12)
    ][::-1]


def random_point_in_polygon(polygon, max_att=1000):
    """Generates a random point within a Shapely polygon."""
    min_x, min_y, max_x, max_y = polygon.bounds
    attempts = 0
    while attempts < max_att:
        x = random.uniform(min_x, max_x)
        y = random.uniform(min_y, max_y)
        random_point = Point(x, y)
        if polygon.contains(random_point):
            return random_point
        attempts += 1
    return None  # Se fallisce dopo max_att tentativi


def get_rand_point_radius(original_point,d): #point in meters, distance in meters
	theta = np.random.uniform(0, 2 * np.pi)
	x_new = original_point.x + d * np.cos(theta)
	y_new = original_point.y + d * np.sin(theta)
	new_point = Point(x_new, y_new)
	lon, lat = transformer.transform(new_point.x, new_point.y)
	return (lat,lon)


def point_to_cell(point_x, point_y, cellx, celly, xmin, ymax):
    col = int((point_x - xmin) / cellx)
    row = int((point_y - ymax) / -celly)
    return row,col


crs_to = "EPSG:4326"
metric_crs = "EPSG:3857"  # Web Mercator (o scegli un UTM locale)


transformer = Transformer.from_crs(metric_crs, crs_to, always_xy=True)

# realm_pols = gpd.read_file("./MEOW/meow_ecos.shp")
# realm_pols = realm_pols.to_crs(metric_crs)


giants_sites = list(csv.reader(open('MTG_data_16_04.csv','r')))
#giants_sites[0]
#['lat', 'lon', 'depth', 'size', 'Bleaching', 'Genus', 'health']

giants_sites = giants_sites[1:]
mean_d = -np.mean([float(i[2]) for i in giants_sites if i[2]!='NA'])
giants_rep_id_lat_lon_depth = [[0,i,float(giants_sites[i][0]),float(giants_sites[i][1]),-float(giants_sites[i][2]) if giants_sites[i][2]!='NA' else mean_d] for i in range(len(giants_sites))]

gc_points = gpd.GeoDataFrame(geometry=[Point([i[3],i[2]]) for i in giants_rep_id_lat_lon_depth],crs="EPSG:4326")
gc_points["ID"] = range(len(gc_points))
gc_points = gc_points.to_crs(metric_crs)


###get depth of giant colonies; sample random localities at the same depth within a range of x km
cellx,celly = 0.004166666666666667,0.004166666666666667
xmin = -180
ymax = 90

#Download gebco bathymetry data at https://www.bodc.ac.uk/data/open_download/gebco/gebco_2024_sub_ice_topo/geotiff/
depth = rasterio.open('gebco_depth_24.tif','r').read(1)
rand_points = []
for i in range(len(giants_rep_id_lat_lon_depth)):
	r0, id, lat, lon, depth_i = giants_rep_id_lat_lon_depth[i]
	original_point = gc_points['geometry'][i]
	for rep in range(1,1001):
		stop = 'no'
		att = 0
		while stop=='no':
			att+=1
			rand_dist = randrange(20000, 200000)
			r_lat,r_lon = get_rand_point_radius(original_point,rand_dist)
			r_y,r_x = point_to_cell(r_lon, r_lat, cellx, celly, xmin, ymax)
			depth_r = depth[r_y][r_x]
			min_depth,max_depth = depth_i-5,depth_i+5
			if max_depth>0:
			 	max_depth = -1
			if min_depth<=depth_r<=max_depth:
				r_lat,r_lon = cell_to_coord(r_x, r_y, cellx, celly, xmin, ymax)
				if [rep,id,r_lat,r_lon,depth_r] not in rand_points:
					rand_points.append([rep,id,r_lat,r_lon,depth_r])
					stop = 'yes'
		print(rep,att)



giants_rep_id_lat_lon_depth+=rand_points
R,C = Dataset('./dhw_mon_layers/ct5km_dhw-max_v3.1_199801.nc').variables['degree_heating_week'][:][0].shape


months = ["{:02d}".format(i) for i in range(1,13)]
####

cellx,celly = 0.05,0.05 ###about 5km
xmin = -180
ymax = 90
months = ["{:02d}".format(i) for i in range(1,13)]
giant_sites_time_series = []
for year in range(1986,2025):
	for month in months:
		data = Dataset('./dhw_mon_layers/ct5km_dhw-max_v3.1_'+str(year)+str(month)+'.nc')
		dhw = data.variables['degree_heating_week'][:][0]
		rows,cols = dhw.shape
		for rep,id,y,x,depth_i in giants_rep_id_lat_lon_depth:
			start_row,start_col = point_to_cell(x, y, cellx, celly, xmin, ymax)
			dhw_nn,distance = dhw[start_row][start_col],0
			if is_masked(dhw_nn):
				dhw_nn = 0
				distance = 1
				while distance<10:
					found = False
					for i in range(-distance, distance + 1):
						for j in range(-distance, distance + 1):
							if (abs(i) == distance or abs(j) == distance):  # only border of square
								row = start_row + i
								col = start_col + j
								if 0 <= row < rows and 0 <= col < cols:
									if not np.ma.is_masked(dhw[row, col]) and not np.isnan(dhw[row, col]):
										dhw_nn,distance= dhw[row, col], distance
										found = True
										break
						if found:
							break
					if found:
						break
					distance += 1
			giant_sites_time_series.append([rep,id,y,x,year,month,dhw_nn,depth_i])
	print (giant_sites_time_series[-1])



out = open('giant_sites_time_series.csv','w')
out.write('rep,site_id,lat,lon,year,month,dhw_max,depth\n')
for i in giant_sites_time_series:
	o = out.write(','.join(map(str,i))+'\n')


out.close()



####################################################################################
###identify p function to make values comparable between the high res historical data and the lower res future data

#####Future projections
####first convert nc file to monthly dhw max tiffs
try:
	os.mkdir('./future_dhw_max_mon/')
except:
	pass


data = Dataset('ens5_ssp245_DHW_1985_2100.nc') #download the dataset here: https://adelaide.figshare.com/articles/dataset/Global_projections_of_sea_surface_temperature_and_coral_bleaching_risk_in_the_21st_century/25143128?file=44413187
rows, cols = data.variables['DHW'][0].shape
cellx,celly = 0.5,0.5
xmin = -180
ymax = 35
nodata_value = -999
transform = from_origin(xmin, ymax, cellx, celly)

start_day = int(data.variables['time'][0].data) ####days since 1985-01-01

start_date = datetime(1985, 1, 1)
y_m = [tuple(map(int,(start_date + timedelta(days=int(i))).strftime("%Y-%m").split('-'))) for i in data.variables['time']]


###make dictionary with ids of layers per month
m_ids_dict = dict()
for i in range(len(y_m)):
	m_ids_dict[y_m[i]] = [i]+m_ids_dict.get(y_m[i],[])


#####get month max DHW
rows, cols = data.variables['DHW'][0].shape
months = list(range(1,13))

for year in range(1986,2101):
	for month in months:
		lays = [data.variables['DHW'][lay_id] for lay_id in m_ids_dict[(year,month)]]
		max_matrix = np.ma.maximum.reduce(lays,0)
		max_matrix = np.flipud(max_matrix)
		with rasterio.open(
			'./future_dhw_max_mon/'+str(year)+'_'+str(month)+'.tif',
			"w",
			driver="GTiff",
			height=rows,
			width=cols,
			count=1,
			dtype=max_matrix.dtype,
			crs="EPSG:4326",
			transform=transform,
			nodata=nodata_value
		) as dst:
			dst.write(max_matrix, 1)
	print (year,month)




months = list(range(1,13))
giant_sites_time_series = []
for year in range(1995,2101):
	for month in months:
		dhw = rasterio.open('./future_dhw_max_mon/'+str(year)+'_'+str(month)+'.tif').read(1)
		for rep,id,y,x,depth_i in giants_rep_id_lat_lon_depth:
			if rep==0:
				start_row,start_col = point_to_cell(x, y, cellx, celly, xmin, ymax)
				dhw_nn,distance = dhw[start_row][start_col],0
				if dhw_nn==-999:
					dhw_nn == 0
					distance = 1
					while distance<100:
						found = False
						for i in range(-distance, distance + 1):
							for j in range(-distance, distance + 1):
								if (abs(i) == distance or abs(j) == distance):  # only border of square
									row = start_row + i
									col = start_col + j
									if 0 <= row < rows and 0 <= col < cols:
										if dhw[row, col]!=-999 and not np.isnan(dhw[row, col]):
											dhw_nn,distance= dhw[row, col], distance
											found = True
											break
							if found:
								break
						if found:
							break
						distance += 1
				giant_sites_time_series.append([rep, id, y, x, year, month, dhw_nn, depth_i])
	print (year)




out = open('giant_sites_time_series_future.csv','w')
out.write('rep,site_id,lat,lon,year,month,dhw_max,depth\n')
for i in giant_sites_time_series:
	o = out.write(','.join(map(str,i))+'\n')


out.close()

