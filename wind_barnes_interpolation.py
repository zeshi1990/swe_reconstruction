__author__ = 'zeshi'
import numpy as np
from matplotlib import pyplot as plt
import pickle
from sklearn.neighbors import NearestNeighbors
from numpy.core.defchararray import rstrip
import gdal
from gdal import GA_ReadOnly
from datetime import datetime


##################################################################################
# This function implements the wind interpolation and correction
# function from NLDAS, ASOS, RAWS data. Interpolation is implemented
# from Barnes Objective Interpolation Function, the correction is implemented
# by the Liston et al. 2006 MicroMet model
##################################################################################
DELTA_N = 0.21215704869671478
ALL_STATIONS_ID = np.array([ 720267.,  720614.,  720645.,  725845.,  725846.,  725847.,
        749486.,   42603.,   41901.,   41908.,   42608.,   41909.])
ALL_STATIONS_LOC = np.array([[  38.955  , -121.082  ],
       [  38.909  , -121.351  ],
       [  38.7205 , -120.7515 ],
       [  39.277  , -120.71   ],
       [  39.32   , -120.136  ],
       [  38.8905 , -119.9975 ],
       [  39.224  , -121.003  ],
       [  38.90556, -120.69722],
       [  39.14389, -120.50889],
       [  39.09139, -120.73167],
       [  39.07167, -120.42167],
       [  39.08361, -120.17111]])


# Find the average distance of observations, delta_n
def cal_delta_n():
    from pymongo import MongoClient
    client = MongoClient()
    db = client.amr_proj
    collection = db.wind
    query_list = collection.find({"date": datetime(2000, 1, 1, 0, 0, 0)})
    station_id = np.zeros(12)
    station_loc = np.zeros((12,2))
    for i, item in enumerate(query_list):
        station_id[i] = item['site_id']
        station_loc[i] = [item['lat'], item['lon']]
    nbrs = NearestNeighbors(n_neighbors=2, algorithm='auto').fit(station_loc)
    distances, indices = nbrs.kneighbors(station_loc)
    delta_n = np.average(distances[:, 1])
    return delta_n, station_id, station_loc

# Start implementing the Barnes Objective Function
def barnes_obj_interpolate(res=500):
    dem_fn = "DEM/" + str(res) + "m_dem.tif"
    dem_ds = gdal.Open(dem_fn, GA_ReadOnly)
    dem_georef = dem_ds.GetGeoTransform()
    dem_raster = dem_ds.ReadAsArray()
    [row_idx, col_idx] = np.indices(dem_raster.shape)
    pixel_lat = (row_idx + 0.5) * dem_georef[5] + dem_georef[3]
    pixel_lon = (col_idx + 0.5) * dem_georef[1] + dem_georef[0]

    weight_list = []
    total_weight = 0

    for station_loc in ALL_STATIONS_LOC:
        print station_loc
        lat_diff = pixel_lat - station_loc[0]
        lon_diff = pixel_lon - station_loc[1]
        distance = np.square(lat_diff) + np.square(lon_diff)
        k = 5.052 * (2 * DELTA_N / np.pi) ** 2
        weight = np.exp(-distance / k)
        weight_list.append(weight)

    for weight in weight_list:
        total_weight += weight

    for i, weight in enumerate(weight_list):
        divided_weight = weight / total_weight
        station_id = str(int(ALL_STATIONS_ID[i]))
        print station_id, ALL_STATIONS_LOC[i]
        file_name = "wind/" + station_id + "_" + str(res) + "m_weight.npy"
        np.save(file_name, divided_weight)

def main():
    pass

if __name__ == "__main__":
    main()
