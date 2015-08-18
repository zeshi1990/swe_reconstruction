__author__ = 'zeshi'
import numpy as np
from matplotlib import pyplot as plt
import pickle
from sklearn.neighbors import NearestNeighbors
from numpy.core.defchararray import rstrip
import gdal
from gdal import GA_ReadOnly


##################################################################################
# This function implements the wind interpolation and correction
# function from NLDAS, ASOS, RAWS data. Interpolation is implemented
# from Barnes Objective Interpolation Function, the correction is implemented
# by the Liston et al. 2006 MicroMet model
##################################################################################

# Find the average distance of observations, delta_n
def cal_delta_n():
    num_link = 0
    stations_info = np.load("data/ASOS_wind/station_info.npy")
    stations_id = stations_info[:, 0]
    stations_loc = stations_info[:, 1:]
    id_string = stations_id.astype(int).astype(str)
    # print id_string
    id_string_fp = np.char.rjust(id_string, 6)
    unique_id_fp = np.unique(id_string_fp)
    unique_id_idx = []
    for id_fp in unique_id_fp:
        idx = np.where(id_string_fp == id_fp)
        # print idx
        unique_id_idx.append(idx[0][0])
    unique_id_string = id_string[unique_id_idx]
    unique_loc = stations_loc[unique_id_idx, :]
    unique_id = unique_id_string.astype(float)
    unique_station = np.column_stack((unique_id, unique_loc))
    # print unique_loc
    raws_stations_info = np.load("data/ASOS_wind/RAWS_station_info.npy")
    all_stations_info = np.append(unique_station, raws_stations_info, axis=0)
    all_stations_loc = all_stations_info[:, 1:]
    # print all_stations_loc
    nbrs = NearestNeighbors(n_neighbors=2, algorithm='auto').fit(all_stations_loc)
    distances, indices = nbrs.kneighbors(all_stations_loc)
    # print indices
    # print distances
    delta_n = np.average(distances[:, 1])
    return delta_n, all_stations_info

# Start implementing the Barnes Objective Function
def barnes_obj_interpolate(res=500):
    dem_fn = "data/toporad/radiation/DEM/" + str(res) + "m_dem.tif"
    dem_ds = gdal.Open(dem_fn, GA_ReadOnly)
    dem_georef = dem_ds.GetGeoTransform()
    dem_raster = dem_ds.ReadAsArray()
    [row_idx, col_idx] = np.indices(dem_raster.shape)
    pixel_lat = (row_idx + 0.5) * dem_georef[5] + dem_georef[3]
    pixel_lon = (col_idx + 0.5) * dem_georef[1] + dem_georef[0]

    delta_n, all_stations_info = cal_delta_n()
    weight_list = []
    total_weight = 0
    print all_stations_info[:, 0]
    for station_info in all_stations_info:
        lat_diff = pixel_lat - station_info[1]
        lon_diff = pixel_lon - station_info[2]
        distance = np.square(lat_diff) + np.square(lon_diff)
        k = 5.052 * (2 * delta_n / np.pi) ** 2
        weight = np.exp(-distance / k)
        weight_list.append(weight)
    for weight in weight_list:
        total_weight += weight
    unit = 0
    for i, weight in enumerate(weight_list):
        divided_weight = weight / total_weight
        unit += divided_weight
        # print divided_weight
        station_id = str(int(all_stations_info[i, 0]))
        # print station_id
        file_name = station_id + "_weight.npy"
        print len(pickle.dumps(divided_weight, -1))
        yield station_id, divided_weight

def main():
    barnes_obj_interpolate()

if __name__ == "__main__":
    main()
