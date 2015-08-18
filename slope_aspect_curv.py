__author__ = 'zeshi'
import gdal
from gdal import GA_ReadOnly
import numpy as np
from scipy.signal import convolve2d

def cal_slope_aspect(res = 500):
    dem_fn = "DEM/" + str(res) + "m_dem.tif"
    dem_ds = gdal.Open(dem_fn, GA_ReadOnly)
    dem_array = dem_ds.ReadAsArray()
    # dem_array = np.zeros((10, 10))
    for i in range(0, len(dem_array)):
        dem_array[i] = 10 - i
    x, y = np.gradient(dem_array)
    slope = np.arctan(np.sqrt(x * x + y * y))
    aspect = 3 * np.pi / 2 - np.arctan2(-x, y)
    return slope, aspect

def cal_curvature(res = 500, scale = 500):
    dem_fn = "DEM/" + str(res) + "m_dem.tif"
    dem_ds = gdal.Open(dem_fn, GA_ReadOnly)
    dem_array = dem_ds.ReadAsArray()
    grid_num = 2 * int(float(scale) / float(res)) + 1
    matrix = np.zeros((grid_num, grid_num))
    matrix[grid_num/2, grid_num/2] = float(1) / float(scale) + float(1) / (np.sqrt(2) * float(scale))
    matrix[0, grid_num/2] = float(-1) / (4. * float(scale))
    matrix[grid_num/2, 0] = float(-1) / (4. * float(scale))
    matrix[grid_num/2, -1] = float(-1) / (4. * float(scale))
    matrix[-1, grid_num/2] = float(-1) / (4. * float(scale))
    matrix[0, 0] = float(-1) / (4. * np.sqrt(2) * float(scale))
    matrix[0, -1] = float(-1) / (4. * np.sqrt(2) * float(scale))
    matrix[-1, 0] = float(-1) / (4. * np.sqrt(2) * float(scale))
    matrix[-1, -1] = float(-1) / (4. * np.sqrt(2) * float(scale))
    matrix *= 1.0 / 4
    curv = convolve2d(dem_array, matrix, mode='full')
    abs_max, abs_min = abs(np.max(curv)), abs(np.min(curv))
    max_scaler = np.max([abs_max, abs_min])
    curv *= 0.5/max_scaler
    return curv

def wind_ts_query(station_id, datetime, db):
    collection = db.wind_ts
    query_json = {"time": datetime, "station_id": station_id}
    return collection.get



def wind_redistribution(datetime, db, res=500, scale=500):
    from wind_barnes_interpolation import barnes_obj_interpolate
    wind_station_array = [0,0,0,0,0]
    u_wind = 0
    v_wind = 0
    for station_id, weight in barnes_obj_interpolate(res=res):
        u_wind_temp, v_wind_temp = wind_ts_query(station_id, datetime, db) #pseudo code for querying the wind from time series
        u_wind += u_wind_temp * weight
        v_wind += v_wind_temp * weight
    wind_magnitude = np.sqrt(u_wind ** 2 + v_wind ** 2)
    theta = 3 * np.pi / 2 - np.arctan2(v_wind, u_wind)
    slope, aspect = cal_slope_aspect(res=res)
    omega_s = slope * np.cos(theta - aspect)
    omega_c = cal_curvature(res=res, scale=scale)
    wind_correction_weight =  1 + 0.5 * omega_s + 0.5 * omega_c
    theta_t = theta - 0.5 * omega_s * np.sin(2*(aspect - theta))
    return wind_correction_weight * wind_magnitude, theta_t




def main():
    slope, aspect = cal_slope_aspect()
    curv = cal_curvature()
    from nldas_temp_lw_ds import array2raster

if __name__ == "__main__":
    main()
