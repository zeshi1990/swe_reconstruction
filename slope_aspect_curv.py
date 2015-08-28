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
    slope_name = "DEM/slope_" + str(res) + "m.npy"
    aspect_name = "DEM/aspect_" + str(res) + "m.npy"
    np.save(slope_name, slope)
    np.save(aspect_name, aspect)
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
    curv_name = "DEM/curvature_" + str(res) + "m.npy"
    np.save(curv_name, curv)
    return curv


class wind_redistribution():

    ALL_STATIONS_ID = np.array([ 720267.,  720614.,  720645.,  725845.,  725846.,  725847.,
        749486.,   42603.,   41901.,   41908.,   42608.,   41909.]).astype(int)

    def __init__(self, collection, res=500):
        self.slope = gdal.Open("DEM/slope_" + str(res) + "m.npy")
        self.aspect = gdal.Open("DEM/aspect_" + str(res) + "m.npy")
        self.curvature = gdal.Open("DEM/curvature_" + str(res) + "m.npy")
        self.res = res
        self.collection = collection

    def wind_ts_query(self, station_id, datetime):
        query_json = {"date": datetime, "site_id": station_id}
        result = self.collection.find(query_json)
        return result[0]["wind speed"], result[0]["wind direction"]

    def wind_redistribution(self, datetime):
        datetime_obj = datetime
        u_wind = 0.
        v_wind = 0.
        for station_id in self.ALL_STATIONS_ID:
            weight = np.load("wind/" + str(int(station_id)) + "_" + str(self.res) + "m_weight.npy")
            wind_speed, wind_dir = self.wind_ts_query(station_id, datetime_obj)
            wind_dir_rad = wind_dir / 180.0 * np.pi
            u = - wind_speed * np.cos(3.0 * np.pi / 2.0 - wind_dir_rad)
            v = - wind_speed * np.sin(3.0 * np.pi / 2.0 - wind_dir_rad)
            u_wind += u * weight
            v_wind += v * weight
        wind_magnitude = np.sqrt(u_wind ** 2 + v_wind ** 2)
        theta = 3 * np.pi / 2 - np.arctan2(v_wind, u_wind)
        omega_s = self.slope * np.cos(theta - self.aspect)
        omega_c = self.curvature
        wind_correction_weight =  1. + 0.5 * omega_s + 0.5 * omega_c
        return wind_correction_weight * wind_magnitude

def main():
    pass

if __name__ == "__main__":
    main()
