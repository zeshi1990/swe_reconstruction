__author__ = 'zeshi'

import gdal
from gdal import GA_ReadOnly
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime, timedelta, date
from nldas_sw_ds import short_wave_ds
from nldas_temp_lw_ds import nldas_attr_ds
from nldas_wind_ds import wind_redistribution
from snow_surface_temp import snow_surface_temp
from snow_17 import snow17, calc_rh
from geo_spatial_tool import array2raster
from multiprocessing import Pool, cpu_count
from functools import partial
import logging
import time
import gc

class delta_swe_hourly():
    L = 334000                          # J/kg, Latent heat of fusion of water at 0 degree C
    RHO_WATER = 1000                    # kg/m^3, density of water
    DELTA_T = 3600                      # Second, Time step, 1 hour = 3600 second
    """
    This is the class for doing SWE estimation of each time step
    """
    def __init__(self, year, month, day, hour, shortwave_ds, wind_redistribution, res=500):
        self.datetime = datetime(year, month, day, hour, 0, 0)
        previous_step_datetime = self.datetime - timedelta(hours=1)
        self.DSW = shortwave_ds.nldas_solar_ds(year, month, day, hour, res)
        self.DLW = nldas_attr_ds(year, month, day, hour, "DLWRF", res)
        self.AIR_TEMP = nldas_attr_ds(year, month, day, hour, "TMP", res)
        self.PREV_STEP_AIR_TEMP = nldas_attr_ds(previous_step_datetime.year, previous_step_datetime.month,
                                                previous_step_datetime.day, previous_step_datetime.hour, "TMP", res)
        self.SPFH =nldas_attr_ds(year, month, day, hour, "SPFH", res)
        self.SPFH[np.where(self.SPFH < 0)] = 0.001
        self.PRES = nldas_attr_ds(year, month, day, hour, "PRES", res)
        self.APCP = nldas_attr_ds(year, month, day, hour, "APCP", res)
        self.RH = calc_rh(self.AIR_TEMP, self.PRES, self.SPFH)
        self.SNOW_TEMP = snow_surface_temp(self.AIR_TEMP, self.PREV_STEP_AIR_TEMP, self.RH)
        self.U_WIND, self.V_WIND = nldas_attr_ds(year, month, day, hour, "WIND", res)
        self.WIND = wind_redistribution.wind_redistribution(self.U_WIND, self.V_WIND)

        self.abo_fn = "SCA/" + str(year) + "/abo" + str(year) + str(month).zfill(2) + str(day).zfill(2) + ".tif"
        self.veg = gdal.Open("vegetation/amr_nlcd2011_" + str(res) + "m.tif", GA_ReadOnly).ReadAsArray().astype(float)
        if res == 500:
            self.ABO = gdal.Open(self.abo_fn, GA_ReadOnly).ReadAsArray()
        else:
            pass

        self.snow17 = snow17(self.AIR_TEMP, self.SNOW_TEMP, self.PRES, self.SPFH, self.WIND, self.DLW, self.APCP, res=res)
        self.LH = self.snow17.LH
        self.SH = self.snow17.SH
        self.CDLW = self.snow17.CDLW
        self.ULW = self.snow17.ULW
        self.Qm = self.snow17.Qm
        self.PRECIP_SNOW = self.APCP
        self.PRECIP_SNOW[self.AIR_TEMP > 0.] = 0.

    def run(self, reset=False):
        if reset == False:
            delta_swe = self.DSW * (1 - self.ABO) + self.CDLW - self.ULW + self.LH + self.SH
        else:
            delta_swe = self.DSW * (1 - self.ABO) + self.CDLW - self.ULW + self.LH + self.SH
            delta_swe[delta_swe < 0.] = 0.
        return delta_swe, self.PRECIP_SNOW

    def visual_net_rad(self):
        sw_net = self.DSW * (1 - self.ABO)
        lw_net = self.CDLW - self.ULW
        abo = gdal.Open('SCA/2006/abo20060524.tif').ReadAsArray()
        dsw = self.DSW
        dsw[abo < 0.5] = float("nan")
        sw_net[abo < 0.5] = float("nan")
        lw_net[abo < 0.5] = float("nan")
        plt.imshow(dsw)
        plt.colorbar()
        plt.show()
        # plt.imshow(lw_net)
        # plt.colorbar()
        # plt.show()
        # plt.imshow(sw_net + lw_net)
        # plt.colorbar()
        # plt.show()

    def energy_return(self):
        return self.DSW * (1 - self.ABO), self.CDLW - self.ULW, self.LH + self.SH

def delta_swe_hourly_wrapper(datetime, shortwave_ds, wind_redistribution, res=500):
    try:
        year = datetime.year
        month = datetime.month
        day = datetime.day
        hour = datetime.hour
        print str(year) + "-" + str(month).zfill(2) + "-" + str(day).zfill(2) + "-" + str(hour).zfill(2), \
            " delta SWE computation starts"
        delta_swe = delta_swe_hourly(year, month, day, hour, shortwave_ds, wind_redistribution, res=res)
        energy, precip = delta_swe.run(reset=True)
        delta_swe = None
        return (energy, precip)
    except AttributeError, e:
        print e
        return None
    except RuntimeError, e:
        print e
        return None

class delta_swe_daily():
    def __init__(self, date, sw_ds, wind_rds, res=500):
        datetime_list = []
        for hour in range(0, 24):
            datetime_list.append(datetime(date.year, date.month, date.day, hour))
        part_hourly_energy = partial(delta_swe_hourly_wrapper,
                                     shortwave_ds = sw_ds,
                                     wind_redistribution = wind_rds,
                                     res=500)
        pool = Pool(processes=cpu_count())
        daily_energy_list = pool.map(part_hourly_energy, datetime_list)
        pool.close()
        pool.join()
        self.daily_energy = 0.
        self.daily_precip_snow = 0.
        for i in range(0, len(daily_energy_list)):
            self.daily_energy += daily_energy_list[i][0]
            self.daily_precip_snow += daily_energy_list[i][1] / 1000.
        self.daily_energy[self.daily_energy < 0.] = 0.
        sca_fn = "SCA/" + str(date.year) + "/sca" + str(date.year) + \
                 str(date.month).zfill(2) + str(date.day).zfill(2) + ".tif"
        self.FSCA = gdal.Open(sca_fn, GA_ReadOnly).ReadAsArray()
        self.FRA = gdal.Open("vegetation/amr_nlcd2011_500m.tif", GA_ReadOnly).ReadAsArray().astype(float)
        self.FRA[self.FRA >= 100.0] = 99.0
        self.FSCA = self.FSCA / (1.0 - self.FRA / 100.0)
        self.FSCA[self.FSCA > 1.0] = 1.0
        self.delta_swe = self.FSCA * self.daily_energy * 3600.0 / (343000.0 * 1000.0) - self.daily_precip_snow
        export_fn = "output/swe_change_reset_zero_hourly/" + str(date.year) + str(date.month).zfill(2) + \
                    str(date.day).zfill(2) + "delta_swe.tif"
        array2raster(sca_fn, export_fn, self.delta_swe, 0, 0)
        daily_energy_list = None
        self.daily_energy = None
        self.FSCA = None
        self.FRA = None
        self.delta_swe = None
        gc.collect()

def series_run(wy):
    wind_redist_obj = wind_redistribution(res=500)
    sw_ds = short_wave_ds(res=500)
    end_month_list = [6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 4]
    end_day_list = [30, 31, 31, 31, 31, 31, 31, 31, 31, 31, 30, 30]
    i = wy - 2001
    start_date = date(wy, 4, 1)
    end_date = date(wy, end_month_list[i], end_day_list[i])
    temp_date = start_date
    while temp_date <= end_date:
        swe = delta_swe_daily(temp_date, sw_ds, wind_redist_obj, res=500)
        temp_date += timedelta(days=1)
        swe = None
        gc.collect()

def sum_daily_delta_swe(wy):
    end_month_list = [6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 4]
    end_day_list = [30, 31, 31, 31, 31, 31, 31, 31, 31, 31, 30, 30]
    i = wy - 2001
    print "Start summing up year of", wy
    start_date = date(wy, 4, 1)
    end_date = date(wy, end_month_list[i], end_day_list[i])
    temp_date = end_date
    accumulate_swe = 0.
    while temp_date >= start_date:
        delta_swe_fn = "output/swe_change_reset_zero_hourly/" + str(temp_date.year) + str(temp_date.month).zfill(2) + \
                str(temp_date.day).zfill(2) + "delta_swe.tif"
        delta_swe = gdal.Open(delta_swe_fn, GA_ReadOnly).ReadAsArray()
        accumulate_swe += delta_swe
        export_fn = "output/swe_500m_reset_zero_hourly/" + str(wy) + "/" + str(wy) + str(temp_date.month).zfill(2) + \
                str(temp_date.day).zfill(2) + "_swe.tif"
        array2raster(delta_swe_fn, export_fn, accumulate_swe, 0., 0.)
        temp_date -= timedelta(days=1)

def test():
    for hour in range(12, 19):
        wind_redist_obj = wind_redistribution(res=500)
        sw_ds = short_wave_ds(res=500)
        test_forcing = delta_swe_hourly(2006, 5, 24, hour, sw_ds, wind_redist_obj, res=500)
        test_forcing.visual_net_rad()

def main():
    gdal.UseExceptions()
    for year in range(2001, 2013):
        series_run(year)
        sum_daily_delta_swe(year)
    # test()

if __name__ == "__main__":
    main()

# def parallel_running():
#     logging.basicConfig(filename="reconstructionError.log", level=logging.INFO)
#     logger = logging.getLogger()
#     wind_redist_obj = wind_redistribution(res=500)
#     sw_ds = short_wave_ds(res=500)
#     end_month_list = [6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 4]
#     end_day_list = [30, 31, 31, 31, 31, 31, 31, 31, 31, 31, 30, 30]
#     datetime_list = []
#     for i, wy in enumerate(range(2001, 2013)):
#         start_datetime = datetime(wy, 4, 1, hour = 0)
#         end_datetime = datetime(wy, end_month_list[i], end_day_list[i], hour=23)
#         temp_datetime = start_datetime
#         while temp_datetime <= end_datetime:
#             datetime_list.append(temp_datetime)
#             temp_datetime += timedelta(hours=1)
#     partial_delta_swe_hourly_wrapper = partial(delta_swe_hourly_wrapper,
#                                                shortwave_ds=sw_ds,
#                                                wind_redistribution=wind_redist_obj,
#                                                res=500)
#     pool = Pool(processes=cpu_count())
#     Error = pool.map(partial_delta_swe_hourly_wrapper, datetime_list)
#     for err in Error:
#         logger.info(err)

# def sum_delta_swe(end_datetime, res=500):
#     wy = end_datetime.year
#     end_month = end_datetime.month
#     end_day = end_datetime.day
#     if res == 500:
#         ref_fn = "SCA/2001/sca20010101.tif"
#     else:
#         ref_fn = ""
#     start_date = date(wy, 4, 1)
#     end_date = date(wy, end_month, end_day)
#     temp_date = end_date
#     hour = range(0, 24)
#     total_sum_swe = 0.
#     while temp_date >= start_date:
#         year = temp_date.year
#         month = temp_date.month
#         day = temp_date.day
#         print "Starting to sum over the swe change of " + str(year) + "-" + str(month) + "-" + str(day)
#         day_sum_swe = 0.
#         for hr in hour:
#             input_fn = "output/swe_change/" + str(year) + str(month).zfill(2) + str(day).zfill(2) + \
#                         str(hr).zfill(2) + "delta_swe.tif"
#             rasterArray = gdal.Open(input_fn, GA_ReadOnly).ReadAsArray()
#             rasterArray[rasterArray<0] = 0.
#             day_sum_swe += rasterArray
#         total_sum_swe += day_sum_swe
#         output_fn = "output/swe_500m/" + str(year) + "/" + str(year) + str(month).zfill(2) + \
#                     str(day).zfill(2) + "_swe.tif"
#         array2raster(ref_fn, output_fn, total_sum_swe, 0, 0)
#         temp_date -= timedelta(days=1)
#     return

# def parallel_sum_delta_swe():
#     end_month_list = [6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 4]
#     end_day_list = [30, 31, 31, 31, 31, 31, 31, 31, 31, 31, 30, 30]
#     datetime_list = []
#     for i, wy in enumerate(range(2001, 2013)):
#         datetime_list.append(datetime(wy, end_month_list[i], end_day_list[i], 23))
#     partial_sum_delta_swe = partial(sum_delta_swe, res=500)
#     pool = Pool(processes=len(datetime_list))
#     pool.map(partial_sum_delta_swe, datetime_list)
#
# def parallel_test():
#     wind_redist_obj = wind_redistribution(res=30)
#     sw_ds = short_wave_ds(res=30)
#     datetime_list = []
#     for hour in range(1, 19):
#         datetime_list.append(datetime(2001, 1, 1, hour, 0, 0))
#     partial_delta_swe_hourly_wrapper = partial(delta_swe_hourly_wrapper, shortwave_ds=sw_ds,
#                                     wind_redistribution=wind_redist_obj, res=30)
#     pool = Pool(processes=cpu_count())
#     t = time.time()
#     pool.map(partial_delta_swe_hourly_wrapper, datetime_list)
#     print "Parallel computation time is:", time.time() - t