__author__ = 'zeshi'

import gdal
from datetime import datetime, timedelta
from nldas_sw_ds import short_wave_ds
from nldas_temp_lw_ds import nldas_attr_ds
from nldas_wind_ds import wind_redistribution
from snow_surface_temp import snow_surface_temp
from snow_17 import snow17, calc_rh
from pymongo import MongoClient
from multiprocessing import Pool, cpu_count
from functools import partial
import time

class nldas():
    L = 334000                          # J/kg, Latent heat of fusion of water at 0 degree C
    RHO_WATER = 1000                    # kg/m^3, density of water
    DELTA_T = 3600                      # Second, Time step, 1 hour = 3600 second
    """
    This is the class for doing SWE estimation of each time step
    """
    def __init__(self, year, month, day, hour, shortwave_ds, wind_redistribution, res=500):
        self.datetime = datetime(year, month, day, 0, 0, 0)
        previous_step_datetime = self.datetime - timedelta(hours=1)
        self.DSW = shortwave_ds.nldas_solar_ds(year, month, day, hour, res)
        self.DLW = nldas_attr_ds(year, month, day, hour, "DLWRF", res)
        self.AIR_TEMP = nldas_attr_ds(year, month, day, hour, "TMP", res)
        self.PREV_STEP_AIR_TEMP = nldas_attr_ds(previous_step_datetime.year, previous_step_datetime.month,
                                                previous_step_datetime.day, previous_step_datetime.hour, "TMP", res)
        self.SPFH =nldas_attr_ds(year, month, day, hour, "SPFH", res)
        self.PRES = nldas_attr_ds(year, month, day, hour, "PRES", res)
        self.RH = calc_rh(self.AIR_TEMP, self.PRES, self.SPFH)
        self.SNOW_TEMP = snow_surface_temp(self.AIR_TEMP, self.PREV_STEP_AIR_TEMP, self.RH)
        self.U_WIND, self.V_WIND = nldas_attr_ds(year, month, day, hour, "WIND", res)
        self.WIND = wind_redistribution.wind_redistribution(self.U_WIND, self.V_WIND)

        sca_fn = "SCA/" + str(year) + "/sca" + str(year) + str(month).zfill(2) + str(day).zfill(2) + ".tif"
        abo_fn = "SCA/" + str(year) + "/abo" + str(year) + str(month).zfill(2) + str(day).zfill(2) + ".tif"
        if res == 500:
            self.FSCA = gdal.Open(sca_fn).ReadAsArray()
            self.ABO = gdal.Open(abo_fn).ReadAsArray()
        else:
            pass
        self.snow17 = snow17(self.AIR_TEMP, self.SNOW_TEMP, self.PRES, self.SPFH, self.WIND, self.DLW, res=res)
        self.LH = self.snow17.LH
        self.SH = self.snow17.SH
        self.CDLW = self.snow17.CDLW
        self.ULW = self.snow17.ULW




def nldas_wrapper(datetime, shortwave_ds, wind_redistribution, res=500):
    year = datetime.year
    month = datetime.day
    day = datetime.day
    hour = datetime.hour
    nldas(year, month, day, hour, shortwave_ds, wind_redistribution, res=res)

def parallel_test():
    wind_redist_obj = wind_redistribution(res=30)
    wind = MongoClient().amr_proj.wind
    sw_ds = short_wave_ds(res=30)
    datetime_list = []
    for hour in range(1, 19):
        datetime_list.append(datetime(2001, 1, 1, hour, 0, 0))
    partial_nldas_wrapper = partial(nldas_wrapper, shortwave_ds=sw_ds,
                                    wind_redistribution=wind_redist_obj, res=30)
    pool = Pool(processes=cpu_count())
    t = time.time()
    pool.map(partial_nldas_wrapper, datetime_list)
    print "Parallel computation time is:", time.time() - t

def test():
    t = time.time()
    wind_redist_obj = wind_redistribution(res=30)
    wind = MongoClient().amr_proj.wind
    sw_ds = short_wave_ds(res=30)
    print "Algorithm initialization cost :", time.time() - t, " sec"
    t = time.time()
    test_forcing = nldas(2001, 1, 1, 18, sw_ds, wind, wind_redist_obj, res=30)
    print "Algorithm initialization cost :", time.time() - t, " sec"

def main():
    parallel_test()

if __name__ == "__main__":
    main()