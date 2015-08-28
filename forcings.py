__author__ = 'zeshi'

import numpy as np
import gdal
from gdal import GA_ReadOnly
from datetime import datetime
from nldas_sw_ds import nldas_solar_ds
from nldas_temp_lw_ds import nldas_attr_ds

class nldas_forcings():

    def __init__(self, year, month, day, hour, collection, wind_redistribution, res=500):
        self.datetime = datetime(year, month, day, 0, 0, 0)
        self.collection = collection
        self.DSW = nldas_solar_ds(year, month, day, hour, res)
        self.DLW = nldas_attr_ds(year, month, day, hour, "DLWRF", res)
        self.TMP = nldas_attr_ds(year, month, day, hour, "TMP", res)
        self.SPFH =nldas_attr_ds(year, month, day, hour, "SPFH", res)
        self.PRES = nldas_attr_ds(year, month, day, hour, "PRES", res)
        self.WIND = wind_redistribution.wind_redistribution(self.datetime)
        sca_fn = "SCA/" + str(year) + "/sca" + str(year) + str(month).zfill(2) + str(day).zfill(2) + ".tif"
        abo_fn = "SCA/" + str(year) + "/abo" + str(year) + str(month).zfill(2) + str(day).zfill(2) + ".tif"
        if res == 500:
            self.FSCA = gdal.Open(sca_fn).ReadAsArray()
            self.ABO = gdal.Open(abo_fn).ReadAsArray()
        else:
            pass