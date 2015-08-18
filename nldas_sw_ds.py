__author__ = 'zeshi'

###################################################################################
#
# This function is for downscaling the NLDAS Solar Radiation to toporad resolution
# Note: Still need to average daily solar radiation if possible
#
###################################################################################

import gdal
import numpy as np
from osgeo.gdal import GA_ReadOnly
from datetime import date
import pickle

# This function will find the band number of the variable that is of interest from the
# .GRB file
def find_band_raster(dataset, variable):
    '''
    Finds the band number inside the GRIB file, given the variable and the level names
    '''
    for i in range(1,dataset.RasterCount + 1):
        band = dataset.GetRasterBand(i)
        metadata = band.GetMetadata()
        band_variable = metadata['GRIB_ELEMENT']
        if (variable == band_variable):
            return dataset.GetRasterBand(i).ReadAsArray()
    return None

def day_of_year(year, month, day):
    Date = date(year, month, day)
    day_of_year = str(Date.timetuple().tm_yday).zfill(3)
    return str(day_of_year).zfill(3)

def smoothing():
    pass

def nldas_solar_ds(year, month, day, hour, res=500):
    # Load NLDAS solar radiation
    Date = date(year, month, day)
    day_of_year = str(Date.timetuple().tm_yday).zfill(3)
    filename = "NLDAS_data/" + str(year) + "/" + day_of_year + "/" + \
               "NLDAS_FORA0125_H.A" + str(year) + str(month).zfill(2) + \
               str(day).zfill(2) + "." + str(hour).zfill(2) + "00.002.grb"
    nldas_dataset = gdal.Open(filename, GA_ReadOnly)
    nldas_sw_raster= find_band_raster(nldas_dataset, "DSWRF")
    # scaler is the toporad weighting matrix calculated before
    # idx_array is the index array of the data that is going to be queried in the original NLDAS data
    scaler_fn = "toporad/" + str(res) + "m_" + str(month).zfill(2) + "_scale.npy"
    idx_array_fn = "toporad/" + str(res) + "m_nldas_idx.p"
    scaler = np.load(scaler_fn)
    idx_array = pickle.load(open(idx_array_fn, "rb"))
    nldas_sw_crop_raster = nldas_sw_raster[idx_array]
    nldas_sw_crop_scaled_raster = scaler * nldas_sw_crop_raster
    return nldas_sw_crop_scaled_raster


def main():
    nldas_solar_ds(2013, 4, 6, 18)
    return None

if __name__ == "__main__":
    main()