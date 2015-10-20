__author__ = 'zeshi'

###################################################################################
#
# This function is for downscaling all the altitudinal dependent NLDAS altributes
# DLWRF, TMP, UGRD, VGRD, APCP, SPFH, PRES
# Wind information may need to be corrected by using MicroMet model
# We found that TMP and DLWRF are valid for altitudinal interpolation
#
###################################################################################

import gdal
import numpy as np
from osgeo.gdal import GA_ReadOnly
from datetime import date
import pylab as plt
from sklearn.linear_model import LinearRegression, RANSACRegressor
from skimage.transform import resize

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

# This function proceed bilinear interpolation from nldas to 1/8 degree dem boundary
# since the cell size are the same, no need to consider the difference
def bilinear_interpolation(attribute, attribute_georef):
    dem_filename = "DEM/8dge_dem.tif"
    dem_dataset = gdal.Open(dem_filename, GA_ReadOnly)
    dem = dem_dataset.ReadAsArray()
    dem_georef = dem_dataset.GetGeoTransform()

    (dem_ySize, dem_xSize) = dem.shape
    dem_xRef = dem_georef[0]
    dem_yRef = dem_georef[3]
    dem_xCell = dem_georef[1]
    dem_yCell = dem_georef[5]

    # print dem_georef
    # print attribute_georef

    (attribute_ySize, attribute_xSize) = attribute.shape
    attribute_xRef = attribute_georef[0]
    attribute_yRef = attribute_georef[3]

    idx_upperLeftX = int(np.floor((dem_xRef - attribute_xRef) / dem_xCell))
    idx_upperLeftY = int(np.floor((dem_yRef - attribute_yRef) / dem_yCell))

    x_x1 = (dem_xRef - attribute_xRef) / dem_xCell - float(idx_upperLeftX)
    x2_x = 1.0 - x_x1
    y_y1 = (dem_yRef - attribute_yRef) / dem_yCell - float(idx_upperLeftY)
    y2_y = 1.0 - y_y1

    # print x_x1, x2_x, y_y1, y2_y
    attribute_12 = attribute[idx_upperLeftY:idx_upperLeftY+dem_ySize, idx_upperLeftX:idx_upperLeftX+dem_xSize]
    attribute_21 = attribute[idx_upperLeftY+1:idx_upperLeftY+dem_ySize+1, idx_upperLeftX+1:idx_upperLeftX+dem_xSize+1]
    attribute_11 = attribute[idx_upperLeftY+1:idx_upperLeftY+dem_ySize+1, idx_upperLeftX:idx_upperLeftX+dem_xSize]
    attribute_22 = attribute[idx_upperLeftY:idx_upperLeftY+dem_ySize, idx_upperLeftX+1:idx_upperLeftX+dem_xSize+1]

    # print x2_x * y2_y + x_x1 * y2_y + x2_x * y_y1 + x_x1 * y_y1

    bilinear_interpolation_result = attribute_11 * x2_x * y2_y + attribute_21 * x_x1 * y2_y + \
        attribute_12 * x2_x * y_y1 + attribute_22 * x_x1 * y_y1

    return dem, bilinear_interpolation_result

# This function does linear regression on 1/8 degree data and return the linear model as well as the residual in higher
# resolution dem shape
def regression_information(dem, bilinear_interpolation_results):
    dem_shape = dem.shape
    # print dem_shape
    dem = dem.flatten()
    bilinear_interpolation_results = bilinear_interpolation_results.flatten()
    alt_data = np.column_stack((dem, bilinear_interpolation_results))
    alt_data = alt_data[np.where(alt_data[:, 0] > 0)]
    RANSAC_lr = RANSACRegressor(LinearRegression())
    RANSAC_lr.fit(alt_data[:, 0:1], alt_data[:, 1])
    predict_result = RANSAC_lr.predict(alt_data[:, 0:1]).transpose()[0]
    # print predict_result
    # print predict_result.shape
    residual = bilinear_interpolation_results - predict_result
    residual = np.reshape(residual, dem_shape)
    return RANSAC_lr, residual

# This function check linear relationship of feature
def linear_check(year, month, day, hour, attr="TMP"):
    yday = day_of_year(year, month, day)
    nldas_fn = "NLDAS_data/" + str(year) + "/" + yday + "/" + \
               "NLDAS_FORA0125_H.A" + str(year) + str(month).zfill(2) + \
               str(day).zfill(2) + "." + str(hour).zfill(2) + "00.002.grb"
    nldas_ds = gdal.Open(nldas_fn, GA_ReadOnly)
    nldas_gt = nldas_ds.GetGeoTransform()
    attr_raster = find_band_raster(nldas_ds, attr)
    dem, bilinear_result = bilinear_interpolation(attr_raster, nldas_gt)
    plt.plot(dem.flatten(), bilinear_result.flatten(), '.b')
    plt.show()

# This function check linear relationship of wind
def linear_check_wind(year, month, day, hour):
    yday = day_of_year(year, month, day)
    nldas_fn = "NLDAS_data/" + str(year) + "/" + yday + "/" + \
               "NLDAS_FORA0125_H.A" + str(year) + str(month).zfill(2) + \
               str(day).zfill(2) + "." + str(hour).zfill(2) + "00.002.grb"
    nldas_ds = gdal.Open(nldas_fn, GA_ReadOnly)
    nldas_gt = nldas_ds.GetGeoTransform()
    u_raster = find_band_raster(nldas_ds, "UGRD")
    v_raster = find_band_raster(nldas_ds, "VGRD")
    wind_raster = np.sqrt(u_raster ** 2 + v_raster ** 2)
    dem, bilinear_result = bilinear_interpolation(wind_raster, nldas_gt)
    # plt.imshow(bilinear_result)
    plt.plot(dem.flatten(), bilinear_result.flatten(), '.b')
    plt.show()
# This function apply the linear regression model implemented in regression_information function
# to the DEM that is of interest (either 30m or 500m) and the residual is interpolated across the area
def apply_regression_on_dem(lr, new_dem, residual):
    new_dem_flat_row = new_dem.flatten()
    new_dem_flat_column = new_dem_flat_row[:, np.newaxis]
    lr_result = lr.predict(new_dem_flat_column).transpose()[0]
    lr_result = np.reshape(lr_result, new_dem.shape)
    interpolated_residual = resize(residual, new_dem.shape)
    result = lr_result + interpolated_residual
    return result

def nldas_attr_ds(year, month, day, hour, attr="TMP", res=500):
    assert attr=="TMP" or attr=="DLWRF" or attr=="SPFH" or attr=="PRES" or attr=="WIND" or attr=="APCP",\
        "The attribute %r is not supported in linear interpolation" % attr
    if attr=="WIND":
        yday = day_of_year(year, month, day)
        nldas_fn = "NLDAS_data/" + str(year) + "/" + yday + "/" + \
                   "NLDAS_FORA0125_H.A" + str(year) + str(month).zfill(2) + \
                   str(day).zfill(2) + "." + str(hour).zfill(2) + "00.002.grb"
        nldas_ds = gdal.Open(nldas_fn, GA_ReadOnly)
        nldas_gt = nldas_ds.GetGeoTransform()
        u_wind_raster = find_band_raster(nldas_ds, "UGRD")
        v_wind_raster = find_band_raster(nldas_ds, "VGRD")
        dem, u_wind_bilinear_result = bilinear_interpolation(u_wind_raster, nldas_gt)
        dem, v_wind_bilinear_result = bilinear_interpolation(v_wind_raster, nldas_gt)
        u_wind_lr, u_wind_residual = regression_information(dem, u_wind_bilinear_result)
        v_wind_lr, v_wind_residual = regression_information(dem, v_wind_bilinear_result)
        if res==500:
            new_dem_fn = "DEM/500m_dem.tif"
        else:
            new_dem_fn = "DEM/30m_dem.tif"
        new_dem_ds = gdal.Open(new_dem_fn, GA_ReadOnly)
        new_dem = new_dem_ds.ReadAsArray()
        u_wind_downscale_raster = apply_regression_on_dem(u_wind_lr, new_dem, u_wind_residual)
        v_wind_downscale_raster = apply_regression_on_dem(v_wind_lr, new_dem, v_wind_residual)
        # attr_imshow(u_wind_bilinear_result, u_wind_downscale_raster)
        # attr_imshow(v_wind_bilinear_result, v_wind_downscale_raster)
        return u_wind_downscale_raster, v_wind_downscale_raster

    elif attr=='APCP':
        yday = day_of_year(year, month, day)
        nldas_fn = "NLDAS_data/" + str(year) + "/" + yday + "/" + \
                   "NLDAS_FORA0125_H.A" + str(year) + str(month).zfill(2) + \
                   str(day).zfill(2) + "." + str(hour).zfill(2) + "00.002.grb"
        nldas_ds = gdal.Open(nldas_fn, GA_ReadOnly)
        nldas_gt = nldas_ds.GetGeoTransform()
        attr_raster = find_band_raster(nldas_ds, attr)
        dem, bilinear_result = bilinear_interpolation(attr_raster,nldas_gt)
        if res==500:
            new_dem_fn = "DEM/500m_dem.tif"
        else:
            new_dem_fn = "DEM/30m_dem.tif"
        new_dem_ds = gdal.Open(new_dem_fn, GA_ReadOnly)
        new_dem = new_dem_ds.ReadAsArray()
        downscale_raster = resize(bilinear_result, new_dem.shape)
        return downscale_raster

    else:
        yday = day_of_year(year, month, day)
        nldas_fn = "NLDAS_data/" + str(year) + "/" + yday + "/" + \
                   "NLDAS_FORA0125_H.A" + str(year) + str(month).zfill(2) + \
                   str(day).zfill(2) + "." + str(hour).zfill(2) + "00.002.grb"
        nldas_ds = gdal.Open(nldas_fn, GA_ReadOnly)
        nldas_gt = nldas_ds.GetGeoTransform()
        attr_raster = find_band_raster(nldas_ds, attr)
        dem, bilinear_result = bilinear_interpolation(attr_raster,nldas_gt)
        lr, residual = regression_information(dem, bilinear_result)
        if res==500:
            new_dem_fn = "DEM/500m_dem.tif"
        else:
            new_dem_fn = "DEM/30m_dem.tif"
        new_dem_ds = gdal.Open(new_dem_fn, GA_ReadOnly)
        new_dem = new_dem_ds.ReadAsArray()
        downscale_raster = apply_regression_on_dem(lr, new_dem, residual)
        # attr_imshow(bilinear_result, downscale_raster)
        return downscale_raster

def attr_imshow(original, scaled):
        min = np.min([np.min(original), np.min(scaled)])
        max = np.max([np.max(original), np.max(scaled)])
        plt.subplot(2, 1, 1)
        plt.imshow(original, vmin=min, vmax=max)
        plt.subplot(2, 1, 2)
        plt.imshow(scaled, vmin=min, vmax=max)
        plt.colorbar(orientation="horizontal")
        plt.show()

def main():
    for day in range(1, 32):
        linear_check(2001, 1, day, 18, "APCP")
if __name__ == "__main__":
    main()