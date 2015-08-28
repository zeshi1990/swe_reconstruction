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
from nldas_sw_ds import day_of_year, find_band_raster
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression, RANSACRegressor
from skimage.transform import resize

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
    assert attr=="TMP" or attr=="DLWRF", "The attribute %r is not supported in linear interpolation" % attr
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
    return bilinear_result, downscale_raster

# This is an example to validate the bilinear interpolation code could work
def example_bilinear_interpolation():
    bilinear_result, upscale_raster = nldas_attr_ds(2001, 6, 6, 18, attr="TMP", res=500)
    f, ax = plt.subplots(1, 2)
    im = ax[0].imshow(bilinear_result, vmin=np.min(upscale_raster), vmax=np.max(upscale_raster), interpolation='none')
    ax[1].imshow(upscale_raster, vmin=np.min(upscale_raster), vmax=np.max(upscale_raster))
    plt.colorbar(im)
    plt.show()

def main():
    linear_check(2001, 1, 1, 18, attr="PRES")
if __name__ == "__main__":
    main()