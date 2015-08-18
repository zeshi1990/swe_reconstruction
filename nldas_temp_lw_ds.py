__author__ = 'zeshi'

###################################################################################
#
# This function is for downscaling all the altitudinal dependent NLDAS altributes
# DLWRF, TMP, UGRD, VGRD, APCP
# Wind information may need to be corrected by using MicroMet model
# We found that TMP and DLWRF are valid for altitudinal interpolation
#
###################################################################################

import gdal, osr
import numpy as np
from osgeo.gdal import GA_ReadOnly
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression, RANSACRegressor
from skimage.transform import resize

# This function will find the band number of the variable that is of interest from the
# .GRB file
def find_band_number(dataset, variable):
    '''
    Finds the band number inside the GRIB file, given the variable and the level names
    '''
    for i in range(1,dataset.RasterCount + 1):
        band = dataset.GetRasterBand(i)
        metadata = band.GetMetadata()
        band_variable = metadata['GRIB_ELEMENT']
        if (variable == band_variable):
            return i
    return None


# This function proceed bilinear interpolation from nldas to 1/8 degree dem boundary
# since the cell size are the same, no need to consider the difference
def bilinear_interpolation(dem, attribute, dem_georef, attribute_georef):
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

    return bilinear_interpolation_result

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

# This is an example to validate the bilinear interpolation code could work
def example_bilinear_interpolation():
    dem_filename = "DEM/8dge_dem.tif"
    dem_dataset = gdal.Open(dem_filename, GA_ReadOnly)
    dem_raster = dem_dataset.ReadAsArray()
    dem_georef = dem_dataset.GetGeoTransform()

    nldas_filename = "NLDAS/NLDAS_FORA0125_H.A20150620.1400.002.grb"
    attribute_name = "TMP"
    nldas_dataset = gdal.Open(nldas_filename, GA_ReadOnly)
    nldas_georef = nldas_dataset.GetGeoTransform()
    tmp_band_id = find_band_number(nldas_dataset, attribute_name)
    band_tmp = nldas_dataset.GetRasterBand(tmp_band_id)
    tmp_raster = band_tmp.ReadAsArray()
    bilinear_result = bilinear_interpolation(dem_raster, tmp_raster, dem_georef, nldas_georef)
    lr, residual = regression_information(dem_raster, bilinear_result)
    new_dem_fn = "data/toporad/radiation/DEM/500m_dem.tif"
    new_dem_dataset = gdal.Open(new_dem_fn, GA_ReadOnly)
    new_dem = new_dem_dataset.ReadAsArray()
    upscale_raster = apply_regression_on_dem(lr, new_dem, residual)
    f, ax = plt.subplots(1, 2)
    im = ax[0].imshow(bilinear_result, vmin=np.min(upscale_raster), vmax=np.max(upscale_raster), interpolation='none')
    ax[1].imshow(upscale_raster, vmin=np.min(upscale_raster), vmax=np.max(upscale_raster))
    plt.colorbar(im)
    plt.show()

# This function creates a new raster data (.tif) that of different georeference
def array2raster(rasterfn,newRasterfn,array,new_y_start_idx, new_x_start_idx):
    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    [rows, cols] = array.shape
    originX = originX + pixelWidth * new_x_start_idx
    originY = originY + pixelHeight * new_y_start_idx

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt(raster.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

def main():
    example_bilinear_interpolation()

if __name__ == "__main__":
    main()


# dge_8 = "snsrt_8dge"
# m_30 = "snsrt_30m"
# m_500 = "snsrt_500m"
# dem_filename = "data/toporad/radiation/DEM/snsrt_8dge/hdr.adf"
# new_dem_filename = "data/toporad/radiation/DEM/snsrt_8dge/dem.tif"
# dem_dataset = gdal.Open(dem_filename, GA_ReadOnly)
# raster = dem_dataset.ReadAsArray()
# print raster
# [nodata_idx_y, nodata_idx_x] = np.where(raster < 0)
#
# for (i,j) in zip(nodata_idx_y, nodata_idx_x):
#     raster[i, j] = np.round((raster[i-1, j-1] + raster[i, j-1] + raster[i-1, j]) / 3)
# print raster
# array2raster(dem_filename, new_dem_filename, raster, 0, 0)
#

# for code in [dge_8, m_30, m_500]:
#     dem_filename = "data/toporad/radiation/DEM/" + code + "/hdr.adf"
#     dem_dataset = gdal.Open(dem_filename, GA_ReadOnly)
#     dem_raster = dem_dataset.ReadAsArray()
#     # print dem_raster
#     [idx_y, idx_x] = np.where(dem_raster < 0)
#     count_idx_y = np.bincount(idx_y)
#     ii = np.nonzero(count_idx_y)[0]
#     print zip(ii, count_idx_y[ii])
#     count_idx_x = np.bincount(idx_x)
#     jj = np.nonzero(count_idx_x)[0]
#     print zip(jj, count_idx_x[jj])

# dem_30m_filename = "data/toporad/radiation/DEM/snsrt_30m/hdr.adf"
# dem_dataset = gdal.Open(dem_30m_filename, GA_ReadOnly)
# dem_raster = dem_dataset.ReadAsArray()
# dem_raster = dem_raster[26:2630, 49:4931]
# new_dem_30m_filename = "data/toporad/radiation/DEM/snsrt_30m/dem.tif"
# array2raster(dem_30m_filename, new_dem_30m_filename, dem_raster, 26, 49)
#
#
#
# dem_500m_filename = "data/toporad/radiation/DEM/snsrt_500m/hdr.adf"
# dem_dataset = gdal.Open(dem_500m_filename, GA_ReadOnly)
# dem_raster = dem_dataset.ReadAsArray()
# dem_raster = dem_raster[2:158, 3:]
# new_dem_500m_filename = "data/toporad/radiation/DEM/snsrt_500m/dem.tif"
# array2raster(dem_500m_filename, new_dem_500m_filename, dem_raster, 2, 3)
# dem_30m_filename = "data/toporad/radiation/DEM/snsrt_30m/hdr.adf"
# new_dem_30m_filename = "data/toporad/radiation/DEM/snsrt_30m/dem.tif"
# dem_old = gdal.Open(dem_30m_filename)
# dem_new = gdal.Open(new_dem_30m_filename)
# print dem_old.GetGeoTransform()
# print dem_new.GetGeoTransform()
# band = dem_new.GetRasterBand(1)
# print band.XSize
# print band.YSize
# print band.ReadAsArray().shape

    # originally under def apply_regression_on_dem
    # Code below might be valuable
    # new_dem_xRef = new_dem_georef[0]
    # new_dem_yRef = new_dem_georef[3]
    # new_dem_xCell = new_dem_georef[1]
    # new_dem_yCell = new_dem_georef[5]
    #
    # new_dem_idx_row, new_dem_idx_col = np.indices(new_dem.shape)
    # new_dem_idx_row = new_dem_idx_row.flatten()
    # new_dem_idx_col = new_dem_idx_col.flatten()
    #
    # new_dem_pixel_coord_x = new_dem_idx_col * new_dem_xCell + new_dem_xRef + 0.5 * new_dem_xCell
    # new_dem_pixel_coord_y = new_dem_idx_row * new_dem_yCell + new_dem_yRef + 0.5 * new_dem_yCell
    #
    # old_dem_idx_row = np.floor((new_dem_pixel_coord_y - old_dem_yRef) / old_dem_yCell).astype(int)
    # old_dem_idx_col = np.floor((new_dem_pixel_coord_x - old_dem_xRef) / old_dem_xCell).astype(int)
    #
    # old_dem_pixel_coord_x = old_dem_idx_col * old_dem_xCell + old_dem_xRef + 0.5 * old_dem_xCell
    # old_dem_pixel_coord_y = old_dem_idx_row * old_dem_yCell + old_dem_yRef + 0.5 * old_dem_yCell
    #
    # old_new_dem_coord_x_diff = new_dem_pixel_coord_x - old_dem_pixel_coord_x
    # old_new_dem_coord_y_diff = new_dem_pixel_coord_y - old_dem_pixel_coord_y
    #
    # x_diff_g_zero = old_new_dem_coord_x_diff > 0
    # x_diff_leq_zero = old_new_dem_coord_x_diff <= 0
    # y_diff_g_zero = old_new_dem_coord_y_diff > 0
    # y_diff_leq_zero = old_new_dem_coord_y_diff <= 0
    #
    # bottom_left_idx_col = np.zeros(len(new_dem_idx_row))
    # bottom_left_idx_row = np.zeros(len(new_dem_idx_row))
    # bottom_right_idx_col = np.zeros(len(new_dem_idx_row))
    # bottom_right_idx_row = np.zeros(len(new_dem_idx_row))
    # top_left_idx_col = np.zeros(len(new_dem_idx_row))
    # top_left_idx_row = np.zeros(len(new_dem_idx_row))
    # top_right_idx_col = np.zeros(len(new_dem_idx_row))
    # top_right_idx_row = np.zeros(len(new_dem_idx_row))
    #
    # # Need to address the four quadrants of the differences
    # # 1st quadrants
    # bottom_left_idx_col[np.where(np.logical_and(x_diff_g_zero, y_diff_g_zero))] = \
    #     old_dem_idx_col[np.where(np.logical_and(x_diff_g_zero, y_diff_g_zero))]
    # bottom_left_idx_row[np.where(np.logical_and(x_diff_g_zero, y_diff_g_zero))] = \
    #     old_dem_idx_row[np.where(np.logical_and(x_diff_g_zero, y_diff_g_zero))]
    # bottom_right_idx_col[np.where(np.logical_and(x_diff_g_zero, y_diff_g_zero))] = \
    #     old_dem_idx_col[np.where(np.logical_and(x_diff_g_zero, y_diff_g_zero))] + 1
    # bottom_right_idx_row[np.where(np.logical_and(x_diff_g_zero, y_diff_g_zero))] = \
    #     old_dem_idx_row[np.where(np.logical_and(x_diff_g_zero, y_diff_g_zero))]
    # top_left_idx_col[np.where(np.logical_and(x_diff_g_zero, y_diff_g_zero))] = \
    #     old_dem_idx_col[np.where(np.logical_and(x_diff_g_zero, y_diff_g_zero))]
    # top_left_idx_row[np.where(np.logical_and(x_diff_g_zero, y_diff_g_zero))] = \
    #     old_dem_idx_row[np.where(np.logical_and(x_diff_g_zero, y_diff_g_zero))] - 1
    # top_right_idx_col[np.where(np.logical_and(x_diff_g_zero, y_diff_g_zero))] = \
    #     old_dem_idx_col[np.where(np.logical_and(x_diff_g_zero, y_diff_g_zero))] + 1
    # top_right_idx_row[np.where(np.logical_and(x_diff_g_zero, y_diff_g_zero))] = \
    #     old_dem_idx_row[np.where(np.logical_and(x_diff_g_zero, y_diff_g_zero))] - 1
    #
    # # 2nd quadrants
    # bottom_left_idx_col[np.where(np.logical_and(x_diff_leq_zero, y_diff_g_zero))] = \
    #     old_dem_idx_col[np.where(np.logical_and(x_diff_leq_zero, y_diff_g_zero))] - 1
    # bottom_left_idx_row[np.where(np.logical_and(x_diff_leq_zero, y_diff_g_zero))] = \
    #     old_dem_idx_row[np.where(np.logical_and(x_diff_leq_zero, y_diff_g_zero))]
    # bottom_right_idx_col[np.where(np.logical_and(x_diff_leq_zero, y_diff_g_zero))] = \
    #     old_dem_idx_col[np.where(np.logical_and(x_diff_leq_zero, y_diff_g_zero))]
    # bottom_right_idx_row[np.where(np.logical_and(x_diff_leq_zero, y_diff_g_zero))] = \
    #     old_dem_idx_row[np.where(np.logical_and(x_diff_leq_zero, y_diff_g_zero))]
    # top_left_idx_col[np.where(np.logical_and(x_diff_leq_zero, y_diff_g_zero))] = \
    #     old_dem_idx_col[np.where(np.logical_and(x_diff_leq_zero, y_diff_g_zero))] - 1
    # top_left_idx_row[np.where(np.logical_and(x_diff_leq_zero, y_diff_g_zero))] = \
    #     old_dem_idx_row[np.where(np.logical_and(x_diff_leq_zero, y_diff_g_zero))] - 1
    # top_right_idx_col[np.where(np.logical_and(x_diff_leq_zero, y_diff_g_zero))] = \
    #     old_dem_idx_col[np.where(np.logical_and(x_diff_leq_zero, y_diff_g_zero))]
    # top_right_idx_row[np.where(np.logical_and(x_diff_leq_zero, y_diff_g_zero))] = \
    #     old_dem_idx_row[np.where(np.logical_and(x_diff_leq_zero, y_diff_g_zero))] - 1
    #
    # # 3rd quadrants
    # bottom_left_idx_col[np.where(np.logical_and(x_diff_leq_zero, y_diff_leq_zero))] = \
    #     old_dem_idx_col[np.where(np.logical_and(x_diff_leq_zero, y_diff_leq_zero))] - 1
    # bottom_left_idx_row[np.where(np.logical_and(x_diff_leq_zero, y_diff_leq_zero))] = \
    #     old_dem_idx_row[np.where(np.logical_and(x_diff_leq_zero, y_diff_leq_zero))] + 1
    # bottom_right_idx_col[np.where(np.logical_and(x_diff_leq_zero, y_diff_leq_zero))] = \
    #     old_dem_idx_col[np.where(np.logical_and(x_diff_leq_zero, y_diff_leq_zero))]
    # bottom_right_idx_row[np.where(np.logical_and(x_diff_leq_zero, y_diff_leq_zero))] = \
    #     old_dem_idx_row[np.where(np.logical_and(x_diff_leq_zero, y_diff_leq_zero))] + 1
    # top_left_idx_col[np.where(np.logical_and(x_diff_leq_zero, y_diff_leq_zero))] = \
    #     old_dem_idx_col[np.where(np.logical_and(x_diff_leq_zero, y_diff_leq_zero))] - 1
    # top_left_idx_row[np.where(np.logical_and(x_diff_leq_zero, y_diff_leq_zero))] = \
    #     old_dem_idx_row[np.where(np.logical_and(x_diff_leq_zero, y_diff_leq_zero))]
    # top_right_idx_col[np.where(np.logical_and(x_diff_leq_zero, y_diff_leq_zero))] = \
    #     old_dem_idx_col[np.where(np.logical_and(x_diff_leq_zero, y_diff_leq_zero))]
    # top_right_idx_row[np.where(np.logical_and(x_diff_leq_zero, y_diff_leq_zero))] = \
    #     old_dem_idx_row[np.where(np.logical_and(x_diff_leq_zero, y_diff_leq_zero))]
    #
    # # 4th quadrants
    # bottom_left_idx_col[np.where(np.logical_and(x_diff_g_zero, y_diff_leq_zero))] = \
    #     old_dem_idx_col[np.where(np.logical_and(x_diff_g_zero, y_diff_leq_zero))]
    # bottom_left_idx_row[np.where(np.logical_and(x_diff_g_zero, y_diff_leq_zero))] = \
    #     old_dem_idx_row[np.where(np.logical_and(x_diff_g_zero, y_diff_leq_zero))] + 1
    # bottom_right_idx_col[np.where(np.logical_and(x_diff_g_zero, y_diff_leq_zero))] = \
    #     old_dem_idx_col[np.where(np.logical_and(x_diff_g_zero, y_diff_leq_zero))] + 1
    # bottom_right_idx_row[np.where(np.logical_and(x_diff_g_zero, y_diff_leq_zero))] = \
    #     old_dem_idx_row[np.where(np.logical_and(x_diff_g_zero, y_diff_leq_zero))] + 1
    # top_left_idx_col[np.where(np.logical_and(x_diff_g_zero, y_diff_leq_zero))] = \
    #     old_dem_idx_col[np.where(np.logical_and(x_diff_g_zero, y_diff_leq_zero))]
    # top_left_idx_row[np.where(np.logical_and(x_diff_g_zero, y_diff_leq_zero))] = \
    #     old_dem_idx_row[np.where(np.logical_and(x_diff_g_zero, y_diff_leq_zero))]
    # top_right_idx_col[np.where(np.logical_and(x_diff_g_zero, y_diff_leq_zero))] = \
    #     old_dem_idx_col[np.where(np.logical_and(x_diff_g_zero, y_diff_leq_zero))] + 1
    # top_right_idx_row[np.where(np.logical_and(x_diff_g_zero, y_diff_leq_zero))] = \
    #     old_dem_idx_row[np.where(np.logical_and(x_diff_g_zero, y_diff_leq_zero))]
    #
    # # Need to address the edge pixel problems of the indices
    # bottom_left_idx_col[np.where(bottom_left_idx_col < 0)] = 0
    # bottom_left_idx_col[np.where(bottom_left_idx_col >= old_dem_xSize)] = old_dem_xSize
    # bottom_left_idx_row[np.where(bottom_left_idx_row < 0)] = 0
    # bottom_left_idx_row[np.where(bottom_left_idx_row >= old_dem_ySize)] = old_dem_ySize
    # bottom_right_idx_col[np.where(bottom_right_idx_col < 0)] = 0
    # bottom_right_idx_col[np.where(bottom_right_idx_col >= old_dem_xSize)] = old_dem_xSize
    # bottom_right_idx_row =
    # bottom_right_idx_row =
    # top_left_idx_col =
    # top_left_idx_col =
    # top_left_idx_row =
    # top_left_idx_row =
    # top_right_idx_col =
    # top_right_idx_col =
    # top_right_idx_row =
    # top_right_idx_row =


    # idx_y_array = np.arange(new_dem_ySize)
    # idx_x_array = np.arange(new_dem_xSize)
    #
    # new_dem_xCord = idx_x_array * new_dem_xCell + new_dem_xRef + 0.5 * new_dem_xCell
    # new_dem_yCord = idx_y_array * new_dem_yCell + new_dem_yRef + 0.5 * new_dem_yCell
    # old_dem_idx_x = np.floor((new_dem_xCord - old_dem_xRef) / old_dem_xCell).astype(int)
    # old_dem_idx_y = np.floor((new_dem_yCord - old_dem_yRef) / old_dem_yCell).astype(int)
    # old_dem_xCord_left = old_dem_idx_x * old_dem_xCell + old_dem_xRef
    # old_dem_xCord_right = (old_dem_idx_x + 1) * old_dem_xCell + old_dem_xRef
    # old_dem_yCord_up = old_dem_idx_y * old_dem_yCell + old_dem_yRef
    # old_dem_yCord_low = (old_dem_idx_y + 1) * old_dem_yCell + old_dem_yRef
