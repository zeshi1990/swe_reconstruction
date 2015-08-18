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
from matplotlib import pyplot as plt
from nldas_temp_lw_ds import array2raster
import pickle

# Finds the band number inside the GRIB file, given the variable and the level names
def find_band_number(dataset, variable):
    for i in range(1,dataset.RasterCount + 1):
        band = dataset.GetRasterBand(i)
        metadata = band.GetMetadata()
        band_variable = metadata['GRIB_ELEMENT']
        if (variable == band_variable):
            return i
    return None

# Returns raster data and georeference of toporad data
def read_toporad(month, res=500):
    filename=""
    if res == 500:
        if month <= 5:
            filename = "data/toporad/radiation/reproject_toporad/net_rad" + "050" + str(month) + "15/hdr.adf"
        elif month == 11 or month == 12:
            filename = "data/toporad/radiation/reproject_toporad/net_rad" + "04" + str(month) + "15/hdr.adf"
        else:
            raise NameError("Month should be from 11 - 5")
    elif res == 30:
        if month <= 5:
            filename = "data/toporad/radiation/reproject_toporad/rad" + "050" + str(month) + "15/hdr.adf"
        elif month == 11 or month == 12:
            filename = "data/toporad/radiation/reproject_toporad/rad" + "04" + str(month) + "15/hdr.adf"
        else:
            raise NameError("Month should be from 11 - 5")
    dataset = gdal.Open(filename, GA_ReadOnly)
    dataset_raster = dataset.ReadAsArray()
    dataset_georef = dataset.GetGeoTransform()
    return dataset_raster, dataset_georef, filename

def nldas_ds(year, month, day, hour, res=500):
    # Load NLDAS solar radiation
    filename = "data/NLDAS/NLDAS_FORA0125_H.A" + str(year) + str(month).zfill(2) + str(day).zfill(2) + "." + str(hour).zfill(2) + "00.002.grb"
    nldas_dataset = gdal.Open(filename, GA_ReadOnly)
    nldas_sw_band_id = find_band_number(nldas_dataset, "DSWRF")
    nldas_sw_raster = nldas_dataset.GetRasterBand(nldas_sw_band_id).ReadAsArray()
    scaler_fn = "data/toporad/" + str(res) + "m_" + str(month).zfill(2) + "_scale.npy"
    idx_array_fn = "data/toporad/" + str(res) + "m_nldas_idx.p"
    scaler = np.load(scaler_fn)
    idx_array = pickle.load(open(idx_array_fn, "rb"))
    nldas_sw_crop_raster = nldas_sw_raster[idx_array]
    nldas_sw_crop_scaled_raster = scaler * nldas_sw_crop_raster
    f, ax = plt.subplots(3, 1)
    idx_array_min_row = np.min(idx_array[0])
    idx_array_max_row = np.max(idx_array[0])
    idx_array_min_col = np.min(idx_array[1])
    idx_array_max_col = np.max(idx_array[1])
    row = [idx_array_min_row, idx_array_max_row]
    print row
    col = [idx_array_min_col, idx_array_max_col]
    print col
    ax[0].imshow(nldas_sw_raster, vmin=np.min(nldas_sw_crop_scaled_raster),
                      vmax=np.max(nldas_sw_crop_scaled_raster))
    ax[0].plot(col, row, ".r")
    im = ax[1].imshow(nldas_sw_crop_raster, vmin=np.min(nldas_sw_crop_scaled_raster),
                      vmax=np.max(nldas_sw_crop_scaled_raster))
    im = ax[2].imshow(nldas_sw_crop_scaled_raster, vmin=np.min(nldas_sw_crop_scaled_raster),
                      vmax=np.max(nldas_sw_crop_scaled_raster))
    plt.colorbar(im)
    plt.show()


def main():
    nldas_ds(2013, 4, 6, 18)
    return None

if __name__ == "__main__":
    main()


'''
    # Not efficient version BELOW!!!
    # for idx_y in range(toporad_ySize):
    #     for idx_x in range(toporad_xSize):
    #         if toporad_raster[idx_y, idx_x] > 0:
    #             toporad_xCord = idx_x * toporad_xCell + toporad_xRef + 0.5 * toporad_xCell
    #             toporad_yCord = idx_y * toporad_yCell + toporad_yRef + 0.5 * toporad_yCell
    #             nldas_idx_x = int(np.floor((toporad_xCord - nldas_xRef) / nldas_xCell))   # need to ground to integer
    #             nldas_idx_y = int(np.floor((toporad_yCord - nldas_yRef) / nldas_yCell))   # need to ground to integer
    #             nldas_xCord_left = nldas_idx_x * nldas_xCell + nldas_xRef
    #             nldas_xCord_right = (nldas_idx_x + 1) * nldas_xCell + nldas_xRef
    #             nldas_yCord_up = nldas_idx_y * nldas_yCell + nldas_yRef
    #             nldas_yCord_low = (nldas_idx_y + 1) * nldas_yCell + nldas_yRef
    #             toporad_idx_xLeft = int(np.around((nldas_xCord_left - toporad_xRef) / toporad_xCell))
    #             toporad_idx_xRight = int(np.around((nldas_xCord_right - toporad_xRef) / toporad_xCell))
    #             toporad_idx_yUp = int(np.around((nldas_yCord_up - toporad_yRef) / toporad_yCell))
    #             toporad_idx_yLow = int(np.around((nldas_yCord_low - toporad_yRef) / toporad_yCell))
    #             toporad_raster_temp = toporad_raster[toporad_idx_yUp:toporad_idx_yLow, toporad_idx_xLeft:toporad_idx_xRight]
    #             toporad_raster_temp = toporad_raster_temp[np.where(toporad_raster_temp >= 0)]
    #             toporad_local = toporad_raster[idx_y, idx_x]
    #             nldas_ds_local = toporad_local / np.average(toporad_raster_temp) * nldas_sw_raster[nldas_idx_y, nldas_idx_x]
    '''

# for res in [30, 500]:
#     raster, georef = read_toporad(5, res=res)
#     # print dem_raster
#     [idx_y, idx_x] = np.where(raster < 0)
#     count_idx_y = np.bincount(idx_y)
#     ii = np.nonzero(count_idx_y)[0]
#     print zip(ii, count_idx_y[ii])
#     count_idx_x = np.bincount(idx_x)
#     jj = np.nonzero(count_idx_x)[0]
#     print zip(jj, count_idx_x[jj])
# for month in [11, 12, 1, 2, 3, 4, 5]:
#     raster_30m, georef, filename_30m = read_toporad(month, res=30)
#     raster_30m = raster_30m[25:2584, 48:4875]
#     print np.where(raster_30m < 0)
#     new_toporad_filename_30m = "data/toporad/radiation/reproject_toporad/30m_" + str(month).zfill(2) + ".tif"
#     array2raster(filename_30m, new_toporad_filename_30m, raster_30m, 25, 48)
#
#     raster, georef, filename = read_toporad(month, res=500)
#     raster = raster[2:155, 3:]
#     print np.where(raster < 0)
#     new_toporad_filename = "data/toporad/radiation/reproject_toporad/500m_" + str(month).zfill(2) + ".tif"
#     array2raster(filename, new_toporad_filename, raster, 2, 3)

# new_toporad_filename = "data/toporad/radiation/reproject_toporad/500m_05.tif"
# ds = gdal.Open(new_toporad_filename)
# raster = ds.ReadAsArray()
# georef_topo = ds.GetGeoTransform()
# # print georef_topo
# print raster.shape
#
# new_toporad_filename = "data/toporad/radiation/DEM/500m_dem.tif"
# ds = gdal.Open(new_toporad_filename)
# raster = ds.ReadAsArray()
# georef_dem = ds.GetGeoTransform()
# # print georef_dem
# print raster.shape
#
# print (georef_topo[0] - georef_dem[0])/georef_dem[1]
# print (georef_topo[3] - georef_dem[3])/georef_dem[5]
#
# new_toporad_filename = "data/toporad/radiation/reproject_toporad/30m_05.tif"
# ds = gdal.Open(new_toporad_filename)
# raster = ds.ReadAsArray()
# georef_topo = ds.GetGeoTransform()
# # print georef_dem
# print raster.shape
#
# new_toporad_filename = "data/toporad/radiation/DEM/30m_dem.tif"
# ds = gdal.Open(new_toporad_filename)
# raster = ds.ReadAsArray()
# georef_dem = ds.GetGeoTransform()
# # print georef_dem
# print raster.shape
#
# print (georef_topo[0] - georef_dem[0])/georef_dem[1]
# print (georef_topo[3] - georef_dem[3])/georef_dem[5]


# for res in [30, 500]:
#     for month in [1]:
#         print res, month
#         toporad_fn = "data/toporad/radiation/reproject_toporad/" + str(res) + "m_" + str(month).zfill(2) + ".tif"
#         ds = gdal.Open(toporad_fn, GA_ReadOnly)
#         topo_georef = ds.GetGeoTransform()
#         topo_raster = ds.ReadAsArray()
#         topLeft_lon = topo_georef[0]
#         topLeft_lat = topo_georef[3]
#         width = topo_georef[1]
#         height = topo_georef[5]
#
#         nldas_fn = "data/NLDAS/NLDAS_FORA0125_H.A20140606.1400.002.grb"
#         nldas_dataset = gdal.Open(nldas_fn, GA_ReadOnly)
#         nldas_georef = nldas_dataset.GetGeoTransform()
#         sw_band = find_band_number(nldas_dataset, "DSWRF")
#         sw_raster = nldas_dataset.GetRasterBand(sw_band).ReadAsArray()
#         [nldas_topLeft_lon, nldas_width, rot_1, nldas_topLeft_lat, rot_2, nldas_height] = nldas_georef
#
#         [topo_idx_row, topo_idx_col] = np.indices(topo_raster.shape)
#         # print topo_idx_row
#         # print topo_idx_col
#         topo_lat = (topo_idx_row.astype(float) + 0.5) * height + topLeft_lat
#         topo_lon = (topo_idx_col.astype(float) + 0.5) * width + topLeft_lon
#         # print topo_lat
#         # print topo_lon
#         nldas_topo_idx_row = np.floor((topo_lat - nldas_topLeft_lat) / nldas_height).astype(int)
#         nldas_topo_idx_col = np.floor((topo_lon - nldas_topLeft_lon) / nldas_width).astype(int)
#         # print nldas_topo_idx_row
#         # print nldas_topo_idx_col
#         nldas_idx = [nldas_topo_idx_row, nldas_topo_idx_col]
#         print sw_raster[nldas_idx]
#         fn = "data/toporad/" + str(res) + "m_nldas_idx.p"
#         pickle.dump(nldas_idx, open(fn, "wb"))

        # sw_subset = sw_raster[nldas_topo_idx_row, nldas_topo_idx_col]
        # # print sw_subset
        # # print len(np.unique(nldas_topo_idx_row)) * len(np.unique(nldas_topo_idx_col))
        # # print len(np.unique(sw_subset))
        # sw_subset_unique = np.unique(sw_subset)
        # # print sw_subset_unique
        # toporad_scale = np.zeros(sw_subset.shape).astype(float)
        # for item in sw_subset_unique:
        #     indices = np.where(sw_subset == item)
        #     toporad_scale[indices] = topo_raster[indices].astype(float) / np.average(topo_raster[indices].astype(float))
        # savename = "data/toporad/" + str(res) + "m_" + str(month).zfill(2) + "_scale.npy"
        # np.save(savename, toporad_scale)
        # print toporad_scale
        # plt.imshow(topo_raster)
        # plt.colorbar()
        # plt.show()
        # plt.imshow(toporad_scale)
        # plt.show()

# Need to save toporad_scale as raster data and also need to save