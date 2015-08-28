__author__ = 'zeshi'
import gdal, osr
from gdal import GA_ReadOnly, GDT_Float32, GRA_Bilinear

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

def reproject_clip(src_fn, match_fn, dst_fn):
    # Source
    src_filename = src_fn
    src = gdal.Open(src_filename, GA_ReadOnly)
    src_proj = src.GetProjection()
    src_geotrans = src.GetGeoTransform()

    # We want a section of source that matches this:
    match_filename = match_fn
    match_ds = gdal.Open(match_filename, GA_ReadOnly)
    match_proj = match_ds.GetProjection()
    match_geotrans = match_ds.GetGeoTransform()
    wide = match_ds.RasterXSize
    high = match_ds.RasterYSize

    # Output / destination
    dst_filename = dst_fn
    dst = gdal.GetDriverByName('GTiff').Create(dst_filename, wide, high, 1, GDT_Float32)
    dst.SetGeoTransform( match_geotrans )
    dst.SetProjection( match_proj)

    # Do the work
    gdal.ReprojectImage(src, dst, src_proj, match_proj, GRA_Bilinear)