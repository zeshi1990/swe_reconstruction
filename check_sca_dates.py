__author__ = 'zeshi'

import gdal
from osgeo.gdal import GA_ReadOnly, GDT_Float32, GRA_Bilinear
from calendar import monthrange
from datetime import date, timedelta
import pandas as pd

def check_sca_dates(pre = "sca"):
    dates = []
    for wy in range(2001, 2013):
            if wy == 2001:
                month_list = range(1, 10)
            elif wy == 2009:
                month_list = range(1,8) + range(10, 13)
            else:
                month_list = range(1, 13)
            for month in month_list:
                if month >= 10:
                    year = wy - 1
                else:
                    year = wy
                num_days = monthrange(year, month)[1]
                for day in range(1, num_days + 1):
                    chk_fn = "SCA/" + str(wy) + "/" + pre  +str(year) + \
                             str(month).zfill(2) + str(day).zfill(2) + ".tif"
                    try:
                        gdal.Open(chk_fn, GA_ReadOnly).ReadAsArray()
                    except AttributeError:
                        dates.append(date(year, month, day))
                        pass
    ts = pd.Series(dates)
    fn = "SCA/" + pre + "_noData_Dates.csv"
    ts.to_csv(fn)

import gdal
from osgeo.gdal import GA_ReadOnly, GDT_Float32, GRA_Bilinear
from calendar import monthrange
from datetime import date, timedelta
import pandas as pd
# from nldas_temp_lw_ds import array2raster
from matplotlib import pyplot as plt

def reproject_clip_sca_abo(src_fn, match_fn, dst_fn, src_proj_boolean = False):
    # Source
    src_filename = src_fn
    src = gdal.Open(src_filename, GA_ReadOnly)
    src_proj = src.GetProjection()
    src_geotrans = src.GetGeoTransform()

    abo_20120101 = gdal.Open("SCA/2012/abo2012_rs/abo20120101/hdr.adf")
    proj = abo_20120101.GetProjection()

    if src_proj_boolean == True:
        proj = src_proj

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
    gdal.ReprojectImage(src, dst, proj, match_proj, GRA_Bilinear)


def reproject_all(type = "sca"):
    if type == "sca":
        pre = "sca"
    else:
        pre = "abo"
    for wy in [2009]:
        if wy == 2001:
            month_list = range(1, 10)
        elif wy == 2009:
            month_list = range(5,8) + range(10, 13)
        else:
            month_list = range(1, 13)
        for month in month_list:
            if month >= 10:
                year = wy - 1
            else:
                year = wy
            num_days = monthrange(year, month)[1]
            for day in range(1, num_days + 1):
                src_fn = "SCA/" + str(wy) + "/" + pre + str(wy) + "_rs/" + pre + str(year) + \
                         str(month).zfill(2) + str(day).zfill(2) + "/hdr.adf"
                match_fn = "DEM/500m_dem.tif"
                dst_fn = "SCA/" + str(wy) + "/" + pre + str(year) + str(month).zfill(2) + \
                         str(day).zfill(2) + ".tif"
                try:
                    reproject_clip_sca_abo(src_fn, match_fn, dst_fn)
                except AttributeError:
                    print str(year) + "/" + str(month) + "/" + str(day) + " " + pre + " data not in file system!"
                    pass

def visual_sca():
    fn = "SCA/2001/sca20010630.tif"
    data = gdal.Open(fn, GA_ReadOnly).ReadAsArray()
    plt.imshow(data)
    plt.show()


def main():
    # reproject_all(type = "abo")
    # reproject_all(type = "sca")
    # check_sca_dates("abo")
    # check_sca_dates("sca")
    visual_sca()
if __name__ == "__main__":
    main()
