__author__ = 'zeshi'

import gdal
from gdal import GA_ReadOnly
from datetime import datetime, date, timedelta
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pandas as pd

def visualize_check(year):
    cdec_station_df = pd.read_csv("cdec/cdec_stations.csv", header=0, delimiter=",")
    start_date = date(year, 4, 1)
    end_date = date(year, 4, 1)
    temp_date = start_date
    lat = cdec_station_df['lat'].as_matrix()
    lon = cdec_station_df['long'].as_matrix()
    while temp_date <= end_date:
        check_fn = "output/swe_500m_reset_zero_hourly/" + str(year) + "/" + \
                   str(year) + str(temp_date.month).zfill(2) + \
                   str(temp_date.day).zfill(2) + "_swe.tif"
        array = gdal.Open(check_fn, GA_ReadOnly).ReadAsArray()
        georef = gdal.Open(check_fn, GA_ReadOnly).GetGeoTransform()
        x_idx = (lon - georef[0])/georef[1]
        y_idx = (lat - georef[3])/georef[5]
        plt.imshow(array, interpolation='None')
        plt.colorbar()
        plt.scatter(x_idx, y_idx, s=20)
        # plt.scatter(254, 122, s=20, c='r')
        plt.title("Zeshi")
        plt.show()
        print array[119, 251]
        temp_date += timedelta(days=1)

def visualize_check_noah(year):
    start_date = date(year, 4, 1)
    end_date = date(year, 4, 1)
    temp_date = start_date
    while temp_date <= end_date:
        check_fn = "2012_2014_500m/" + str(year) + "/" + str(temp_date.day).zfill(2) + \
                   temp_date.strftime("%B")[0:3].upper() + str(year) + ".tif"
        array = gdal.Open(check_fn, GA_ReadOnly).ReadAsArray()[480:668, 240:598]
        plt.imshow(array)
        plt.colorbar()
        plt.title("Noah")
        plt.show()
        temp_date += timedelta(days=1)


def visualize_both(year):
    zeshi_fn = "output/swe_500m_reset_zero_hourly/" + str(2001) + "/" + str(2001) + str(4).zfill(2) + \
                   str(1).zfill(2) + "_swe.tif"
    noah_fn = "2012_2014_500m/" + str(year) + "/" + str(1).zfill(2) + \
                   "APR" + str(year) + ".tif"
    zeshi = gdal.Open(zeshi_fn, GA_ReadOnly).ReadAsArray()
    noah = gdal.Open(noah_fn, GA_ReadOnly).ReadAsArray()[480:668, 240:598]
    min = np.min([np.min(zeshi), np.min(noah)])
    max = np.max([np.max(zeshi), np.max(noah)])
    plt.subplot(2,1,1)
    plt.imshow(zeshi, vmin = min, vmax = max)
    plt.subplot(2,1,2)
    plt.imshow(noah, vmin = min, vmax = max)
    plt.colorbar()
    plt.show()



def visualize_compare(year):
    zeshi_prev_fn = "output/swe_500m_reset_zero_hourly/" + str(year) + "/" + str(year) + str(4).zfill(2) + \
                   str(1).zfill(2) + "_swe.tif"
    zeshi_post_fn = "output/swe_500m_reset_zero_hourly/" + str(year) + "/" + str(year) + str(4).zfill(2) + \
                   str(2).zfill(2) + "_swe.tif"
    noah_prev_fn = "2012_2014_500m/" + str(year) + "/" + str(1).zfill(2) + \
                   "APR" + str(year) + ".tif"
    noah_post_fn = "2012_2014_500m/" + str(year) + "/" + str(2).zfill(2) + \
                   "APR" + str(year) + ".tif"
    from skimage.transform import resize
    zeshi_prev = gdal.Open(zeshi_prev_fn).ReadAsArray().astype(float)
    zeshi_post = gdal.Open(zeshi_post_fn).ReadAsArray().astype(float)
    noah_prev = gdal.Open(noah_prev_fn).ReadAsArray()[480:668, 240:598].astype(float)
    noah_post = gdal.Open(noah_post_fn).ReadAsArray()[480:668, 240:598].astype(float)
    zeshi_diff = zeshi_prev - zeshi_post
    noah_diff = noah_prev - noah_post
    noah = resize(noah_diff, zeshi_diff.shape)
    scale = np.zeros(zeshi_diff.shape)
    scale[zeshi_diff != 0.] = noah[zeshi_diff != 0.] / zeshi_diff[zeshi_diff != 0.]
    scale[scale > 100] = 0
    plt.imshow(scale)
    plt.colorbar()
    plt.show()

def calc_noah_ar_idx():
    noah = gdal.Open("2012_2014_500m/2012/01APR2012.tif", GA_ReadOnly)
    noah_array = noah.ReadAsArray()
    noah_geotrans = noah.GetGeoTransform()
    zeshi = gdal.Open("output/swe_500m/2012/20120401_swe.tif", GA_ReadOnly)
    zeshi_array = zeshi.ReadAsArray()
    zeshi_geotrans = zeshi.GetGeoTransform()
    zeshi_ulx = zeshi_geotrans[0]
    zeshi_uly = zeshi_geotrans[3]
    zeshi_llx = zeshi_ulx + zeshi_geotrans[1] * zeshi_array.shape[1]
    zeshi_lly = zeshi_uly + zeshi_geotrans[5] * zeshi_array.shape[0]
    print noah_geotrans
    print zeshi_ulx, zeshi_uly, zeshi_llx, zeshi_lly
    print int((zeshi_uly - noah_geotrans[3]) / noah_geotrans[5]), int((zeshi_ulx - noah_geotrans[0]) / noah_geotrans[1]), \
        int((zeshi_lly - noah_geotrans[3]) / noah_geotrans[5]),int((zeshi_llx - noah_geotrans[0]) / noah_geotrans[1])

def plot():
    data = gdal.Open("vegetation/amr_nlcd2011_500m.tif", GA_ReadOnly).ReadAsArray()
    plt.imshow(data)
    plt.colorbar()
    plt.show()

class Formatter(object):
    def __init__(self, im):
        self.im = im
    def __call__(self, x, y):
        z = self.im.get_array()[int(y), int(x)]
        return 'x={:.01f}, y={:.01f}, z={:.01f}'.format(x, y, z)

def cdec_compare(start_year, end_year, idx, elev=None):
    cdec_station_df = pd.read_csv("cdec/cdec_stations.csv", header=0, delimiter=",")
    cdec_ts_df = pd.read_csv("cdec/cdec_all.csv")
    if elev != None:
        cdec_station_df = cdec_station_df.loc[cdec_station_df['elevation']==elev]
    else:
        cdec_station_df = cdec_station_df.loc[cdec_station_df['elevation']>=2300.0]
    cdec_ts_df = cdec_ts_df.loc[cdec_ts_df['Station'].isin(cdec_station_df['station_id'].tolist())]
    lat = cdec_station_df['lat'].as_matrix()
    lon = cdec_station_df['long'].as_matrix()
    date_conv_fun = lambda date_obj: str(date_obj.month).lstrip('0') + "/" + \
                                     str(date_obj.day).lstrip('0') + "/" + \
                                     str(date_obj.year)[2:]
    swe_gr = gdal.Open("output/swe_500m_reset_zero_hourly/2001/20010401_swe.tif", GA_ReadOnly).GetGeoTransform()
    x_idx = np.round(((lon - swe_gr[0]) / swe_gr[1])).astype(int)
    y_idx = np.round(((lat - swe_gr[3]) / swe_gr[5])).astype(int)
    compare_ts_list = []
    for year in range(start_year, end_year + 1):
        temp_date = date(year, 4, 1)
        recon_swe = []
        cdec_swe = []
        date_list = []
        while temp_date <= date(year, 6, 30):
            temp_recon_fn = "output/swe_500m_reset_zero_hourly/" + str(year) + "/" + \
                            temp_date.strftime("%Y%m%d") + "_swe.tif"
            temp_date_str = date_conv_fun(temp_date)
            temp_swe_cdec = cdec_ts_df.loc[cdec_ts_df['Dates'] == temp_date_str]["SWE(inches)"].as_matrix() * 0.025

            # Just SCN
            if idx != False:
                query_idx = idx
                temp_x_idx = x_idx[idx]
                temp_y_idx = y_idx[idx]
            else:
            # All sites
                query_idx = np.isnan(temp_swe_cdec)
                temp_x_idx = x_idx[query_idx == False]
                temp_y_idx = y_idx[query_idx == False]

            temp_swe_recon = gdal.Open(temp_recon_fn, GA_ReadOnly).ReadAsArray()
            temp_swe_recon_cdec = temp_swe_recon[temp_y_idx, temp_x_idx]
            if idx != False:
                temp_swe_cdec = temp_swe_cdec[idx]
            else:
                temp_swe_cdec = temp_swe_cdec[query_idx == False]
            recon_swe.append(np.nanmean(temp_swe_recon_cdec))
            cdec_swe.append(np.nanmean(temp_swe_cdec))
            date_list.append(temp_date)
            temp_date += timedelta(days=1)
        compare_ts_list.append(pd.DataFrame({'dates': date_list, 'recon': recon_swe, 'cdec': cdec_swe}))
    nrows = int(round(float(len(compare_ts_list)) / 2.0))
    fig, axarr = plt.subplots(nrows = nrows, ncols = 2, figsize=(15, 20))
    if len(compare_ts_list) % 2 == 1:
        fig.delaxes(axarr[nrows-1, 1])
    for i, temp_ts in enumerate(compare_ts_list):
        row_idx = i % nrows
        col_idx = i / nrows
        l1, l2 = axarr[row_idx, col_idx].plot(temp_ts["dates"], temp_ts["recon"], 'b', temp_ts["dates"], temp_ts["cdec"], '--k')
        axarr[row_idx, col_idx].set_ylim([0, 1.2])
        axarr[row_idx, col_idx].set_xlim([datetime(year=i+2001, month=4, day=1), date(year=i+2001, month=6, day=30)])
        months = mdates.MonthLocator(bymonthday=15, interval=1)
        monthsFmt = mdates.DateFormatter('%b')
        axarr[row_idx, col_idx].xaxis.set_major_locator(months)
        axarr[row_idx, col_idx].xaxis.set_major_formatter(monthsFmt)
        axarr[row_idx, col_idx].grid(True)
        axarr[row_idx, col_idx].text(temp_ts["dates"][3], 1, str(i+2001))
        # if row_idx < nrows - 1 and i != len(compare_ts_list) - 1:
        #     plt.setp(axarr[row_idx, col_idx].get_xticklabels(), visible=False)
        if row_idx >= 1:
            plt.setp(axarr[row_idx, col_idx].get_yticklabels()[-1], visible=False)
        if i == 0:
            axarr[row_idx, col_idx].legend((l1, l2), ("Reconstruction", "Snow pillow"), frameon=False, loc = "upper right")
    fig.subplots_adjust(hspace=0.)
    # fig.legend((l1, l2), ("Reconstruction", "Snow pillow"), "lower right")
    plt.show()

def cdec_compare_site_wise(start_year, end_year, idx):
    cdec_station_df_all = pd.read_csv("cdec/cdec_stations.csv", header=0, delimiter=",")
    cdec_ts_df_all = pd.read_csv("cdec/cdec_all.csv")
    elev = [2667.0, 2621.28, 2438.4, 2316.48, 2316.48, 2164.08, 2042.16, 2011.68, 1798.32, 1706.88, 1609.344, 1569.72]
    nrows = end_year - start_year + 1
    ncols = len(elev)
    fig, axarr = plt.subplots(nrows = nrows, ncols = ncols, figsize=(50, 30))
    for j, temp_elev in enumerate(elev):
        cdec_station_df = cdec_station_df_all.loc[cdec_station_df_all['elevation']==temp_elev]
        cdec_ts_df = cdec_ts_df_all.loc[cdec_ts_df_all['Station'].isin(cdec_station_df['station_id'].tolist())]
        lat = cdec_station_df['lat'].as_matrix()
        lon = cdec_station_df['long'].as_matrix()
        date_conv_fun = lambda date_obj: str(date_obj.month).lstrip('0') + "/" + \
                                         str(date_obj.day).lstrip('0') + "/" + \
                                         str(date_obj.year)[2:]
        swe_gr = gdal.Open("output/swe_500m_reset_zero_hourly/2001/20010401_swe.tif", GA_ReadOnly).GetGeoTransform()
        x_idx = np.round(((lon - swe_gr[0]) / swe_gr[1])).astype(int)
        y_idx = np.round(((lat - swe_gr[3]) / swe_gr[5])).astype(int)
        compare_ts_list = []
        for year in range(start_year, end_year + 1):
            temp_date = date(year, 4, 1)
            recon_swe = []
            cdec_swe = []
            date_list = []
            while temp_date <= date(year, 6, 30):
                temp_recon_fn = "output/swe_500m_reset_zero_hourly/" + str(year) + "/" + \
                                temp_date.strftime("%Y%m%d") + "_swe.tif"
                temp_date_str = date_conv_fun(temp_date)
                temp_swe_cdec = cdec_ts_df.loc[cdec_ts_df['Dates'] == temp_date_str]["SWE(inches)"].as_matrix() * 0.025

                # Just SCN
                if idx != False:
                    query_idx = idx
                    temp_x_idx = x_idx[idx]
                    temp_y_idx = y_idx[idx]
                else:
                # All sites
                    query_idx = np.isnan(temp_swe_cdec)
                    temp_x_idx = x_idx[query_idx == False]
                    temp_y_idx = y_idx[query_idx == False]

                temp_swe_recon = gdal.Open(temp_recon_fn, GA_ReadOnly).ReadAsArray()
                temp_swe_recon_cdec = temp_swe_recon[temp_y_idx, temp_x_idx]
                if idx != False:
                    temp_swe_cdec = temp_swe_cdec[idx]
                else:
                    temp_swe_cdec = temp_swe_cdec[query_idx == False]
                recon_swe.append(np.nanmean(temp_swe_recon_cdec))
                cdec_swe.append(np.nanmean(temp_swe_cdec))
                date_list.append(temp_date)
                temp_date += timedelta(days=1)
            compare_ts_list.append(pd.DataFrame({'dates': date_list, 'recon': recon_swe, 'cdec': cdec_swe}))
        for i, temp_ts in enumerate(compare_ts_list):
            row_idx = i
            col_idx = j
            l1, l2 = axarr[row_idx, col_idx].plot(temp_ts["dates"], temp_ts["recon"], 'b', temp_ts["dates"], temp_ts["cdec"], '--k')
            axarr[row_idx, col_idx].set_ylim([0, 1.2])
            axarr[row_idx, col_idx].set_xlim([datetime(year=i+2001, month=4, day=1), date(year=i+2001, month=6, day=30)])
            months = mdates.MonthLocator(bymonthday=15, interval=1)
            monthsFmt = mdates.DateFormatter('%b')
            axarr[row_idx, col_idx].xaxis.set_major_locator(months)
            axarr[row_idx, col_idx].xaxis.set_major_formatter(monthsFmt)
            axarr[row_idx, col_idx].grid(True)
            axarr[row_idx, col_idx].text(temp_ts["dates"][3], 1, str(i+2001))
            # if row_idx < nrows - 1 and i != len(compare_ts_list) - 1:
            #     plt.setp(axarr[row_idx, col_idx].get_xticklabels(), visible=False)
            if row_idx >= 1:
                plt.setp(axarr[row_idx, col_idx].get_yticklabels()[-1], visible=False)
            if row_idx < nrows - 1:
                plt.setp(axarr[row_idx, col_idx].get_xticklabels(), visible=False)
            if col_idx >= 1:
                plt.setp(axarr[row_idx, col_idx].get_yticklabels(), visible=False)
            if i == 0 and j == 0:
                axarr[row_idx, col_idx].legend((l1, l2), ("Reconstruction", "Snow pillow"), frameon=False, loc = "upper right", prop={'size':10})
            if i == 0:
                axarr[row_idx, col_idx].set_title(str(temp_elev))
        fig.subplots_adjust(hspace=0.)
        # fig.legend((l1, l2), ("Reconstruction", "Snow pillow"), "lower right")
    plt.show()

def cdec_correlate_compare():
    cdec_station_df = pd.read_csv("cdec/cdec_stations.csv", header=0, delimiter=",")
    cdec_ts_df = pd.read_csv("cdec/cdec_all.csv")
    lat = cdec_station_df['lat'].as_matrix()
    lon = cdec_station_df['long'].as_matrix()
    date_conv_fun = lambda date_obj: str(date_obj.month).lstrip('0') + "/" + \
                                     str(date_obj.day).lstrip('0') + "/" + \
                                     str(date_obj.year)[2:]
    swe_gr = gdal.Open("output/swe_500m_reset_zero_hourly/2001/20010401_swe.tif", GA_ReadOnly).GetGeoTransform()
    x_idx = np.round(((lon - swe_gr[0]) / swe_gr[1])).astype(int)
    y_idx = np.round(((lat - swe_gr[3]) / swe_gr[5])).astype(int)
    recon_swe_array = np.array([0])
    cdec_swe_array = np.array([0])
    year_array = np.array([0])
    for year in range(2001, 2012):
        temp_date = date(year, 4, 1)
        while temp_date <= date(year, 6, 30):
            temp_recon_fn = "output/swe_500m_reset_zero_hourly/" + str(year) + "/" + \
                            temp_date.strftime("%Y%m%d") + "_swe.tif"
            temp_date_str = date_conv_fun(temp_date)
            temp_swe_cdec = cdec_ts_df.loc[cdec_ts_df['Dates'] == temp_date_str]["SWE(inches)"].as_matrix() * 0.025
            query_idx = np.isnan(temp_swe_cdec)
            temp_x_idx = x_idx[query_idx == False]
            temp_y_idx = y_idx[query_idx == False]
            temp_swe_recon = gdal.Open(temp_recon_fn, GA_ReadOnly).ReadAsArray()
            temp_swe_recon_cdec = temp_swe_recon[temp_y_idx, temp_x_idx]
            temp_swe_cdec = temp_swe_cdec[query_idx == False]
            recon_swe_array = np.append(recon_swe_array, temp_swe_recon_cdec)
            cdec_swe_array = np.append(cdec_swe_array, temp_swe_cdec)
            year_array = np.append(year_array, np.ones(len(temp_swe_cdec))*year)
            temp_date += timedelta(days=1)
    recon_swe_array = recon_swe_array[cdec_swe_array > 0]
    year_array = year_array[cdec_swe_array > 0]
    cdec_swe_array = cdec_swe_array[cdec_swe_array > 0]
    cdec_swe_array = cdec_swe_array[recon_swe_array > 0]
    year_array = year_array[recon_swe_array > 0]
    recon_swe_array = recon_swe_array[recon_swe_array > 0]
    year_2009 = year_array[year_array == 2011]
    cdec_2009 = cdec_swe_array[year_array == 2011]
    recon_2009 = recon_swe_array[year_array == 2011  ]
    plt.scatter(cdec_swe_array, recon_swe_array, c=year_array)
    plt.xlim([0, 2.5])
    plt.ylim([0, 2.5])
    plt.show()


def cdec_elevBand_compare():
    elevBands = [1490, 1650, 1930, 2060, 2256, 2550]
    dem_ds = gdal.Open("DEM/ws_dem_500m.tif", GA_ReadOnly)
    dem = dem_ds.ReadAsArray()
    cdec_station_df = pd.read_csv("cdec/cdec_stations.csv", header=0, delimiter=",")
    cdec_ts_df = pd.read_csv("cdec/cdec_all.csv")
    date_conv_fun = lambda date_obj: str(date_obj.month).lstrip('0') + "/" + \
                                     str(date_obj.day).lstrip('0') + "/" + \
                                     str(date_obj.year)[2:]

    f, axarr = plt.subplots(nrows=6, ncols=11, sharex=False, figsize=(40, 20))
    for year in range(2001, 2012):
        swe_df_list = []
        for temp_elev in elevBands:
            temp_elev_cdec_stations = cdec_station_df.loc[(cdec_station_df['elevation'] > temp_elev) &
                                                          (cdec_station_df['elevation'] < (temp_elev + 200))]
            temp_cdec_id = temp_elev_cdec_stations['station_id'].tolist()
            temp_elev_cdec_ts = cdec_ts_df.loc[cdec_ts_df['Station'].isin(temp_cdec_id)]
            temp_date = date(year, 4, 1)
            end_date = date(year, 6, 30)
            dates = []
            cdec_SWE_mean = []
            recon_SWE_mean = []
            recon_SWE_std = []
            while (temp_date <= end_date):
                temp_recon_fn = "output/swe_500m_reset_zero_hourly/" + str(year) + "/" + \
                                temp_date.strftime("%Y%m%d") + "_swe.tif"
                temp_date_str = date_conv_fun(temp_date)
                temp_swe_cdec = temp_elev_cdec_ts.loc[temp_elev_cdec_ts['Dates'] == temp_date_str]["SWE(inches)"].as_matrix() * 0.025
                temp_swe_recon = gdal.Open(temp_recon_fn, GA_ReadOnly).ReadAsArray()
                temp_swe_recon = temp_swe_recon[(dem > temp_elev) & (dem < temp_elev + 200)]
                recon_SWE_mean.append(np.nanmean(temp_swe_recon))
                recon_SWE_std.append(np.nanstd(temp_swe_recon))
                cdec_SWE_mean.append(np.nanmean(temp_swe_cdec))
                dates.append(temp_date)
                temp_date += timedelta(days=1)
            swe_df_list.append(pd.DataFrame({'dates': dates, 'recon_mean': recon_SWE_mean, 'recon_std': recon_SWE_std,
                                             'cdec': cdec_SWE_mean}))
        num_subplots = len(swe_df_list)

        for i in range(0, num_subplots):
            months = mdates.MonthLocator(bymonthday=15, interval=1)
            monthsFmt = mdates.DateFormatter('%b')

            temp_df = swe_df_list[i]
            l1 = axarr[i, year-2001].plot(temp_df['dates'], temp_df['recon_mean'], 'b', label="Reconstruction")
            l2 = axarr[i, year-2001].plot(temp_df['dates'], temp_df['cdec'], '--k', label="Snow pillows")
            axarr[i, year-2001].xaxis.set_major_locator(months)
            axarr[i, year-2001].xaxis.set_major_formatter(monthsFmt)
            dates_list = temp_df['dates'].tolist()
            y = temp_df['recon_mean'].as_matrix()
            err = temp_df['recon_std'].as_matrix()
            axarr[i, year-2001].fill_between(dates_list, y-2*err, y+2*err, color='b', alpha=0.2)
            ylab = str(elevBands[i]) + "m-" + str(elevBands[i]+200) + "m"
            if year == 2001:
                axarr[i, year-2001].set_ylabel(ylab)
            else:
                plt.setp(axarr[i, year-2001].get_yticklabels(), visible=False)
            if i < num_subplots-1:
                plt.setp(axarr[i, year-2001].get_xticklabels(), visible=False)
            if i == 0:
                axarr[i, year-2001].set_ylim([0, 1.2])
                axarr[i, year-2001].set_title(str(year))
            elif i == 1:
                axarr[i, year-2001].set_ylim([0, 1.4])
            elif i == 2:
                axarr[i, year-2001].set_ylim([0, 2.0])
            elif i == 3:
                axarr[i, year-2001].set_ylim([0, 2.0])
            elif i == 4:
                axarr[i, year-2001].set_ylim([0, 2.5])
            elif i == 5:
                axarr[i, year-2001].set_ylim([0, 2.5])
            if i >= 1 and year == 2001:
                plt.setp(axarr[i, year-2001].get_yticklabels()[-1], visible=False)
            if i == 0 and year == 2001:
                axarr[i, year-2001].legend(frameon=False, loc = "upper right", prop={'size': 10})
            axarr[i, year-2001].grid(True)
    f.subplots_adjust(hspace=0.)
    plt.show()

def peak_swe_elev_analysis(year, bin_size=10):
        print year
        dem = gdal.Open("DEM/ws_dem_500m.tif", GA_ReadOnly).ReadAsArray()
        dem[dem < 0] = float('nan')
        dem_min = np.nanmin(dem)
        dem_max = np.nanmax(dem)
        swe_fn = "output/swe_500m_reset_zero_hourly/" + str(year) + "/" + str(year) + "0401_swe.tif"
        swe = gdal.Open(swe_fn, GA_ReadOnly).ReadAsArray()
        result = np.empty((0, 3))
        for temp_elev in range(int(dem_min), int(dem_max) + bin_size, bin_size):
            temp_swe = swe[(dem >= temp_elev) & (dem <= (temp_elev + bin_size))]
            temp_swe_avg = np.nanmean(temp_swe)
            temp_swe_std = np.nanstd(temp_swe)
            result = np.vstack((result, [temp_elev + int(float(bin_size) / 2.), temp_swe_avg, temp_swe_std]))
        return result

def peak_swe_elev_plot(start_year, end_year, bin_size):
    def generate_plot(year_list, analysis_results):
        num_years = len(year_list)
        if num_years <= 5:
            nrows = num_years
            ncols = 1
        else:
            nrows = 5
            ncols = int(np.ceil(num_years / 5.))
        fig, axarr = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, figsize = (10 * ncols, 15))
        for year, result in zip(year_list, analysis_results):
            row_idx = (year - start_year) % 5
            col_idx = (year - start_year) / 5
            axarr[row_idx, col_idx].plot(result[:, 0], result[:, 1], 'b', linewidth = 3)
            axarr[row_idx, col_idx].fill_between(result[:, 0], result[:, 1]-result[:, 2], result[:, 1]+result[:, 2],
                                                 color='b', alpha=0.2)
            axarr[row_idx, col_idx].set_ylim([0, 2])
            if col_idx >= 1:
                plt.setp(axarr[row_idx, col_idx].get_yticklabels(), visible=False)
            if row_idx <= 3:
                plt.setp(axarr[row_idx, col_idx].get_xticklabels(), visible=False)
            if row_idx >= 1:
                plt.setp(axarr[row_idx, col_idx].get_yticklabels()[-1], visible=False)

            axarr[row_idx, col_idx].text(np.min(result[:, 0]), 1.5, str(year))
            axarr[row_idx, col_idx].grid(True)
        fig.subplots_adjust(hspace=0.)
        fig.text(0.5, 0.06, "Elevation, m", ha='center')
        fig.text(0.08, 0.5, "Snow water equivalent, m", va='center', rotation='vertical')
        plt.show()

    from multiprocessing import Pool, cpu_count
    from functools import partial
    print "partial"

    part_peak_swe_elev_analysis = partial(peak_swe_elev_analysis, bin_size=bin_size)
    year_list = range(start_year, end_year + 1)
    pool = Pool(processes=cpu_count())
    print "parallel"
    peak_swe_elev_analysis_results = pool.map(part_peak_swe_elev_analysis, year_list)
    pool.close()
    pool.join()
    generate_plot(year_list, peak_swe_elev_analysis_results)


# def plot_heat():
#     wind_redist_obj = wind_redistribution(res=500)
#     sw_ds = short_wave_ds(res=500)
#     obj = delta_swe_hourly(2012, 4, 1, 6, sw_ds, wind_redist_obj, res=500)
#     for i in range(0, 3):
#         if i == 0:
#             array = obj.LH
#         elif i == 1:
#             array = obj.SH
#         else:
#             array = obj.CDLW
#         plt.imshow(array)
#         plt.colorbar()
#         plt.show()

def main():
    peak_swe_elev_plot(2001, 2010, 10)
    # for year in range(2001, 2012):
    # #     cdec_compare(year, year + 1, False)
    #     visualize_check(year)
    # cdec_elevBand_compare()
    # elev = [2667.0, 2621.28, 2438.4, 2316.48, 2316.48, 2164.08, 2042.16, 2011.68, 1798.32, 1706.88, 1609.344, 1569.72]
    # for temp_elev in elev:
    #     cdec_compare(2001, 2011, False, elev=temp_elev)
    # cdec_compare(2001, 2011, False)
    # cdec_compare_site_wise(2001, 2011, False)
    # cdec_correlate_compare()
    # for year in range(2001, 2013):
    #     visualize_check(year)
    # for i in range(0, 10):
    #     cdec_compare(i)
    # lat = 38.747
    # lon = -120.068
    # dem_500m_ds = gdal.Open("DEM/500m_dem.tif")
    # dem_500m_gr = dem_500m_ds.GetGeoTransform()
    # dem_500m = dem_500m_ds.ReadAsArray()
    # print dem_500m_gr
    # idx_y = (lat - dem_500m_gr[3]) / dem_500m_gr[5]
    # idx_x = (lon - dem_500m_gr[0]) / dem_500m_gr[1]
    # print idx_x, idx_y
    # print 119.0 * dem_500m_gr[5] + dem_500m_gr[3]
    # print 251.0 * dem_500m_gr[1] + dem_500m_gr[0]
    # plt.imshow(dem_500m)
    # plt.scatter(idx_x, idx_y, s=1)
    # plt.show()
    # for year in range(2001, 2012):
    #     visualize_check(year)
    # visualize_check_noah(2012)
    # visualize_compare(2012)
    # visualize_check(2006)

if __name__ == "__main__":
    main()

