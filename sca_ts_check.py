__author__ = 'zeshi'
import numpy as np
import gdal
from gdal import GA_ReadOnly
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
from datetime import date, datetime, timedelta
import pandas as pd
import pickle

def form_cdec_sca_ts(start_year, end_year):
    cdec_station_df = pd.read_csv("cdec/cdec_stations.csv", header=0, delimiter=",")
    cdec_ts_df = pd.read_csv("cdec/cdec_all.csv")
    elev = [2667.0, 2621.28, 2438.4, 2316.48, 2316.48, 2164.08, 2042.16, 2011.68, 1798.32, 1706.88, 1609.344, 1569.72]
    veg = gdal.Open("vegetation/amr_nlcd2011_500m.tif", GA_ReadOnly).ReadAsArray()
    fveg = veg.astype(float) / 100.
    nrows = end_year - start_year + 1
    ncols = len(elev)

    # fig, axarr = plt.subplots(nrows = nrows, ncols = ncols, figsize=(50, 30))
    lat = cdec_station_df['lat'].as_matrix()
    lon = cdec_station_df['long'].as_matrix()

    date_conv_fun = lambda date_obj: str(date_obj.month).lstrip('0') + "/" + \
                                     str(date_obj.day).lstrip('0') + "/" + \
                                     str(date_obj.year)[2:]

    sca_gr = gdal.Open("SCA/2001/sca20010101.tif", GA_ReadOnly).GetGeoTransform()
    x_idx = np.round(((lon - sca_gr[0]) / sca_gr[1])).astype(int)
    y_idx = np.round(((lat - sca_gr[3]) / sca_gr[5])).astype(int)
    compare_ts_list = []
    for year in range(start_year, end_year + 1):
        temp_date = date(year, 4, 1)
        cdec_sca = np.empty([0, ncols])
        cdec_swe = np.empty([0, ncols])
        date_list = []
        while temp_date <= date(year, 6, 30):
            temp_sca_fn = "SCA/" + str(year) + "/sca" + \
                            temp_date.strftime("%Y%m%d") + ".tif"
            temp_date_str = date_conv_fun(temp_date)
            temp_swe_cdec = cdec_ts_df.loc[cdec_ts_df['Dates'] == temp_date_str]["SWE(inches)"].as_matrix() * 0.025

            temp_sca = gdal.Open(temp_sca_fn, GA_ReadOnly).ReadAsArray()
            temp_sca = temp_sca / (1 - fveg)
            temp_sca[temp_sca > 1.0] = 1.0
            temp_sca_cdec = temp_sca[y_idx, x_idx]
            cdec_sca = np.vstack((cdec_sca, temp_sca_cdec))
            cdec_swe = np.vstack((cdec_swe, temp_swe_cdec))
            date_list.append(temp_date)
            temp_date += timedelta(days=1)
        temp_df = pd.DataFrame({'dates': date_list})
        for idx in range(0, ncols):
            sca_header = "sca_" + str(idx+1).zfill(2)
            cdec_header = "cdec_" + str(idx+1).zfill(2)
            temp_df[sca_header] = cdec_sca[:, idx]
            temp_df[cdec_header] = cdec_swe[:, idx]
        print temp_df
        compare_ts_list.append(temp_df)
    pickle.dump(compare_ts_list, open("cdec_sca/cdec_sca_ts.p", 'wb'))

def cdec_sca_plot():
    compare_ts_list = pickle.load(open("cdec_sca/cdec_sca_ts.p", 'rb'))
    elev = [2667.0, 2621.28, 2438.4, 2316.48, 2316.48, 2164.08, 2042.16, 2011.68, 1798.32, 1706.88, 1609.344, 1569.72]
    nrows = len(compare_ts_list)
    ncols = len(elev)
    fig, axarr = plt.subplots(nrows=nrows, ncols=ncols, figsize=(50, 30))
    for row_idx in range(0, nrows):
        temp_df = compare_ts_list[row_idx]
        for col_idx in range(0, ncols):
            sca_header = "sca_" + str(col_idx+1).zfill(2)
            cdec_header = "cdec_" + str(col_idx+1).zfill(2)
            dates = temp_df["dates"]
            sca = temp_df[sca_header]
            cdec = temp_df[cdec_header]
            (l1, l2) = axarr[row_idx, col_idx].plot(dates, sca, "b", dates, cdec, "--k")
            axarr[row_idx, col_idx].set_ylim([0, 1.2])
            axarr[row_idx, col_idx].set_xlim([datetime(year=row_idx+2001, month=4, day=1), date(year=row_idx+2001, month=6, day=30)])
            months = mdates.MonthLocator(bymonthday=15, interval=1)
            monthsFmt = mdates.DateFormatter('%b')
            axarr[row_idx, col_idx].xaxis.set_major_locator(months)
            axarr[row_idx, col_idx].xaxis.set_major_formatter(monthsFmt)
            axarr[row_idx, col_idx].grid(True)
            axarr[row_idx, col_idx].text(dates[3], 1, str(row_idx+2001))
            if row_idx >= 1:
                plt.setp(axarr[row_idx, col_idx].get_yticklabels()[-1], visible=False)
            if row_idx < nrows - 1:
                plt.setp(axarr[row_idx, col_idx].get_xticklabels(), visible=False)
            if col_idx >= 1:
                plt.setp(axarr[row_idx, col_idx].get_yticklabels(), visible=False)
            if row_idx == 0 and col_idx == 0:
                axarr[row_idx, col_idx].legend((l1, l2), ("SCA", "CDEC SWE"), frameon=False, loc = "upper right", prop={'size':10})
            if row_idx == 0:
                axarr[row_idx, col_idx].set_title(str(elev[col_idx]))
    fig.subplots_adjust(hspace=0.)
    plt.show()


def form_sca_elev_band_ts():
    pass

form_cdec_sca_ts(2001, 2011)
cdec_sca_plot()