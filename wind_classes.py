__author__ = 'zeshi'

import numpy as np
from datetime import date, timedelta
from sklearn.linear_model import LinearRegression as lr
from matplotlib import pyplot as plt
import pickle

class Wind:
    """
    Super class of all wind related data
    has function calculating the yearly average
    anomaly through time series
    and gap filling function inside
    """
    def __init__(self, lat, lon, station_id, type, date, wind_speed, wind_dir):
        self.lat = lat
        self.lon = lon
        self.station_id = station_id
        self.type = type
        self.date = date
        self.wind_speed = wind_speed
        self.wind_dir = wind_dir
        self.wind_speed_yday_avg = None
        self.wind_dir_yday_avg = None
        self.wind_speed_anomaly = None
        self.wind_dir_anomaly = None
        self.yday = None

    def set_yday(self):
        self.yday = np.zeros(len(self.date)).astype(int)
        for i, temp_date in enumerate(self.date):
            self.yday[i] = temp_date.timetuple().tm_yday

    def cal_yearly_avg(self):
        assert self.yday is not None, "Need to calcualte yday first!"
        yday = range(1, 367)
        self.wind_speed_yday_avg = np.zeros(len(yday)).astype(float)
        self.wind_dir_yday_avg = np.zeros(len(yday)).astype(float)
        for i,temp_yday in enumerate(yday):
            idx = np.where(self.yday == temp_yday)
            self.wind_speed_yday_avg[i] = np.nanmean(self.wind_speed[idx])
            self.wind_dir_yday_avg[i] = np.nanmean(self.wind_dir[idx])

    def cal_anomaly(self):
        assert self.wind_speed_yday_avg is not None and self.wind_dir_yday_avg is not None, "Need to calculate yday " \
                                                                                    "average first!"
        self.wind_speed_anomaly = np.zeros(len(self.yday))
        self.wind_dir_anomaly = np.zeros(len(self.yday))
        yday = range(1, 367)
        for i, temp_yday in enumerate(yday):
            idx = np.where(self.yday == temp_yday)
            self.wind_speed_anomaly[idx] = self.wind_speed[idx] - self.wind_speed_yday_avg[i]
            self.wind_dir_anomaly[idx] = self.wind_dir[idx] - self.wind_dir_yday_avg[i]

    def ts_plot(self):
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        ax1.plot(self.date, self.wind_speed)
        ax1.set_ylabel("Wind speed ,m/s")
        ax2.plot(self.date, self.wind_dir)
        ax2.set_ylabel("Wind direction, degree")
        ax1.grid()
        ax2.grid()
        plt.show()

class Nldas(Wind):
    """
    sub class of Wind, specifically for nldas wind data
    """
    def __init__(self, lat, lon, station_id, date, wind_speed, wind_dir):
        Wind.__init__(self, lat, lon, station_id, "nldas", date, wind_speed, wind_dir)

class Asos(Wind):
    """
    sub class of Wind, specifically for Asos wind data
    has gap filling method with it
    """
    def __init__(self, lat, lon, station_id, date, wind_speed, wind_dir):
        Wind.__init__(self, lat, lon, station_id, "asos", date, wind_speed, wind_dir)

    def set_new_yday(self, date_array):
        yday = np.zeros(len(date_array)).astype(int)
        for i, temp_date in enumerate(date_array):
            yday[i] = temp_date.timetuple().tm_yday
        return yday

    def nldas_correlate(self):
        nldas_list = pickle.load(open("wind/nldas.p", "rb"))
        nldas = None
        for temp in nldas_list:
            if self.station_id == temp.station_id:
                nldas = temp
                break
        nldas_idx = np.where(np.logical_and(nldas.date >= np.min(self.date), nldas.date <= np.max(self.date)))
        nldas_wind_speed_anomaly = nldas.wind_speed_anomaly[nldas_idx]
        nldas_wind_dir_anomaly = nldas.wind_dir_anomaly[nldas_idx]

        fit_lr = lr()
        # mask1 = self.reject_outliers(self.wind_speed_anomaly)
        mask1 = ~np.isnan(self.wind_speed_anomaly)
        fit_lr.fit(nldas_wind_speed_anomaly[mask1].reshape((len(nldas_wind_speed_anomaly[mask1]), 1)), self.wind_speed_anomaly[mask1])
        result1 = fit_lr.predict(nldas_wind_speed_anomaly[mask1].reshape((len(nldas_wind_speed_anomaly[mask1]), 1)))
        std = np.sqrt(np.sum((self.wind_speed_anomaly[mask1] - result1) ** 2) / (len(result1) - 2))
        print "Standard deviation of the wind speed estimate is", std

        fit_lr = lr()
        # mask2 = self.reject_outliers(self.wind_dir_anomaly)
        mask2 = ~np.isnan(self.wind_dir_anomaly)
        fit_lr.fit(nldas_wind_dir_anomaly[mask2].reshape((len(nldas_wind_dir_anomaly[mask2]), 1)), self.wind_dir_anomaly[mask2])
        result2 = fit_lr.predict(nldas_wind_dir_anomaly[mask2].reshape((len(nldas_wind_dir_anomaly[mask2]), 1)))
        std = np.sqrt(np.sum((self.wind_dir_anomaly[mask2] - result2) ** 2) / (len(result2) - 2))
        print "Standard deviation of the wind direction estimate is", std

        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax1.plot(nldas_wind_speed_anomaly[mask1], self.wind_speed_anomaly[mask1], '.b')
        ax1.plot(nldas_wind_speed_anomaly[mask1], result1, '-r')
        ax2 = fig.add_subplot(212)
        ax2.plot(nldas_wind_dir_anomaly[mask2], self.wind_dir_anomaly[mask2], '.g')
        ax2.plot(nldas_wind_dir_anomaly[mask2], result2, '-r')
        plt.show()

    def linear_model(self, nldas_wind, type = 'speed'):
        X = nldas_wind
        if type == 'speed':
            y = self.wind_speed_anomaly
        else:
            y = self.wind_dir_anomaly
        mask = ~np.isnan(y)
        X = X[mask].reshape((len(X[mask]), 1))
        y = y[mask]
        lr_model = lr()
        lr_model.fit(X, y)
        est_y = lr_model.predict(X)
        std = np.sqrt(np.sum((est_y - y) ** 2) / (len(y) - 2))
        return lr_model, std

    def gap_filling(self):
        nldas_list = pickle.load(open("wind/nldas.p", "rb"))
        nldas = None
        for temp in nldas_list:
            if self.station_id == temp.station_id:
                nldas = temp
                break

        start = np.min(nldas.date)
        end = np.max(nldas.date)

        date_delta = (end - start).days
        new_date_list = [start + timedelta(days=x) for x in range(0, date_delta + 1)]
        new_date_array = np.array(new_date_list)
        new_yday = self.set_new_yday(new_date_array)
        new_wind_speed = np.zeros(len(new_date_array))
        new_wind_dir = np.zeros(len(new_date_array))
        nldas_wind_speed_anomaly_new_date = nldas.wind_speed_anomaly[np.where(np.logical_and(nldas.date >= start,
                                                                                             nldas.date <= end))]
        nldas_wind_dir_anomaly_new_date = nldas.wind_dir_anomaly[np.where(np.logical_and(nldas.date >= start,
                                                                                         nldas.date <= end))]

        min_old_date = np.min(self.date)
        max_old_date = np.max(self.date)
        nldas_wind_speed_anomaly_old_date = nldas.wind_speed_anomaly[np.where(np.logical_and(nldas.date >= min_old_date,
                                                                                    nldas.date <= max_old_date))]
        nldas_wind_dir_anomaly_old_date = nldas.wind_dir_anomaly[np.where(np.logical_and(nldas.date >= min_old_date,
                                                                                    nldas.date <= max_old_date))]

        wind_speed_anomaly_lr, speed_std = self.linear_model(nldas_wind_speed_anomaly_old_date, type='speed')

        wind_dir_anomaly_lr, dir_std = self.linear_model(nldas_wind_dir_anomaly_old_date, type='dir')

        for i, temp_new_date in enumerate(new_date_array):
            if temp_new_date in self.date:
                idx = np.where(self.date == temp_new_date)
                if ~np.isnan(self.wind_speed[idx])  and ~np.isnan(self.wind_dir[idx]):
                    new_wind_speed[i] = self.wind_speed[idx]
                    new_wind_dir[i] = self.wind_dir[idx]
                else:
                    new_wind_speed[i] = np.nan
                    new_wind_dir[i] = np.nan
            else:
                new_wind_speed[i] = np.nan
                new_wind_dir[i] = np.nan

        # substitution
        nan_new_dates_idx = np.where(np.isnan(new_wind_speed))[0]
        nan_new_yday = new_yday[nan_new_dates_idx]
        nan_new_nldas_wind_speed_anomaly = nldas_wind_speed_anomaly_new_date[nan_new_dates_idx]
        nan_new_nldas_wind_dir_anomaly = nldas_wind_dir_anomaly_new_date[nan_new_dates_idx]
        nan_new_wind_speed = np.zeros(len(nan_new_yday))
        nan_new_wind_dir = np.zeros(len(nan_new_yday))
        yday = np.arange(1, 367, 1)

        for temp_yday in yday:
            temp_nan_new_yday_idx = np.where(nan_new_yday == temp_yday)[0]
            if len(temp_nan_new_yday_idx) == 0:
                continue
            temp_length = len(temp_nan_new_yday_idx)
            temp_nan_new_nldas_wind_speed_anomaly = nan_new_nldas_wind_speed_anomaly[temp_nan_new_yday_idx]
            temp_nan_new_nldas_wind_dir_anomaly = nan_new_nldas_wind_dir_anomaly[temp_nan_new_yday_idx]
            X_1 = temp_nan_new_nldas_wind_speed_anomaly.reshape((temp_length, 1))
            X_2 = temp_nan_new_nldas_wind_dir_anomaly.reshape((temp_length, 1))
            nan_new_wind_speed[temp_nan_new_yday_idx] = self.wind_speed_yday_avg[temp_yday - 1] + \
                                                        wind_speed_anomaly_lr.predict(X_1)
                                                        # + np.random.normal(loc=0.0, scale=speed_std, size=temp_length)
            nan_new_wind_dir[temp_nan_new_yday_idx] = self.wind_dir_yday_avg[temp_yday - 1] + \
                                                      wind_dir_anomaly_lr.predict(X_2)
                                                      # + np.random.normal(loc=0.0, scale=dir_std, size=temp_length)

        new_wind_speed[nan_new_dates_idx] = nan_new_wind_speed
        new_wind_dir[nan_new_dates_idx] = nan_new_wind_dir

        self.date = new_date_array
        self.wind_speed = new_wind_speed
        self.wind_dir = new_wind_dir
        self.set_yday()
        self.cal_yearly_avg()
        self.cal_anomaly()

class Raws(Asos):
    """
    sub class of Wind, specifically for Raws wind data
    has gap filling method with it
    """
    def __init__(self, lat, lon, station_id, date, wind_speed, wind_dir):
        Wind.__init__(self, lat, lon, station_id, "raws", date, wind_speed, wind_dir)

def test_asos():
    asos = pickle.load(open("formated_asos_wind/asos_info.p", "rb"))
    for site in asos:
        new_site = Asos(site.lat, site.lon, site.station_id, site.daily_wind.date, site.daily_wind.wind_speed,
                        site.daily_wind.dir)
        new_site.set_yday()
        new_site.cal_yearly_avg()
        new_site.cal_anomaly()
        print np.where(np.isnan(new_site.wind_speed_yday_avg))
        yday_range = np.unique(new_site.yday)
        yday_range.sort()
        fig = plt.figure()
        ax1 = fig.add_subplot(411)
        ax1.plot(yday_range, new_site.wind_speed_yday_avg)
        ax2 = fig.add_subplot(412)
        ax2.plot(yday_range, new_site.wind_dir_yday_avg)
        ax3 = fig.add_subplot(413)
        ax3.plot(new_site.date, new_site.wind_speed_anomaly)
        ax4 = fig.add_subplot(414)
        ax4.plot(new_site.date, new_site.wind_dir_anomaly)
        plt.show()

def test_nldas():
    nldas = pickle.load(open("formated_asos_wind/nldas_raws_daily.p", "rb"))
    for site in nldas:
        new_site = Nldas(site.lat, site.lon, site.station_id, site.date, site.wind_speed,
                        site.wind_dir)
        new_site.set_yday()
        new_site.cal_yearly_avg()
        new_site.cal_anomaly()
        print np.where(np.isnan(new_site.wind_speed_yday_avg))
        yday_range = np.unique(new_site.yday)
        yday_range.sort()
        fig = plt.figure()
        ax1 = fig.add_subplot(411)
        ax1.plot(yday_range, new_site.wind_speed_yday_avg)
        ax2 = fig.add_subplot(412)
        ax2.plot(yday_range, new_site.wind_dir_yday_avg)
        ax3 = fig.add_subplot(413)
        ax3.plot(new_site.date, new_site.wind_speed_anomaly)
        ax4 = fig.add_subplot(414)
        ax4.plot(new_site.date, new_site.wind_dir_anomaly)
        plt.show()

def test_raws():
    raws = pickle.load(open("formated_asos_wind/raws_daily.p", "rb"))
    for site in raws:
        new_site = Raws(site.lat, site.lon, site.station_id, site.date, site.wind_speed,
                        site.wind_dir)
        new_site.set_yday()
        new_site.cal_yearly_avg()
        new_site.cal_anomaly()
        yday_range = np.unique(new_site.yday)
        yday_range.sort()
        print np.where(np.isnan(new_site.wind_speed_yday_avg))
        fig = plt.figure()
        ax1 = fig.add_subplot(411)
        ax1.plot(yday_range, new_site.wind_speed_yday_avg)
        ax2 = fig.add_subplot(412)
        ax2.plot(yday_range, new_site.wind_dir_yday_avg)
        ax3 = fig.add_subplot(413)
        ax3.plot(new_site.date, new_site.wind_speed_anomaly)
        ax4 = fig.add_subplot(414)
        ax4.plot(new_site.date, new_site.wind_dir_anomaly)
        plt.show()

def nldas_save():
    nldas_raws = pickle.load(open("formated_asos_wind/nldas_raws_daily.p", "rb"))
    nldas_asos = pickle.load(open("formated_asos_wind/nldas_00_to_15_ts.p", "rb"))
    nldas = nldas_raws + nldas_asos
    for i, site in enumerate(nldas):
        new_site = Nldas(site.lat, site.lon, site.station_id, site.date, site.wind_speed, site.wind_dir)
        new_site.set_yday()
        new_site.cal_yearly_avg()
        new_site.cal_anomaly()
        nldas[i] = new_site
    pickle.dump(nldas, open("wind/nldas.p", "wb"))

def raws_save():
    raws = pickle.load(open("formated_asos_wind/raws_daily.p", "rb"))
    for i, site in enumerate(raws):
        new_site = Raws(site.lat, site.lon, site.station_id, site.date, site.wind_speed,
                        site.wind_dir)
        new_site.set_yday()
        new_site.cal_yearly_avg()
        new_site.cal_anomaly()
        raws[i] = new_site
    pickle.dump(raws, open("wind/raws.p", "wb"))

def asos_save():
    asos = pickle.load(open("formated_asos_wind/asos_info.p", "rb"))
    for i, site in enumerate(asos):
        new_site = Asos(site.lat, site.lon, site.station_id, site.daily_wind.date, site.daily_wind.wind_speed,
                        site.daily_wind.dir)
        new_site.set_yday()
        new_site.cal_yearly_avg()
        new_site.cal_anomaly()
        asos[i] = new_site
    pickle.dump(asos, open("wind/asos.p", "wb"))

def main():
    raws = pickle.load(open("wind/asos.p", "rb"))
    for temp in raws:
        temp.gap_filling()
        temp.ts_plot()

if __name__ == "__main__":
    main()