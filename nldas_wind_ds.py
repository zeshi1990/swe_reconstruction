__author__ = 'zeshi'

import numpy as np
import pylab as plt

class wind_redistribution():

    def __init__(self, res=500):
        self.slope = np.load("DEM/slope_" + str(res) + "m.npy")
        self.aspect = np.load("DEM/aspect_" + str(res) + "m.npy")
        self.curvature = np.load("DEM/curvature_" + str(res) + "m.npy")
        self.forest_correction = np.load("vegetation/wind_correction_" + str(res) + "m.npy")
        self.res = res

    def wind_redistribution(self, u_wind, v_wind):
        wind_magnitude = np.sqrt(u_wind ** 2 + v_wind ** 2)
        theta = 3 * np.pi / 2 - np.arctan2(v_wind, u_wind)
        omega_s = self.slope * np.cos(theta - self.aspect)
        omega_c = self.curvature
        wind_correction_weight =  1. + 0.5 * omega_s + 0.5 * omega_c
        wind_corrected_magnitude = wind_correction_weight * wind_magnitude
        # wind_corrected_magnitude = wind_correction_weight * self.forest_correction * wind_magnitude
        # self.wind_compare(wind_magnitude, wind_corrected_magnitude)
        return wind_corrected_magnitude

    def wind_compare(self, original, corrected):
        min = np.min([np.min(original), np.min(corrected)])
        max = np.max([np.max(original), np.max(corrected)])
        plt.subplot(2, 1, 1)
        plt.imshow(original, vmin=min, vmax=max)
        plt.subplot(2, 1, 2)
        plt.imshow(corrected, vmin=min, vmax=max)
        plt.colorbar(orientation="horizontal")
        plt.show()