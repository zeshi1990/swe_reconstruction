__author__ = 'zeshi'

import numpy as np

def snow_surface_temp(current_air_temperature, previous_air_temperature, relative_humidity):
    b = 22.587
    c = 273.86
    below_n_ten_idx = np.where(current_air_temperature <= -10.)
    above_n_ten_idx = np.where(current_air_temperature > -10.)
    snow_surface_temperature = np.zeros(current_air_temperature.shape)
    snow_surface_temperature[above_n_ten_idx] = previous_air_temperature[above_n_ten_idx]
    dew_point_temperature = c * (np.log(relative_humidity) + b * current_air_temperature /
                                 (c + current_air_temperature)) / (b - np.log(relative_humidity) -
                                                                   b * current_air_temperature /
                                                                   (c + current_air_temperature))
    snow_surface_temperature[below_n_ten_idx] = dew_point_temperature[below_n_ten_idx]
    snow_surface_temperature[np.where(snow_surface_temperature > 0.)] = 0.
    return snow_surface_temperature