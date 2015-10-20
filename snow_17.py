__author__ = 'zeshi'

import numpy as np
import gdal
from gdal import GA_ReadOnly

def calc_rh(air_tmp, pres, spfh):
    rh = 0.263 * pres * spfh * (1 / np.exp(17.67 * air_tmp / (air_tmp + 243.5)))
    return rh

class snow17():
    """
    This class using the snow-17 model and liston et al 1999 to calculate the
    latent heat and sensible heat from the surface, it is a class, user has to
    insert temperature, pressure, specific humidity and wind and resolution (30 or 500)
    to initiate the class and it will calculate LH and SH automatically
    """
    SURFACE_SATURATE_VAPOR_PRESSURE = 611                                   # Pa
    Z_A = 1000.0                                                            # cm, z_a in snow-17 model
    Z_0 = 0.01                                                              # cm, z_0 in snow-17 model
    L_S = 677.0                                                             # cal / g, latent heat of sublimation
    RHO_W = 1.0                                                             # g / cm^3, density of water
    R_SPECIFIC = 287.058                                                    # J/(kg * K), specific gas constant for dry air
    K = 0.40                                                                # dimensionless, von Karman's constant
    DEM_30 = gdal.Open("DEM/30m_dem.tif", GA_ReadOnly).ReadAsArray()       # 30m DEM
    DEM_500 = gdal.Open("DEM/500m_dem.tif", GA_ReadOnly).ReadAsArray()       # 500m DEM
    P_A_30 = 33.86 * (29.9 - 0.335 * DEM_30 + 0.00022 * DEM_30 ** 2.4)      # 30m standard atmosphere
    P_A_500 = 33.86 * (29.9 - 0.335 * DEM_500 + 0.00022 * DEM_500 ** 2.4)   # 500m standard atmosphere
    PA_2_MB = 0.01                                                          # 0.01 Pa/mb

    def __init__(self, air_tmp, snow_tmp, pres, spfh, wind, dlw, precip, res=500):
        self.air_tmp = air_tmp
        self.snow_tmp = snow_tmp
        self.pres = pres
        self.spfh = spfh
        self.wind = wind
        self.dlw = dlw
        self.precip = precip
        afd_fn = "vegetation/amr_nlcd2011_" + str(res) + "m.tif"
        self.F_c = gdal.Open(afd_fn, GA_ReadOnly).ReadAsArray().astype(float) / 100.0
        self.res = res
        self.calc_relative_humidity()
        self.calc_air_saturate_vapor_pressure()
        self.calc_snow_saturate_vapor_pressure()
        self.calc_air_vapor_pressure()
        self.calc_air_density()
        # self.calc_b()
        # self.calc_FU()
        # self.calc_exchange_coefficient()
        # self.calc_unstable_correction()
        self.calc_atm_emissivity()
        self.calc_weighted_emissivity()
        self.calc_LH()
        self.calc_SH()
        self.calc_DLW()
        self.calc_ULW()
        self.calc_Qm()

    def calc_relative_humidity(self):
        """
        Calculate relative humidity
        :return:
        """
        self.rh = 0.263 * self.pres * self.spfh * (1 / np.exp(17.67 * self.air_tmp / (self.air_tmp + 243.5)))


    def calc_air_saturate_vapor_pressure(self):
        """
        Calculate air saturate vapor pressure
        :return:
        """
        self.e_sat_air = 611.2 * np.exp(17.67 * self.air_tmp / (243.5 + self.air_tmp))

    def calc_snow_saturate_vapor_pressure(self):
        """
        Calculate snow saturate vapor pressure
        :return:
        """
        self.e_sat_snow = 611.2 * np.exp(17.67 * self.snow_tmp / (243.5 + self.snow_tmp))

    def calc_air_vapor_pressure(self):
        """
        Calculate air vapor pressure
        :return:
        """
        self.e_a = (self.rh / 100.0) * self.e_sat_air

    def calc_air_density(self):
        """
        Calculate dry air density
        :return:
        """
        self.rho_a = self.pres / (self.R_SPECIFIC * (self.air_tmp + 273.16))

    def calc_b(self):
        """
        Calculate b parameter in wind corrections
        :return:
        """
        if self.res == 500:
            self.b = (0.622 * self.rho_a) / (self.P_A_500 * self.RHO_W) * 10. ** 6. \
                     * (self.K ** 2. / np.log(self.Z_A / self.Z_0) ** 2.)
        else:
            self.b = (0.622 * self.rho_a) / (self.P_A_30 * self.RHO_W) * 10. ** 6. \
                     * (self.K ** 2. / np.log(self.Z_A / self.Z_0) ** 2.)

    def calc_FU(self):
        """
        Calculate wind function
        :return:
        """
        self.fu = self.b * self.wind * 3.6

    def calc_exchange_coefficient(self):
        """
        Calculate exchange coefficient D_(h,e) in Liston 1999
        :return:
        """
        self.d_h_e = self.K ** 2. * self.wind / np.log(self.Z_A / self.Z_0) ** 2.

    def calc_unstable_correction(self):
        """
        Calculate unstable correction for latent and sensible heat
        :return:
        """
        gamma = 5.3 * 9.4 * self.d_h_e / self.wind * (self.Z_A / self.Z_0) ** (0.5)
        R_i = 2. * 10. * self.Z_A / 100. * self.air_tmp / (self.air_tmp * (27.78 * self.wind) ** 2.)
        zeta = np.zeros(R_i.shape)
        zeta[np.where(R_i > 0)] = 1. / (1. + 4.7 * R_i[np.where(R_i > 0)]) ** 2.
        zeta[np.where(R_i < 0)] = 1. - (9.4 * R_i[np.where(R_i < 0)] / (1. + gamma[np.where(R_i < 0)] *
                                                                       np.abs(R_i[np.where(R_i < 0)]) ** (0.5)))
        self.zeta = zeta

    def calc_atm_emissivity(self):
        """
        Calculate atmospheric emissivity
        :return:
        """
        self.atm_emissivity = self.dlw / (((self.air_tmp + 273.16) ** 4.0) * (5.67 * 10.0 ** (-8)))

    def calc_weighted_emissivity(self):
        """
        Calculate weighted emissivity
        :return:
        """
        self.weighted_emissivity = ((1.0 - self.F_c) * self.atm_emissivity) + self.F_c * 0.98

    def calc_LH(self):
        """
        Calculate latent heat
        :return:
        """
        self.LH = (2500.8 - 2.366 * self.snow_tmp) * self.rho_a * 0.0029 * (62.2 / self.pres) * \
                  ((self.e_a - self.e_sat_snow)) / 100. * self.wind * 1000

    def calc_SH(self):
        """
        Calculate sensible heat
        :return:
        """
        self.SH = self.rho_a * 1005 * 0.0029 * (self.air_tmp - self.snow_tmp) * self.wind

    def calc_DLW(self):
        """
        Calculate canopy corrected down longwave radiation
        :return:
        """
        self.CDLW = self.dlw * self.weighted_emissivity / self.atm_emissivity

    def calc_ULW(self):
        """
        Calculate upwelling longwave radiation
        :return:
        """
        self.ULW = 0.98 * 5.67 * 10.0 ** (-8) * (self.snow_tmp + 273.16) ** 4

    def calc_Qm(self):
        """
        Calculate the heat transfer from the precipitation, only rain is counted
        :return:
        """
        precip = self.precip
        air_tmp = self.air_tmp
        self.Qm = precip / 3600.0 * 4180.0 * air_tmp
