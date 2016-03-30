#!/usr/bin/env python
"""
PAR weather generator
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from math import pi, cos, sin, exp, sqrt, acos, asin
import random

__author__  = "Martin De Kauwe"
__version__ = "1.0 (16.03.2016)"
__email__   = "mdekauwe@gmail.com"

class DiurnalPar(object):

    def __init__(self, lat, lon):

        # conversion from SW (MJ m-2 d-1) to PAR (MJ m-2 d-1)
        self.SW_2_PAR_MJ = 0.5

        # conversion from SW (W m-2) to PAR (umol m-s s-1)
        self.SW_2_PAR = 2.3
        self.PAR_2_SW = 1.0 / self.SW_2_PAR
        self.J_TO_MJ = 1E-6
        self.SEC_2_DAY = 86400.0
        self.MJ_TO_J = 1E6
        self.DAY_2_SEC = 1.0 / self.SEC_2_DAY
        self.J_TO_UMOL = 4.57
        self.UMOL_TO_J = 1.0 / self.J_TO_UMOL
        self.UMOLPERJ = 4.57     # Conversion from J to umol quanta

        self.lat = lat
        self.lon = lon

    def main(self, doy, sw_rad_day):

        # MJ m-2 d-1 -> J m-2 s-1 = W m-2 -> umol m-2 s-1 -> MJ m-2 d-1 #
        #par_day = sw_rad_day * MJ_TO_J * DAY_2_SEC * SW_2_PAR * \
        #          UMOL_TO_J * J_TO_MJ * SEC_2_DAY
        par_day = sw_rad_day * self.SW_2_PAR_MJ

        cos_zenith = self.calculate_solar_geometry(doy)
        diffuse_frac = self.spitters(doy, par_day, cos_zenith)
        par = self.estimate_dirunal_par(par_day, cos_zenith, diffuse_frac)

        par_old = self.estimate_dirunal_par_old(lat, doy, sw_rad_day)

        # check disaggregation
        print sw_rad_day, np.sum(par * self.PAR_2_SW * self.J_TO_MJ * self.SEC_2_DAY), "MJ m-2 d-1"
        print sw_rad_day * self.MJ_TO_J * self.DAY_2_SEC * self.SW_2_PAR, np.sum(par), "umol m-2 s-1"


        plt.plot(np.arange(1,48+1), par, "b-")
        plt.plot(np.arange(1,48+1), par_old, "r-")
        plt.ylabel("par (umol m$^{-2}$ s$^{-1}$)")
        plt.xlabel("Hour of day")
        plt.show()


    def estimate_dirunal_par(self, par_day, cos_zenith, diffuse_frac):
        """
        Calculate daily course of incident PAR from daily totals using routine
        from MAESTRA

        Arguments:
        ----------
        par_day : double
            daily total photosynthetically active radiation (MJ m-2 d-1)
        cos_zenith : float
            cosine of zenith angle (radians)
        diffuse_frac : float
            fraction of par which is diffuse [0-1]
        """
        tau = 0.76 # Transmissivity of atmosphere
        direct_frac = 1.0 - diffuse_frac

        cos_bm = np.zeros(48)
        cos_df = np.zeros(48)
        par = np.zeros(48)

        # daily total beam PAR (MJ m-2 d-1)
        beam_rad = par_day * direct_frac;

        # daily total diffuse PAR (MJ m-2 d-1)
        diffuse_rad = par_day * diffuse_frac;

        sum_bm = 0.0
        sum_df = 0.0
        for i in xrange(1,48+1):

            if cos_zenith[i-1] > 0.0:
                zenith = acos(cos_zenith[i-1])

                # set FBM = 0.0 for ZEN > 80 degrees
                if zenith < 80.0 * pi / 180.0:
                    cos_bm[i-1] = cos_zenith[i-1] * tau**(1.0 / cos_zenith[i-1])
                else:
                    cos_bm[i-1] = 0.0

                cos_df[i-1] = cos_zenith[i-1]
                sum_bm += cos_bm[i-1]
                sum_df += cos_df[i-1]

        for i in xrange(1,48+1):

            if sum_bm > 0.0:
                rdbm = beam_rad * cos_bm[i-1] / sum_bm
            else:
                rdbm = 0.0

            if sum_df > 0.0:
                rddf = diffuse_rad * cos_df[i-1] / sum_df
            else:
                rddf = 0.0

            # MJ m-2 d-1 -> J m-2 s-1 -> umol m-2 s-1
            conv = self.MJ_TO_J * self.DAY_2_SEC * self.UMOLPERJ
            par[i-1] = (rddf + rdbm) * conv

        return par

    def calculate_solar_geometry(self, doy):

        """
        The solar zenith angle is the angle between the zenith and the centre
        of the sun's disc. The solar elevation angle is the altitude of the
        sun, the angle between the horizon and the centre of the sun's disc.
        Since these two angles are complementary, the cosine of either one of
        them equals the sine of the other, i.e. cos theta = sin beta. I will
        use cos_zen throughout code for simplicity.

        Arguments:
        ----------
        doy : double
            day of year
        latitude : float
            latitude (degrees)
        longitude : float
            longitude (degrees)

        References:
        -----------
        * De Pury & Farquhar (1997) PCE, 20, 537-557.
        """

        cos_zenith = np.zeros(48)
        for i in xrange(1, 48+1):
            # need to convert 30 min data, 0-47 to 0-23.5
            hod = i / 2.0

            gamma = self.day_angle(doy)
            dec = self.calculate_solar_declination(doy, gamma)
            et = self.calculate_eqn_of_time(gamma)
            t0 = self.calculate_solar_noon(et)
            h = self.calculate_hour_angle(hod, t0)
            rlat = self.lat * pi / 180.0

            # A13 - De Pury & Farquhar
            sin_beta = sin(rlat) * sin(dec) + cos(rlat) * cos(dec) * cos(h)
            cos_zenith[i-1] = sin_beta; # The same thing, going to use throughout
            if cos_zenith[i-1] > 1.0:
                cos_zenith[i-1] = 1.0
            elif cos_zenith[i-1] < 0.0:
                cos_zenith[i-1] = 0.0

            zenith = 180.0 / pi * acos(cos_zenith[i-1])
            elevation = 90.0 - zenith

        return cos_zenith

    def calculate_solar_noon(self, et):
        """
        Calculation solar noon - De Pury & Farquhar, '97: eqn A16

        Arguments:
        ----------
        et : float
            equation of time (radians)
        longitude : float
            longitude (degrees)

        Returns:
        ---------
        t0 - solar noon (hours).

        Reference:
        ----------
        * De Pury & Farquhar (1997) PCE, 20, 537-557.
        """
        # all international standard meridians are multiples of 15deg east/west of
        # greenwich
        Ls = self.round_to_value(self.lon, 15.)
        t0 = 12.0 + (4.0 * (Ls - self.lon) - et) / 60.0

        return (t0)

    def round_to_value(self, number, roundto):
        return round(number / roundto) * roundto

    def calculate_solar_declination(self, doy, gamma):
        """
        Solar Declination Angle is a function of day of year and is indepenent
        of location, varying between 23deg45' to -23deg45'

        Arguments:
        ----------
        doy : int
            day of year, 1=jan 1
        gamma : double
            fractional year (radians)

        Returns:
        --------
        dec: float
            Solar Declination Angle [radians]

        Reference:
        ----------
        * De Pury & Farquhar (1997) PCE, 20, 537-557.
        * Leuning et al (1995) Plant, Cell and Environment, 18, 1183-1200.
        * J. W. Spencer (1971). Fourier series representation of the position of
          the sun.
        """

        # declination (radians) */
        #decl = 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) - \
        #       0.006758 * cos(2.0 * gamma) + 0.000907 * sin(2.0 * gamma) -\
        #       0.002697 * cos(3.0 * gamma) + 0.00148 * sin(3.0 * gamma);*/

        # (radians) A14 - De Pury & Farquhar
        decl = -23.4 * (pi / 180.) * cos(2.0 * pi * (doy + 10) / 365)

        return (decl)

    def calculate_hour_angle(self, t, t0):
        """
        Calculation solar noon - De Pury & Farquhar, '97: eqn A15
        Reference:
        ----------
        * De Pury & Farquhar (1997) PCE, 20, 537-557.
        Returns:
        ---------
        h - hour angle (radians).
        """
        return (pi * (t - t0) / 12.0)


    def day_angle(self, doy):
        """
        Calculation of day angle - De Pury & Farquhar, '97: eqn A18

        Arguments:
        ----------
        doy : int
            day of year

        Reference:
        ----------
        * De Pury & Farquhar (1997) PCE, 20, 537-557.
        * J. W. Spencer (1971). Fourier series representation of the position of
          the sun.

        Returns:
        ---------
        gamma - day angle in radians.
        """

        return (2.0 * pi * (float(doy) - 1.0) / 365.0)

    def calculate_eqn_of_time(self, gamma):
        """
        Equation of time - correction for the difference btw solar time
        and the clock time.

        Arguments:
        ----------
        gamma : double
            fractional year (radians)

        References:
        -----------
        * De Pury & Farquhar (1997) PCE, 20, 537-557.
        * Campbell, G. S. and Norman, J. M. (1998) Introduction to environmental
          biophysics. Pg 169.
        * J. W. Spencer (1971). Fourier series representation of the position of
          the sun.
        * Hughes, David W.; Yallop, B. D.; Hohenkerk, C. Y. (1989),
          "The Equation of Time", Monthly Notices of the Royal Astronomical
          Society 238: 1529-1535
        """

        # minutes - de Pury and Farquhar, 1997 - A17
        et = (0.017 + 0.4281 * cos(gamma) - 7.351 * sin(gamma) - 3.349 *
              cos(2.0 * gamma) - 9.731  * sin(gamma))

        return (et)

    def calc_extra_terrestrial_rad(self, doy, cos_zen):
        """
        Solar radiation incident outside the earth's atmosphere, e.g.
        extra-terrestrial radiation. The value varies a little with the earths
        orbit.
        Using formula from Spitters not Leuning!

        Arguments:
        ----------
        doy : double
            day of year
        cos_zen : double
            cosine of zenith angle (radians)

        Returns:
        --------
        So : float
            solar radiation normal to the sun's bean outside the Earth's atmosphere
            (J m-2 s-1)

        Reference:
        ----------
        * Spitters et al. (1986) AFM, 38, 217-229, equation 1.
        """

        # Solar constant (J m-2 s-1)
        Sc = 1370.0

        if cos_zen > 0.0:
            # remember sin_beta = cos_zenith; trig funcs are cofuncs of each other
            # sin(x) = cos(90-x) and cos(x) = sin(90-x).
            So = Sc * (1.0 + 0.033 * cos(float(doy) / 365.0 * 2.0 * pi)) * cos_zen
        else:
            So = 0.0

        return So

    def spitters(self, doy, par, cos_zenith):
        """
        Spitters algorithm to estimate the diffuse component from the total daily
        incident radiation.

        NB. Eqns. 2a-d, not 20a-d

        Parameters:
        ----------
        doy : int
            day of year
        par : double
            daily total photosynthetically active radiation (MJ m-2 d-1)
        cos_zenith : float
            cosine of zenith angle (radians)

        Returns:
        -------
        diffuse : double
            diffuse component of incoming radiation

        References:
        ----------
        * Spitters, C. J. T., Toussaint, H. A. J. M. and Goudriaan, J. (1986)
          Separating the diffuse and direct component of global radiation and its
          implications for modeling canopy photosynthesis. Part I. Components of
          incoming radiation. Agricultural Forest Meteorol., 38:217-229.
        """

        # Fraction of global radiation that is PAR
        fpar = 0.5
        SEC_2_HFHR = 1800.0
        J_TO_MJ = 1E-6
        CONV = SEC_2_HFHR * J_TO_MJ

        # Calculate extra-terrestrial radiation
        S0 = 0.0
        for i in xrange(48):
            S0 += self.calc_extra_terrestrial_rad(doy, cos_zenith[i]) * CONV

        # atmospheric transmisivity
        tau = (par / fpar) / S0

        ## Spitter's formula (Eqns. 2a-d)
        if tau < 0.07:
            diffuse_frac = 1.0
        elif tau < 0.35:
            diffuse_frac = 1.0 - 2.3 * (tau - 0.07)**2
        elif tau < 0.75:
            diffuse_frac = 1.33 - 1.46 * tau
        else:
            diffuse_frac = 0.23

        return (diffuse_frac)


    def estimate_dirunal_par_old(self, lat, doy, sw_rad_day):

        # Solar constant [MJ/m2/day]
        solar_constant = 1370.0 * self.SEC_2_DAY / 1E6

        # convert to radians
        rlat  = lat * np.pi / 180.0
        ryear = 2.0 * np.pi * (doy - 1.0) / 365.0

        # Declination in radians (Paltridge and Platt eq [3.7])
        rdec = (0.006918 - 0.399912 * np.cos(ryear) + 0.070257 *
                np.sin(ryear) - 0.006758 * np.cos(2.0 * ryear) +
                0.000907 * np.sin(2.0 * ryear) - 0.002697 *
                np.cos(3.0  *ryear) + 0.001480 * np.sin(3.0 * ryear))

        # Declination in degrees
        dec = (180.0 / np.pi) * rdec

        hour_angle = -np.tan(rlat) * np.tan(rdec)

        # Half Day Length (radians), (dawn:noon = noon:dusk)
        # Paltridge and Platt eq [3.21]
        if hour_angle <= -1.0:
            rhlf_day_length = np.pi            # polar summer: sun never sets
        elif hour_angle >= 1.0:
            rhlf_day_length = 0.0              # polar winter: sun never rises
        else:
            rhlf_day_length = np.arccos(hour_angle)

        # hrs
        day_length = 24.0 * 2.0 * rhlf_day_length / (2.0 * np.pi)

        # Daily solar irradiance without atmosphere, normalised by solar constant,
        # with both energy fluxes in MJ/m2/day, calculated from solar geometry
        # Paltridge and Platt eq [3.22]
        solar_norm = (rhlf_day_length * np.sin(rlat) * np.sin(rdec) +
                      np.cos(rlat) * np.cos(rdec) * np.sin(rhlf_day_length)) / np.pi

        # Observed daily solar irradiance as frac of value from solar geometry
        # (should be < 1), sw_rad (MJ m-2 d-1)
        solar_frac_obs = sw_rad_day / (solar_norm * solar_constant + 1.0E-6)

        sw_rad = np.zeros(48)

        # disaggregate radiation
        for i in xrange(1,48+1):

            hour = i / 2.0

            # Time in day frac (-0.5 to 0.5, zero at noon)
            time_noon = float(hour) / 24.0 - 0.5

            # Time in day frac (-Pi to Pi, zero at noon)
            rtime = 2.0 * np.pi * time_noon

            #print hour, time_noon, rtime
            if (np.abs(rtime) < rhlf_day_length): # day: sun is up
                # Paltridge and Platt eq [3.4]
                sw_rad[i-1] = ((sw_rad_day / solar_norm) * (np.sin(rdec) *
                              np.sin(rlat) + np.cos(rdec) * np.cos(rlat) *
                              np.cos(rtime)))

                # Need to half rad, as we are doing this over 48 timesteps and not
                # 24, otherwise we will end up with a total which is double
                # the sw_rad_tot
                #sw_rad[i-1] /= 2.0
            else:
                sw_rad[i-1] = 0.0

        # check the integration worked?
        difference = sw_rad_day - np.sum(sw_rad / 48.)
        print "*", difference, sw_rad_day, np.sum(sw_rad / 48.)

        # Convert sw_rad (MJ/m2/day to J/m2/s = W/m2) to PAR (umol m-2 s-1)
        par = sw_rad * self.MJ_TO_J * self.DAY_2_SEC * self.SW_2_PAR

        return (par)

if __name__ == "__main__":

    lat = -23.575001
    lon = 152.524994
    doy = 180.0
    sw_rad_day = 12.5 # mj m-2 d-1

    P = DiurnalPar(lat, lon)
    P.main(doy, sw_rad_day)
