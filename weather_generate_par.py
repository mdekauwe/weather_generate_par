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


def main():

    # conversion from SW (MJ m-2 d-1) to PAR (MJ m-2 d-1)
    SW_2_PAR_MJ = 0.5

    # conversion from SW (W m-2) to PAR (umol m-s s-1)
    SW_2_PAR = 2.3
    PAR_2_SW = 1.0 / SW_2_PAR
    J_TO_MJ = 1E-6
    SEC_2_DAY = 86400.0
    MJ_TO_J = 1E6
    DAY_2_SEC = 1.0 / SEC_2_DAY
    J_TO_UMOL = 4.57
    UMOL_TO_J = 1.0 / J_TO_UMOL

    lat = -23.575001
    lon = 152.524994
    doy = 180.0
    sw_rad_day = 12.5 # mj m-2 d-1

    # MJ m-2 d-1 -> J m-2 s-1 = W m-2 -> umol m-2 s-1 -> MJ m-2 d-1 #
    #par_day = sw_rad_day * MJ_TO_J * DAY_2_SEC * SW_2_PAR * \
    #          UMOL_TO_J * J_TO_MJ * SEC_2_DAY
    par_day = sw_rad_day * SW_2_PAR_MJ

    cos_zenith = calculate_solar_geometry(doy, lat, lon)
    diffuse_frac = spitters(doy, par_day, cos_zenith)
    par = estimate_dirunal_par(par_day, cos_zenith, diffuse_frac)

    # check disaggregation
    print sw_rad_day, np.sum(par * PAR_2_SW * J_TO_MJ * SEC_2_DAY), "MJ m-2 d-1"
    print sw_rad_day * MJ_TO_J * DAY_2_SEC * SW_2_PAR, np.sum(par), "umol m-2 s-1"


    plt.plot(np.arange(1,48+1), par, "b-")
    plt.ylabel("par (umol m$^{-2}$ s$^{-1}$)")
    plt.xlabel("Hour of day")
    plt.show()


def estimate_dirunal_par(par_day, cos_zenith, diffuse_frac):
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

    # Transmissivity of atmosphere
    tau = 0.76
    direct_frac = 1.0 - diffuse_frac
    MJ_TO_J = 1E6
    UMOLPERJ = 4.57     # Conversion from J to umol quanta
    SEC_2_DAY = 86400.0
    DAY_2_SEC = 1.0 / SEC_2_DAY

    cos_bm = np.zeros(48)
    cos_df = np.zeros(48)

    sum_bm = 0.0
    sum_df = 0.0
    for i in xrange(1,48+1):

        hrtime = float(i) - 0.5

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

    par = np.zeros(48)
    for i in xrange(1,48+1):

        if sum_bm > 0.0:
            rdbm = par_day * direct_frac * cos_bm[i-1] / sum_bm
        else:
            rdbm = 0.0

        if sum_df > 0.0:
            rddf = par_day * diffuse_frac * cos_df[i-1] / sum_df
        else:
            rddf = 0.0

        # MJ m-2 d-1 -> J m-2 s-1 -> umol m-2 s-1
        par[i-1] = (rddf + rdbm) * MJ_TO_J * DAY_2_SEC * UMOLPERJ

    return par

def calculate_solar_geometry(doy, latitude, longitude):

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

        gamma = day_angle(doy)
        dec = calculate_solar_declination(doy, gamma)
        et = calculate_eqn_of_time(gamma)
        t0 = calculate_solar_noon(et, longitude)
        h = calculate_hour_angle(hod, t0)
        rlat = latitude * pi / 180.0

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

def calculate_solar_noon(et, longitude):
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
    Ls = round_to_value(longitude, 15.)
    t0 = 12.0 + (4.0 * (Ls - longitude) - et) / 60.0

    return (t0)

def round_to_value(number, roundto):
    return round(number / roundto) * roundto

def calculate_solar_declination(doy, gamma):
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

def calculate_hour_angle(t, t0):
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


def day_angle(doy):
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

def calculate_eqn_of_time(gamma):
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

def calc_extra_terrestrial_rad(doy, cos_zen):
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

def spitters(doy, par, cos_zenith):
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
        S0 += calc_extra_terrestrial_rad(doy, cos_zenith[i]) * CONV

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



if __name__ == "__main__":

    main()
