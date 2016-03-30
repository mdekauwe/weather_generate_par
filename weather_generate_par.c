/* Quick test of weather generator code in C to find translation issue */
#include "weather_generate_par.h"

int main(int argc, char **argv) {

    float par[NTIMESTEPS];
    float latitude = -23.575001;
    float longitude = 152.524994;
    int i, j;

    /* MJ m-2 d-1 */
    float sw = 12.5;

    for (i = 0; i < 365; i++) {

        estimate_dirunal_par(latitude, longitude, i+1, sw, &(par[0]));
        for (j = 0; j < NTIMESTEPS; j++) {
            printf("%d %f\n", j+1, par[j]);
        }
    }

    return (EXIT_SUCCESS);
}

void estimate_dirunal_par(float lat, float lon, int doy, float sw_rad_day,
                          float *par) {
    /*
        Calculate daily course of incident PAR from daily totals using routine
        from MAESTRA
    */
    int   i;
    float cos_zenith[NTIMESTEPS];
    float tau = 0.76;            /* Transmissivity of atmosphere */
    float direct_frac, diffuse_frac;
    float cos_bm[NTIMESTEPS], cos_df[NTIMESTEPS], sum_bm, sum_df;
    float zenith, rddf, rdbm, par_day, beam_rad, diffuse_rad;

    /* MJ m-2 d-1 -> J m-2 s-1 = W m-2 -> umol m-2 s-1 -> MJ m-2 d-1 */
    par_day = sw_rad_day * MJ_TO_J * DAY_2_SEC * SW_2_PAR * \
              UMOL_TO_J * J_TO_MJ * SEC_2_DAY;

    calculate_solar_geometry(doy, lat, lon, &(cos_zenith[0]));
    diffuse_frac = spitters(doy, par_day, cos_zenith);
    direct_frac = 1.0 - diffuse_frac;

    /* daily total beam PAR (MJ m-2 d-1) */
    beam_rad = par_day * direct_frac;

    /* daily total diffuse PAR (MJ m-2 d-1) */
    diffuse_rad = par_day * diffuse_frac;

    sum_bm = 0.0;
    sum_df = 0.0;
    for (i = 0; i < NTIMESTEPS; i++) {
        cos_bm[i] = 0.0;
        cos_df[i] = 0.0;

        if (cos_zenith[i] > 0.0) {
            zenith = acos(cos_zenith[i]);

            /* set FBM = 0.0 for ZEN > 80 degrees */
            if (zenith < (80.0 * M_PI / 180.0)) {
                cos_bm[i] = cos_zenith[i] * pow(tau, (1.0 / cos_zenith[i]));
            } else {
                cos_bm[i] = 0.0;
            }
            cos_df[i] = cos_zenith[i];
            sum_bm += cos_bm[i];
            sum_df += cos_df[i];
        }
    }

    for (i = 0; i < NTIMESTEPS; i++) {

        if (sum_bm > 0.0) {
            rdbm = beam_rad * cos_bm[i] / sum_bm;
        } else {
            rdbm = 0.0;
        }

        if (sum_df > 0.0) {
            rddf = diffuse_rad * cos_df[i] / sum_df;
        } else {
            rddf = 0.0;
        }

        /* MJ m-2 d-1 -> J m-2 s-1 -> umol m-2 s-1 */
        *(par+i) = (rddf + rdbm) * MJ_TO_J * J_TO_UMOL * DAY_2_SEC;
    }

    return;
}

void calculate_solar_geometry(int doy, float latitude, float longitude,
                              float *cos_zenith) {
    /*
    The solar zenith angle is the angle between the zenith and the centre
    of the sun's disc. The solar elevation angle is the altitude of the
    sun, the angle between the horizon and the centre of the sun's disc.
    Since these two angles are complementary, the cosine of either one of
    them equals the sine of the other, i.e. cos theta = sin beta. I will
    use cos_zen throughout code for simplicity.

    Arguments:
    ----------
    doy : float
        day of year
    latitude : float
        latitude (degrees)
    longitude : float
        longitude (degrees)

    References:
    -----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.
    */
    int   i;
    float rdec, et, t0, h, gamma, rlat, sin_beta;
    float hod;

    for (i = 1; i < NTIMESTEPS+1; i++) {

        /* need to convert 30 min data, 0-47 to 0-23.5 */
        hod = i / 2.0;

        gamma = day_angle(doy);
        rdec = calculate_solar_declination(doy, gamma);
        et = calculate_eqn_of_time(gamma);
        t0 = calculate_solar_noon(et, longitude);
        h = calculate_hour_angle(hod, t0);
        rlat = latitude * M_PI / 180.0;

        /* A13 - De Pury & Farquhar */
        sin_beta = sin(rlat) * sin(rdec) + cos(rlat) * cos(rdec) * cos(h);
        /* The same thing, going to use throughout */
        *(cos_zenith+(i-1)) = sin_beta;
        if (*(cos_zenith+(i-1)) > 1.0) {
            *(cos_zenith+(i-1)) = 1.0;
        } else if (cos_zenith[i-1] < 0.0) {
            *(cos_zenith+(i-1)) = 0.0;
        }
        /*zenith = 180.0 / M_PI * acos(cos_zenith[i-1]);
        elevation = 90.0 - zenith;*/
    }
    return;
}

float day_angle(int doy) {
    /* Calculation of day angle - De Pury & Farquhar, '97: eqn A18

    Reference:
    ----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.
    * J. W. Spencer (1971). Fourier series representation of the position of
      the sun.

    Returns:
    ---------
    gamma - day angle in radians.
    */

    return (2.0 * M_PI * ((float)doy - 1.0) / 365.0);
}

float calculate_solar_declination(int doy, float gamma) {
    /*
    Solar Declination Angle is a function of day of year and is indepenent
    of location, varying between 23deg45' to -23deg45'

    Arguments:
    ----------
    doy : int
        day of year, 1=jan 1
    gamma : float
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
    */
    float decl;

    /* declination (radians) */
    /*decl = 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) - \
           0.006758 * cos(2.0 * gamma) + 0.000907 * sin(2.0 * gamma) -\
           0.002697 * cos(3.0 * gamma) + 0.00148 * sin(3.0 * gamma);*/


    /* (radians) A14 - De Pury & Farquhar  */
    decl = -23.4 * (M_PI / 180.) * cos(2.0 * M_PI * ((float)doy + 10.) / 365.);

    return (decl);

}

float calculate_eqn_of_time(float gamma) {
    /* Equation of time - correction for the difference btw solar time
    and the clock time.

    Arguments:
    ----------
    doy : int
        day of year
    gamma : float
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
      Society 238: 1529–1535
    */
    float et;

    /* radians */
    et = 0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) -\
         0.014615 * cos(2.0 * gamma) - 0.04089 * sin(2.0 * gamma);

    /* radians to minutes */
    et *= 229.18;

    /* radians to hours */
    /*et *= 24.0 / (2.0 * M_PI);*/

    /* minutes - de Pury and Farquhar, 1997 - A17 */
    /*et = (0.017 + 0.4281 * cos(gamma) - 7.351 * sin(gamma) - 3.349 *
          cos(2.0 * gamma) - 9.731  * sin(gamma));*/

    return (et);
}

float calculate_solar_noon(float et, float longitude) {
    /* Calculation solar noon - De Pury & Farquhar, '97: eqn A16

    Reference:
    ----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.

    Returns:
    ---------
    t0 - solar noon (hours).
    */
    float t0, Ls;

    /* all international standard meridians are multiples of 15deg east/west of
       greenwich */
    Ls = round_to_value(longitude, 15.);
    t0 = 12.0 + (4.0 * (Ls - longitude) - et) / 60.0;

    return (t0);
}

float calculate_hour_angle(float t, float t0) {
    /* Calculation solar noon - De Pury & Farquhar, '97: eqn A15

    Reference:
    ----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.

    Returns:
    ---------
    h - hour angle (radians).
    */
    return (M_PI * (t - t0) / 12.0);

}


float spitters(int doy, float par, float *cos_zenith) {
    /*
    Spitters algorithm to estimate the diffuse component from the total daily
    incident radiation.

    NB. Eqns. 2a-d, not 20a-d

    Parameters:
    ----------
    doy : int
        day of year
    par : float
        daily total photosynthetically active radiation (MJ m-2 d-1)
    cos_zenith : float
        cosine of zenith angle (radians)

    Returns:
    -------
    diffuse : float
        diffuse component of incoming radiation

    References:
    ----------
    * Spitters, C. J. T., Toussaint, H. A. J. M. and Goudriaan, J. (1986)
      Separating the diffuse and direct component of global radiation and its
      implications for modeling canopy photosynthesis. Part I. Components of
      incoming radiation. Agricultural Forest Meteorol., 38:217-229.
    */

    /* Fraction of global radiation that is PAR */
    float fpar = 0.5;
    float conv = SEC_2_HFHR * J_TO_MJ;
    float S0, tau, diffuse_frac;
    int   i;


    /* Calculate extra-terrestrial radiation */
    S0 = 0.0;
    for (i = 1; i < NTIMESTEPS+1; i++) {
        S0 += calc_extra_terrestrial_rad(doy, *(cos_zenith+(i-1))) * conv;
    }

    /* atmospheric transmisivity */
    tau = (par / fpar) / S0;

    /* Spitter's formula (Eqns. 2a-d) */
    if (tau < 0.07) {
        diffuse_frac = 1.0;
    } else if (tau < 0.35) {
        diffuse_frac = 1.0 - 2.3 * (tau - 0.07) * (tau - 0.07);
    } else if (tau < 0.75) {
        diffuse_frac = 1.33 - 1.46 * tau;
    } else {
        diffuse_frac = 0.23;
    }

    return (diffuse_frac);
}

float calc_extra_terrestrial_rad(int doy, float cos_zenith) {
    /* Solar radiation incident outside the earth's atmosphere, e.g.
    extra-terrestrial radiation. The value varies a little with the earths
    orbit.

    Using formula from Spitters not Leuning!

    Arguments:
    ----------
    doy : double
        day of year
    cos_zenith : double
        cosine of zenith angle (radians)

    Returns:
    --------
    So : float
        solar radiation normal to the sun's bean outside the Earth's atmosphere
        (J m-2 s-1)

    Reference:
    ----------
    * Spitters et al. (1986) AFM, 38, 217-229, equation 1.
    */

    float So, Sc;

    /* Solar constant (J m-2 s-1) */
    Sc = 1370.0;

    if (cos_zenith > 0.0) {
        /*
        ** remember sin_beta = cos_zenith; trig funcs are cofuncs of each other
        ** sin(x) = cos(90-x) and cos(x) = sin(90-x).
        */
        So = Sc * (1.0 + 0.033 * cos((float)doy / 365.0 * 2.0 * M_PI)) *\
                cos_zenith;
    } else {
        So = 0.0;
    }

    return (So);

}

float round_to_value(float number, float roundto) {
    return (round(number / roundto) * roundto);
}
