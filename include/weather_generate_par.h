#ifndef WEATHER_GENERATE_PAR_H
#define WEATHER_GENERATE_PAR_H

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <sys/types.h>
#include <dirent.h>
#include <string.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#define STRING_LENGTH 2000
#define TRUE 1
#define FALSE 0

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

#define SW_2_PAR 2.3
#define MJ_TO_J 1E6
#define SEC_2_DAY 86400.0
#define DAY_2_SEC 1.0 / SEC_2_DAY
#define J_TO_UMOL 4.57
#define UMOL_TO_J 1.0 / J_TO_UMOL
#define J_TO_MJ 1E-6
#define hPa_2_kPa 0.1
#define DEG_TO_KELVIN 273.15
#define SEC_2_HFHR 1800.0
#define NTIMESTEPS 48
#define SW_2_PAR_MJ 0.5 /* conversion from SW MJ m-2 d-1 to PAR MJ m-2 d-1 */
#define HLFHR_2_SEC 1.0 / 1800.0

void   estimate_dirunal_par(float, float, int, float, float *);
float  calc_day_length(int, int, float);
float  calc_vpd(float, float);
int    is_leap_year(int);
float  spitters(int, float, float *);
float  day_angle(int);
float  calculate_solar_declination(int, float);
float  calculate_eqn_of_time(float);
float  calculate_solar_noon(float, float);
float  calculate_hour_angle(float, float);
float  calc_extra_terrestrial_rad(int, float);
float  round_to_value(float, float);
void   calculate_solar_geometry(int, float, float, float *);


#endif
