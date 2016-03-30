/* Quick test of weather generator code in C to find translation issue */
#include "weather_generate_par.h"

int main(int argc, char **argv) {

    float par[NTIMESTEPS];
    float latitude = -23.575001;
    float longitude = 152.524994;
    int i;

    float sw = ;

    for (i = 0; i < 365; i++) {

        estimate_dirunal_par(latitude, longitude, i+1, sw, &(par[0]));

    }

    return (EXIT_SUCCESS);
}
