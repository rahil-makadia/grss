#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <string>
#include <vector>

#include "SpiceUsr.h"

// #define LONGDOUBLE // use long double instead of double
#ifndef real
    #ifdef LONGDOUBLE
        #define real long double
    #else 
        #define real double
    #endif
#endif

struct Constants{
    real du2m=149597870700.0L; // default au to m
    real tu2sec = 86400.0; // default day to sec
    real G = 6.6743e-11L/(du2m*du2m*du2m)*tu2sec*tu2sec; // default kg au^3 / day^2
    real clight = 299792458.0L/du2m*tu2sec; // default au/day
    real j2000Jd = 2451545.0;
    real JdMinusMjd = 2400000.5;
};

struct NBodyParameters{
    size_t nInteg;
    size_t nSpice;
    size_t nTotal;
    size_t dim;
    size_t nOde;
};

struct IntegrationParameters{
    real t0;
    real tf;
    real dt0;
    real dtMax;
    real dtChangeFactor;
    bool adaptiveTimestep;
    real tolPC;
    real tolInteg;
};

struct NongravParams{
    real a1;
    real a2;
    real a3;
    real alpha;
    real k;
    real m;
    real n;
    real r0_au;
};

void get_spice_state_lt(int spiceID, real t0_mjd, Constants consts, double state[6], double &lt);

void jd_to_mjd(const real jd, real &mjd);
real jd_to_mjd(const real jd);
void mjd_to_jd(const real mjd, real &jd);
real mjd_to_jd(const real mjd);
void jd_to_et(const real jd, real &et);
real jd_to_et(const real jd);
void et_to_jd(const real et, real &jd);
real et_to_jd(const real et);
void mjd_to_et(const real mjd, real &et);
real mjd_to_et(const real mjd);
void et_to_mjd(const real et, real &mjd);
real et_to_mjd(const real et);

#endif
