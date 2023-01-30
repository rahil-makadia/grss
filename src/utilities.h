#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <cmath>
using std::fmin;
using std::fmax;
using std::fabs;
using std::pow;
using std::sqrt;
using std::sin;
using std::cos;
using std::tan;
using std::asin;
using std::acos;
using std::atan2;

#include "SpiceUsr.h"

#ifndef LONGDOUBLE
    #define LONGDOUBLE  // use long double instead of double
#endif

#ifndef real
    #ifdef LONGDOUBLE
        #define real long double
    #else 
        #define real double
    #endif
#endif
#define PI 3.141592653589793238462643383279502884197169399375105820974944
#define RAD2DEG 180.0/PI
#define DEG2RAD PI/180.0
#define EARTH_OBLIQUITY 84381.448/3600.0*DEG2RAD

struct Constants{
    real du2m=149597870700.0L; // default au to m
    real tu2sec = 86400.0; // default day to sec
    real G = 6.6743e-11L/(du2m*du2m*du2m)*tu2sec*tu2sec; // default kg au^3 / day^2
    real clight = 299792458.0L/du2m*tu2sec; // default au/day
    real j2000Jd = 2451545.0;
    real JdMinusMjd = 2400000.5;
};

struct IntegrationParameters{
    size_t nInteg;
    size_t nSpice;
    size_t nTotal;
    real t0;
    real tf;
    real dt0;
    real dtMax;
    real dtMin;
    real dtChangeFactor;
    bool adaptiveTimestep;
    size_t timestepCounter;
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
// real jd_to_mjd(const real jd);
void mjd_to_jd(const real mjd, real &jd);
// real mjd_to_jd(const real mjd);
void jd_to_et(const real jd, real &et);
// real jd_to_et(const real jd);
void et_to_jd(const real et, real &jd);
// real et_to_jd(const real et);
void mjd_to_et(const real mjd, real &et);
// real mjd_to_et(const real mjd);
void et_to_mjd(const real et, real &mjd);
// real et_to_mjd(const real et);

void wrap_to_2pi(real &angle);
void rad_to_deg(const real &rad, real &deg);
void deg_to_rad(const real &deg, real &rad);

void vdot(const std::vector<real> &v1, const std::vector<real> &v2, real &dot);
void vnorm(const std::vector<real> &v, real &norm);
void vunit(const std::vector<real> &v, std::vector<real> &unit);
void vcross(const std::vector<real> &v1, const std::vector<real> &v2, std::vector<real> &cross);
void vadd(const std::vector<real> &v1, const std::vector<real> &v2, std::vector<real> &sum);
void vsub(const std::vector<real> &v1, const std::vector<real> &v2, std::vector<real> &diff);
void vcmul(const std::vector<real> &v, const real &c, std::vector<real> &vc);
void vvmul(const std::vector<real> &v1, const std::vector<real> &v2, std::vector<real> &v3);
void vabs_max(const std::vector<real> &v, real &max);

void mat_vec_mul(const std::vector<std::vector<real>> &A, const std::vector<real> &v, std::vector<real> &Av);
void mat_mat_mul(const std::vector<std::vector<real>> &A, const std::vector<std::vector<real>> &B, std::vector<std::vector<real>> &AB);
void rot_mat_x(const real &theta, std::vector<std::vector<real>> &R);
void rot_mat_y(const real &theta, std::vector<std::vector<real>> &R);
void rot_mat_z(const real &theta, std::vector<std::vector<real>> &R);

void kepler_solve(const real &M, const real &e, real &E, const real &tol=1.0e-14L, const int &max_iter=100);
void cometary_to_keplerian(const real &epochMjd, const std::vector<real> &cometaryState, std::vector<real> &keplerianState, const real G);
void keplerian_to_cometary(const real &epochMjd, const std::vector<real> &keplerianState, std::vector<real> &cometaryState, const real G);
void keplerian_to_cartesian(const std::vector<real> &keplerianState, std::vector<real> &cartesianState, const real G);
void cartesian_to_keplerian(const std::vector<real> &cartesianState, std::vector<real> &keplerianState, const real G);
void cometary_to_cartesian(const real &epochMjd, const std::vector<real> &cometaryState, std::vector<real> &cartesianState, const real G);
void cartesian_to_cometary(const real &epochMjd, const std::vector<real> &cartesianState, std::vector<real> &cometaryState, const real G);

#endif
