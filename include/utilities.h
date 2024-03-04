#ifndef UTILITIES_H
#define UTILITIES_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
using std::acos;
using std::asin;
using std::asinh;
using std::atan2;
using std::cos;
using std::cosh;
using std::fabs;
using std::fmax;
using std::fmin;
using std::log;
using std::pow;
using std::sin;
using std::sinh;
using std::sqrt;
using std::tan;
using std::tanh;

#include "SpiceUsr.h"
#include "spk.h"

#ifdef LONGDOUBLE
#define real long double
#else
#define real double
#endif
#define PI 3.141592653589793238462643383279502884197169399375105820974944
#define RAD2DEG 180.0 / PI
#define DEG2RAD PI / 180.0
#define EARTH_OBLIQUITY 84381.448 / 3600.0 * DEG2RAD

// forward declarations for simulation.h
struct Body;
class SpiceBody;
class IntegBody;
class Event;
class ImpulseEvent;
class propSimulation;

struct Constants {
    real du2m = 149597870700.0L;  // default au to m
    real tu2s = 86400.0;          // default day to sec
    real duptu2mps = du2m/tu2s;   // default au/day to m/s
    real G = 6.6743e-11L / (du2m * du2m * du2m) * tu2s *
        tu2s;                                  // default kg au^3 / day^2
    real clight = 299792458.0L / du2m * tu2s;  // default au/day
    real j2000Jd = 2451545.0;
    real JdMinusMjd = 2400000.5;
};

struct IntegrationParameters {
    size_t nInteg;
    size_t nSpice;
    size_t nTotal;
    size_t n2Derivs;
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

struct NongravParamaters {
    // from https://ssd.jpl.nasa.gov/horizons/manual.html
    real a1 = 0.0L;
    real a2 = 0.0L;
    real a3 = 0.0L;
    real alpha = 0.1112620426L;
    real k = 4.6142L;
    real m = 2.15L;
    real n = 5.093L;
    real r0_au = 2.808L;
};

struct InterpolationParameters {
    std::vector<real> tStack;
    std::vector<std::vector<real>> xIntegStack;
    std::vector<std::vector<std::vector<real>>> bStack;
    std::vector<std::vector<real>> accIntegStack;
};

struct ForceParameters {
    std::vector<real> masses;
    std::vector<real> radii;
    std::vector<int> spiceIdList;
    std::vector<NongravParamaters> ngParamsList;
    std::vector<bool> isPPNList;
    std::vector<bool> isJ2List;
    std::vector<real> J2List;
    std::vector<real> poleRAList;
    std::vector<real> poleDecList;
    std::vector<bool> isNongravList;
    std::vector<bool> isMajorList;
    std::vector<bool> isThrustingList;
};

struct BPlaneParameters {
    real x;
    real y;
    real z;
    std::vector<real> dx = std::vector<real>(6, 0.0L);
    std::vector<real> dy = std::vector<real>(6, 0.0L);
};

class CloseApproachParameters {
   private:
   public:
    real t;
    std::vector<real> xRel = std::vector<real>(6, 0.0L);
    real tMap;
    std::vector<real> xRelMap = std::vector<real>(6, 0.0L);
    real dist;
    real vel;
    real vInf;
    std::string flybyBody;
    int flybyBodyIdx;
    std::string centralBody;
    int centralBodyIdx;
    int centralBodySpiceId;
    bool impact;
    real tPeri;
    real tLin;
    std::vector<real> bVec = std::vector<real>(3, 0.0L);
    real bMag;
    real gravFocusFactor;
    BPlaneParameters kizner;
    BPlaneParameters opik;
    BPlaneParameters scaled;
    BPlaneParameters mtp;
    std::vector<real> dTLinMinusT = std::vector<real>(6, 0.0L);
    std::vector<real> dt = std::vector<real>(6, 0.0L);
    void get_ca_parameters(propSimulation *propSim, const real &tMap);
    void print_summary(int prec=8);
};

class ImpactParameters : public CloseApproachParameters {
   private:
   public:
    std::vector<real> xRelBodyFixed = std::vector<real>(6, 0.0L);
    real lon;
    real lat;
    real alt;
    void get_impact_parameters(propSimulation *propSim);
    void print_summary(int prec=8);
};

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

void wrap_to_2pi(real &angle);
void rad_to_deg(const real &rad, real &deg);
real rad_to_deg(const real rad);
void deg_to_rad(const real &deg, real &rad);
real deg_to_rad(const real deg);

void sort_vector(std::vector<real> &v, const bool &ascending,
                 std::vector<size_t> &sortedIdx);
void sort_vector_by_idx(std::vector<std::vector<real>> &v,
                            const std::vector<size_t> &sortedIdx);
void vdot(const std::vector<real> &v1, const std::vector<real> &v2, real &dot);
void vdot(const real *v1, const real *v2, const size_t &dim, real &dot);
void vnorm(const std::vector<real> &v, real &norm);
void vnorm(const real *v, const size_t &dim, real &norm);
void vunit(const std::vector<real> &v, std::vector<real> &unit);
void vunit(const real *v, const size_t &dim, real *unit);
void vcross(const std::vector<real> &v1, const std::vector<real> &v2,
            std::vector<real> &cross);
void vcross(const real *v1, const real *v2, real *cross);
void vadd(const std::vector<real> &v1, const std::vector<real> &v2,
          std::vector<real> &sum);
void vsub(const std::vector<real> &v1, const std::vector<real> &v2,
          std::vector<real> &diff);
void vcmul(const std::vector<real> &v, const real &c, std::vector<real> &vc);
void vvmul(const std::vector<real> &v1, const std::vector<real> &v2,
           std::vector<real> &v3);
void vabs_max(const std::vector<real> &v, real &max);
void vabs_max(const real *v, const size_t &dim, real &max);

void mat_vec_mul(const std::vector<std::vector<real>> &A,
                 const std::vector<real> &v, std::vector<real> &Av);
void mat_mat_mul(const std::vector<std::vector<real>> &A,
                 const std::vector<std::vector<real>> &B,
                 std::vector<std::vector<real>> &AB);
void mat3_inv(const std::vector<std::vector<real>> &A,
              std::vector<std::vector<real>> &Ainv);
void mat3_mat3_mul(const real *A, const real *B, real *prod);
void mat3_mat3_add(const real *A, const real *B, real *sum);
void rot_mat_x(const real &theta, std::vector<std::vector<real>> &R);
void rot_mat_y(const real &theta, std::vector<std::vector<real>> &R);
void rot_mat_z(const real &theta, std::vector<std::vector<real>> &R);
void LU_decompose(std::vector<std::vector<real>> &A, const size_t &N,
                  const real &tol, size_t *P);
void LU_inverse(std::vector<std::vector<real>> &A, const size_t *P,
                const size_t &N, std::vector<std::vector<real>> &AInv);
void mat_inv(std::vector<std::vector<real>> mat,
             std::vector<std::vector<real>> &matInv,
             const real &tol = 1.0e-16L);
#endif
