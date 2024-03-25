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
