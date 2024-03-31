#ifndef UTILITIES_H
#define UTILITIES_H

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include "SpiceUsr.h"

#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdexcept>
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

/**
 * @brief Real type to be used for floating point calculations.
 * Choose between double and long double.
 * 
 */
#ifdef LONGDOUBLE
#define real long double
#else
#define real double
#endif
/**
 * @brief Value of PI.
 * 
 */
#define PI 3.141592653589793238462643383279502884197169399375105820974944
/**
* @brief Conversion factor from radians to degrees.
* 
*/
#define RAD2DEG 180.0 / PI
/**
 * @brief Conversion factor from degrees to radians.
 * 
 */
#define DEG2RAD PI / 180.0
/**
 * @brief Value of Earth's obliquity in radians.
 * 
 */
#define EARTH_OBLIQUITY 84381.448 / 3600.0 * DEG2RAD

/**
 * @brief Wrap an angle to the range [0, 2*pi).
 */
void wrap_to_2pi(real &angle);

/**
 * @brief Convert an angle from radians to degrees.
 */
void rad_to_deg(const real &rad, real &deg);

/**
 * @brief Convert an angle from radians to degrees.
 */
real rad_to_deg(const real rad);

/**
 * @brief Convert an angle from degrees to radians.
 */
void deg_to_rad(const real &deg, real &rad);

/**
 * @brief Convert an angle from degrees to radians.
 */
real deg_to_rad(const real deg);

/**
 * @brief Sort a vector in ascending or descending order.
 */
void sort_vector(std::vector<real> &v, const bool &ascending,
                 std::vector<size_t> &sortedIdx);

/**
 * @brief Sort a vector of vectors based on the indices of another vector.
 */
void sort_vector_by_idx(std::vector<std::vector<real>> &v,
                            const std::vector<size_t> &sortedIdx);

/**
 * @brief Dot product of two vectors.
 */
void vdot(const std::vector<real> &v1, const std::vector<real> &v2, real &dot);

/**
 * @brief Dot product of two vectors.
 */
void vdot(const real *v1, const real *v2, const size_t &dim, real &dot);

/**
 * @brief Norm of a vector.
 */
void vnorm(const std::vector<real> &v, real &norm);

/**
 * @brief Norm of a vector.
 */
void vnorm(const real *v, const size_t &dim, real &norm);

/**
 * @brief Unit vector of a vector.
 */
void vunit(const std::vector<real> &v, std::vector<real> &unit);

/**
 * @brief Unit vector of a vector.
 */
void vunit(const real *v, const size_t &dim, real *unit);

/**
 * @brief Cross product of two vectors.
 */
void vcross(const std::vector<real> &v1, const std::vector<real> &v2,
            std::vector<real> &cross);

/**
 * @brief Cross product of two vectors.
 */
void vcross(const real *v1, const real *v2, real *cross);

/**
 * @brief Add two vectors.
 */
void vadd(const std::vector<real> &v1, const std::vector<real> &v2,
          std::vector<real> &sum);

/**
 * @brief Subtract two vectors.
 */
void vsub(const std::vector<real> &v1, const std::vector<real> &v2,
          std::vector<real> &diff);

/**
 * @brief Multiply a vector by a constant.
 */
void vcmul(const std::vector<real> &v, const real &c, std::vector<real> &vc);

/**
 * @brief Multiply two vectors element-wise.
 */
void vvmul(const std::vector<real> &v1, const std::vector<real> &v2,
           std::vector<real> &v3);

/**
 * @brief Find the maximum absolute value in a vector.
 */
void vabs_max(const std::vector<real> &v, real &max);

/**
 * @brief Find the maximum absolute value in a vector.
 */
void vabs_max(const real *v, const size_t &dim, real &max);

/**
 * @brief Multiply a matrix by a vector.
 */
void mat_vec_mul(const std::vector<std::vector<real>> &A,
                 const std::vector<real> &v, std::vector<real> &Av);

/**
 * @brief Multiply two matrices.
 */
void mat_mat_mul(const std::vector<std::vector<real>> &A,
                 const std::vector<std::vector<real>> &B,
                 std::vector<std::vector<real>> &AB);

/**
 * @brief Invert a 3x3 matrix.
 */
void mat3_inv(const std::vector<std::vector<real>> &A,
              std::vector<std::vector<real>> &Ainv);

/**
 * @brief Multiply two 3x3 matrices.
 */
void mat3_mat3_mul(const real *A, const real *B, real *prod);

/**
 * @brief Add two 3x3 matrices.
 */
void mat3_mat3_add(const real *A, const real *B, real *sum);

/**
 * @brief Rotation matrix about the x-axis.
 */
void rot_mat_x(const real &theta, std::vector<std::vector<real>> &R);

/**
 * @brief Rotation matrix about the y-axis.
 */
void rot_mat_y(const real &theta, std::vector<std::vector<real>> &R);

/**
 * @brief Rotation matrix about the z-axis.
 */
void rot_mat_z(const real &theta, std::vector<std::vector<real>> &R);

/**
 * @brief Decompose a matrix into its LU form.
 */
void LU_decompose(std::vector<std::vector<real>> &A, const size_t &N,
                  const real &tol, size_t *P);

/**
 * @brief Invert a matrix using its LU decomposition.
 */
void LU_inverse(std::vector<std::vector<real>> &A, const size_t *P,
                const size_t &N, std::vector<std::vector<real>> &AInv);

/**
 * @brief Invert a matrix.
 */
void mat_inv(std::vector<std::vector<real>> mat,
             std::vector<std::vector<real>> &matInv,
             const real &tol = 1.0e-16L);

#endif
