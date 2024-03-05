#include "utilities.h"

void jd_to_et(const real jd, real &et) {
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    et = (jd - j2000) * day2sec;
}

real jd_to_et(const real jd) {
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    return (jd - j2000) * day2sec;
}

void jd_to_mjd(const real jd, real &mjd) {
    real offset = 2400000.5;
    mjd = jd - offset;
}

real jd_to_mjd(const real jd) {
    real offset = 2400000.5;
    return jd - offset;
}

void et_to_jd(const real et, real &jd) {
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    jd = (et / day2sec) + j2000;
}

real et_to_jd(const real et) {
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    return (et / day2sec) + j2000;
}

void et_to_mjd(const real et, real &mjd) {
    real offset = 2400000.5;
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    mjd = (et / day2sec) - offset + j2000;
}

real et_to_mjd(const real et) {
    real offset = 2400000.5;
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    return (et / day2sec) - offset + j2000;
}

void mjd_to_jd(const real mjd, real &jd) {
    real offset = 2400000.5;
    jd = mjd + offset;
}

real mjd_to_jd(const real mjd) {
    real offset = 2400000.5;
    real jd = mjd + offset;
    return jd;
}

void mjd_to_et(const real mjd, real &et) {
    real offset = 2400000.5;
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    et = (mjd + offset - j2000) * day2sec;
}

real mjd_to_et(const real mjd) {
    real offset = 2400000.5;
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    real et = (mjd + offset - j2000) * day2sec;
    return et;
}

void wrap_to_2pi(real &angle) {
    if (angle < 0) {
        angle += 2 * PI;
    } else if (angle > 2 * PI) {
        angle -= 2 * PI;
    }
}

void rad_to_deg(const real &rad, real &deg) { deg = rad * 180.0 / PI; }

real rad_to_deg(const real rad) { return rad * 180.0 / PI; }

void deg_to_rad(const real &deg, real &rad) { rad = deg * PI / 180.0; }

real deg_to_rad(const real deg) { return deg * PI / 180.0; }

void sort_vector(std::vector<real> &v, const bool &ascending,
                 std::vector<size_t> &sortedIdx) {
    std::iota(sortedIdx.begin(), sortedIdx.end(), 0);
    if (ascending) {
        std::sort(sortedIdx.begin(), sortedIdx.end(),
                  [v](size_t a, size_t b) { return v[a] < v[b]; });
    } else {
        std::sort(sortedIdx.begin(), sortedIdx.end(),
                  [v](size_t a, size_t b) { return v[a] > v[b]; });
    }
    std::vector<real> vCopy = v;
    for (size_t i = 0; i < v.size(); i++) {
        v[i] = vCopy[sortedIdx[i]];
    }
}

void sort_vector_by_idx(std::vector<std::vector<real>> &v,
                            const std::vector<size_t> &sortedIdx) {
    if (v.size() != sortedIdx.size()) {
        throw std::runtime_error(
            "sort_vector_by_idx: v and sortedIdx must be the same size");
    }
    std::vector<std::vector<real>> vCopy = v;
    for (size_t i = 0; i < v.size(); i++) {
        v[i] = vCopy[sortedIdx[i]];
    }
    vCopy.clear();
}

void vdot(const std::vector<real> &v1, const std::vector<real> &v2, real &dot) {
    dot = 0;
    for (size_t i = 0; i < v1.size(); i++) {
        dot += v1[i] * v2[i];
    }
}

void vdot(const real *v1, const real *v2, const size_t &dim, real &dot) {
    dot = 0;
    for (size_t i = 0; i < dim; i++) {
        dot += v1[i] * v2[i];
    }
}

void vnorm(const std::vector<real> &v, real &norm) {
    norm = 0;
    for (size_t i = 0; i < v.size(); i++) {
        norm += v[i] * v[i];
    }
    norm = sqrt(norm);
}

void vnorm(const real *v, const size_t &dim, real &norm) {
    norm = 0;
    for (size_t i = 0; i < dim; i++) {
        norm += v[i] * v[i];
    }
    norm = sqrt(norm);
}

void vunit(const std::vector<real> &v, std::vector<real> &vunit) {
    real norm;
    vnorm(v, norm);
    for (size_t i = 0; i < v.size(); i++) {
        vunit[i] = v[i] / norm;
    }
}

void vunit(const real *v, const size_t &dim, real *unit) {
    real norm;
    vnorm(v, dim, norm);
    for (size_t i = 0; i < dim; i++) {
        unit[i] = v[i] / norm;
    }
}

void vcross(const std::vector<real> &v1, const std::vector<real> &v2,
            std::vector<real> &v3) {
    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void vcross(const real *v1, const real *v2, real *v3) {
    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void vadd(const std::vector<real> &v1, const std::vector<real> &v2,
          std::vector<real> &v3) {
    for (size_t i = 0; i < v1.size(); i++) {
        v3[i] = v1[i] + v2[i];
    }
}

void vsub(const std::vector<real> &v1, const std::vector<real> &v2,
          std::vector<real> &v3) {
    for (size_t i = 0; i < v1.size(); i++) {
        v3[i] = v1[i] - v2[i];
    }
}

void vcmul(const std::vector<real> &v, const real &c, std::vector<real> &vc) {
    for (size_t i = 0; i < v.size(); i++) {
        vc[i] = c * v[i];
    }
}

void vvmul(const std::vector<real> &v1, const std::vector<real> &v2,
           std::vector<real> &v3) {
    for (size_t i = 0; i < v1.size(); i++) {
        v3[i] = v1[i] * v2[i];
    }
}

void vabs_max(const std::vector<real> &v, real &max) {
    max = -1.0e300L;
    for (size_t i = 0; i < v.size(); i++) {
        if (fabs(v[i]) > max) {
            max = fabs(v[i]);
        }
    }
}

void vabs_max(const real *v, const size_t &dim, real &max) {
    max = -1.0e300L;
    for (size_t i = 0; i < dim; i++) {
        if (fabs(v[i]) > max) {
            max = fabs(v[i]);
        }
    }
}

void mat_vec_mul(const std::vector<std::vector<real>> &A,
                 const std::vector<real> &v, std::vector<real> &Av) {
    for (size_t i = 0; i < A.size(); i++) {
        Av[i] = 0;
        for (size_t j = 0; j < A[i].size(); j++) {
            Av[i] += A[i][j] * v[j];
        }
    }
}

void mat_mat_mul(const std::vector<std::vector<real>> &A,
                 const std::vector<std::vector<real>> &B,
                 std::vector<std::vector<real>> &AB) {
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < B[i].size(); j++) {
            AB[i][j] = 0;
            for (size_t k = 0; k < A[i].size(); k++) {
                AB[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void mat3_inv(const std::vector<std::vector<real>> &A,
              std::vector<std::vector<real>> &Ainv) {
    real det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
        A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
        A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    Ainv[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / det;
    Ainv[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) / det;
    Ainv[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / det;
    Ainv[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) / det;
    Ainv[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / det;
    Ainv[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]) / det;
    Ainv[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) / det;
    Ainv[2][1] = (A[0][1] * A[2][0] - A[0][0] * A[2][1]) / det;
    Ainv[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) / det;
}

void mat3_mat3_mul(const real *A, const real *B, real *prod){
    prod[0] = A[0]*B[0] + A[1]*B[3] + A[2]*B[6];
    prod[1] = A[0]*B[1] + A[1]*B[4] + A[2]*B[7];
    prod[2] = A[0]*B[2] + A[1]*B[5] + A[2]*B[8];
    prod[3] = A[3]*B[0] + A[4]*B[3] + A[5]*B[6];
    prod[4] = A[3]*B[1] + A[4]*B[4] + A[5]*B[7];
    prod[5] = A[3]*B[2] + A[4]*B[5] + A[5]*B[8];
    prod[6] = A[6]*B[0] + A[7]*B[3] + A[8]*B[6];
    prod[7] = A[6]*B[1] + A[7]*B[4] + A[8]*B[7];
    prod[8] = A[6]*B[2] + A[7]*B[5] + A[8]*B[8];
}

void mat3_mat3_add(const real *A, const real *B, real *sum){
    for (size_t i = 0; i < 9; i++){
        sum[i] = A[i] + B[i];
    }
}

void rot_mat_x(const real &theta, std::vector<std::vector<real>> &R) {
    R[0][0] = 1;
    R[0][1] = 0;
    R[0][2] = 0;
    R[1][0] = 0;
    R[1][1] = cos(theta);
    R[1][2] = -sin(theta);
    R[2][0] = 0;
    R[2][1] = sin(theta);
    R[2][2] = cos(theta);
}

void rot_mat_y(const real &theta, std::vector<std::vector<real>> &R) {
    R[0][0] = cos(theta);
    R[0][1] = 0;
    R[0][2] = sin(theta);
    R[1][0] = 0;
    R[1][1] = 1;
    R[1][2] = 0;
    R[2][0] = -sin(theta);
    R[2][1] = 0;
    R[2][2] = cos(theta);
}

void rot_mat_z(const real &theta, std::vector<std::vector<real>> &R) {
    R[0][0] = cos(theta);
    R[0][1] = -sin(theta);
    R[0][2] = 0;
    R[1][0] = sin(theta);
    R[1][1] = cos(theta);
    R[1][2] = 0;
    R[2][0] = 0;
    R[2][1] = 0;
    R[2][2] = 1;
}

void LU_decompose(std::vector<std::vector<real>> &A, const size_t &N,
                  const real &tol, size_t *P) {
    for (size_t i = 0; i <= N; i++) {
        P[i] = i;  // Unit permutation matrix, P[N] initialized with N
    }
    size_t j, k, imax;
    real maxA, absA;
    std::vector<real> tmp;
    for (size_t i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;
        for (k = i; k < N; k++) {
            if ((absA = fabs(A[k][i])) > maxA) {
                maxA = absA;
                imax = k;
            }
        }
        // failure, matrix is degenerate
        if (maxA < tol) {
            throw std::runtime_error("LUDecompose: Matrix is degenerate");
        }
        if (imax != i) {
            // pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;
            // pivoting rows of A
            tmp = A[i];
            A[i] = A[imax];
            A[imax] = tmp;
            // counting pivots starting from N (for determinant)
            P[N]++;
        }
        for (j = i + 1; j < N; j++) {
            A[j][i] /= A[i][i];
            for (k = i + 1; k < N; k++) {
                A[j][k] -= A[j][i] * A[i][k];
            }
        }
    }
}

void LU_inverse(std::vector<std::vector<real>> &A, const size_t *P,
                const size_t &N, std::vector<std::vector<real>> &AInv) {
    for (size_t j = 0; j < N; j++) {
        for (size_t i = 0; i < N; i++) {
            AInv[i][j] = (size_t)P[i] == j ? 1 : 0;
            for (size_t k = 0; k < i; k++) {
                AInv[i][j] -= A[i][k] * AInv[k][j];
            }
        }
        for (int i = N - 1; i >= 0; i--) {
            for (size_t k = i + 1; k < N; k++) {
                AInv[i][j] -= A[i][k] * AInv[k][j];
            }
            AInv[i][j] /= A[i][i];
        }
    }
}

void mat_inv(std::vector<std::vector<real>> mat,
             std::vector<std::vector<real>> &matInv,
             const real &tol) {
    const size_t order = mat.size();
    size_t *P = new size_t[order + 1];
    LU_decompose(mat, order, tol, P);
    LU_inverse(mat, P, order, matInv);
    delete[] P;
}
