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

void sort_vector(std::vector<real> &v, const bool &ascending) {
    if (ascending) {
        std::stable_sort(v.begin(), v.end());
    } else {
        std::stable_sort(v.begin(), v.end(), std::greater<real>());
    }
}

void sort_vector_by_another(std::vector<real> &v, const std::vector<real> &vRef,
                            const bool &ascending) {
    if (v.size() != vRef.size()) {
        throw std::runtime_error(
            "sort_vector_by_another: v and vRef must be the same size");
    }
    std::vector<size_t> sortedIdx(v.size());
    std::iota(sortedIdx.begin(), sortedIdx.end(), 0);
    if (ascending) {
        std::sort(sortedIdx.begin(), sortedIdx.end(),
                  [&vRef](size_t a, size_t b) { return vRef[a] < vRef[b]; });
    } else {
        std::sort(sortedIdx.begin(), sortedIdx.end(),
                  [&vRef](size_t a, size_t b) { return vRef[a] > vRef[b]; });
    }
    std::vector<real> vCopy = v;
    for (size_t i = 0; i < v.size(); i++) {
        v[i] = vCopy[sortedIdx[i]];
    }
    vCopy.clear();
}

void sort_vector_by_another(std::vector<std::vector<real>> &v,
                            const std::vector<real> &vRef,
                            const bool &ascending) {
    if (v.size() != vRef.size()) {
        throw std::runtime_error(
            "sort_vector_by_another: v and vRef must be the same size");
    }
    std::vector<size_t> sortedIdx(v.size());
    std::iota(sortedIdx.begin(), sortedIdx.end(), 0);
    if (ascending) {
        std::sort(sortedIdx.begin(), sortedIdx.end(),
                  [&vRef](size_t a, size_t b) { return vRef[a] < vRef[b]; });
    } else {
        std::sort(sortedIdx.begin(), sortedIdx.end(),
                  [&vRef](size_t a, size_t b) { return vRef[a] > vRef[b]; });
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

void kepler_solve(const real &M, const real &e, real &E, const real &tol,
                  const int &max_iter) {
    if (e < 0.8) {
        E = M;
    } else {
        E = PI;
    }
    int iter = 0;
    real F = E - e * sin(E) - M;
    real F_prime = 1 - e * cos(E);
    while ((fabs(F) > tol && iter < max_iter)) {
        E -= (F / F_prime);
        F = E - e * sin(E) - M;
        F_prime = 1 - e * cos(E);
        iter++;
    }
    if (iter == max_iter) {
        std::cout << "utilities.cpp: WARNING: kepler_solve did not converge in "
                  << max_iter << " iterations!!!"
                  << " F: " << F << std::endl;
    }
    // std::cout << "iter: " << iter << std::endl;
    // std::cout << "E: " << E*RAD2DEG << std::endl;
    // std::cout << "M: " << M*RAD2DEG << std::endl;
}

void kepler_solve_hyperbolic(const real &M, const real &e, real &EHyp,
                             const real &tol, const int &max_iter) {
    // EHyp = log(2*M/e+1.8);
    EHyp = M;
    int iter = 0;
    real F = e * sinh(EHyp) - EHyp - M;
    real F_prime = e * cosh(EHyp) - 1;
    while ((fabs(F) > tol && iter < max_iter)) {
        EHyp -= (F / F_prime);
        F = e * sinh(EHyp) - EHyp - M;
        F_prime = e * cosh(EHyp) - 1;
        iter++;
    }
    if (iter == max_iter) {
        std::cout << "utilities.cpp: WARNING: kepler_solve_hyperbolic did not "
                     "converge in "
                  << max_iter << " iterations!!!"
                  << " F: " << F << std::endl;
    }
    // std::cout << "iter: " << iter << std::endl;
    // std::cout << "EHyp: " << EHyp*RAD2DEG << std::endl;
    // std::cout << "M: " << M*RAD2DEG << std::endl;
}

void cometary_to_keplerian(const real &epochMjD,
                           const std::vector<real> &cometaryState,
                           std::vector<real> &keplerianState, const real GM) {
    real a = cometaryState[1] / (1 - cometaryState[0]);
    real e = cometaryState[0];
    real n, M, E, nu;
    if (e < 1) {
        n = sqrt(GM / pow(a, 3.0L));
        M = n * (epochMjD - cometaryState[2]);
        wrap_to_2pi(M);
        kepler_solve(M, cometaryState[0], E);
        nu = 2 * atan2(tan(E / 2) * sqrt(1 + e), sqrt(1 - e));
        wrap_to_2pi(nu);
    } else if (e > 1) {
        n = sqrt(-GM / pow(a, 3.0L));
        M = n * (epochMjD - cometaryState[2]);
        wrap_to_2pi(M);
        kepler_solve_hyperbolic(M, cometaryState[0], E);
        nu = 2 * atan2(tanh(E / 2) * sqrt(e + 1), sqrt(e - 1));
    } else {
        throw std::runtime_error(
            "utilities.cpp: cometary_to_keplerian: Cannot handle e = 1 right "
            "now!!!");
    }
    // std::cout << "a: " << a << std::endl;
    // std::cout << "M: " << M*RAD2DEG << std::endl;
    // std::cout << "E: " << E*RAD2DEG << std::endl;
    // std::cout << "nu: " << nu*RAD2DEG << std::endl;
    keplerianState[0] = a;
    keplerianState[1] = cometaryState[0];
    keplerianState[2] = cometaryState[5];
    keplerianState[3] = cometaryState[3];
    keplerianState[4] = cometaryState[4];
    keplerianState[5] = nu;
}

void keplerian_to_cometary(const real &epochMjD,
                           const std::vector<real> &keplerianState,
                           std::vector<real> &cometaryState, const real GM) {
    real a = keplerianState[0];
    real e = keplerianState[1];
    real nu = keplerianState[5];
    real E = 2 * atan2(tan(nu / 2) * sqrt(1 - e), sqrt(1 + e));
    real M = E - e * sin(E);
    real n = sqrt(GM / pow(a, 3.0L));
    real T0 = epochMjD - (M / n);

    cometaryState[0] = e;
    cometaryState[1] = a * (1 - e);
    cometaryState[2] = T0;
    cometaryState[3] = keplerianState[3];
    cometaryState[4] = keplerianState[4];
    cometaryState[5] = keplerianState[2];
}

void keplerian_to_cartesian(const std::vector<real> &keplerianState,
                            std::vector<real> &cartesianState, const real GM) {
    real a = keplerianState[0];
    real e = keplerianState[1];
    real i = keplerianState[2];
    real Omega = keplerianState[3];
    real omega = keplerianState[4];
    real nu = keplerianState[5];
    real p = a * (1 - pow(e, 2.0L));
    real r = p / (1 + e * cos(nu));
    // real E = 2*atan(tan(nu/2)*sqrt((1-e)/(1+e)));
    // real M = E-e*sin(E);
    // real n = sqrt(GM/pow(a, 3.0L));
    // real T0 = epochMjd-(M/n);
    // ConstSpiceDouble elts[8] = {(double) (a*(1-e)), (double)e, (double)i,
    // (double)Omega, (double)omega, (double)M, 0.0, (double)GM}; SpiceDouble
    // state[6]; conics_c(elts, 0.0, state); std::cout << "spice state:" <<
    // std::endl; std::cout << state[0] << " " << state[1] << " " << state[2] <<
    // " " << state[3] << " " << state[4] << " " << state[5] << std::endl;

    std::vector<std::vector<real>> R1(3, std::vector<real>(3));
    std::vector<std::vector<real>> R2(3, std::vector<real>(3));
    std::vector<std::vector<real>> R3(3, std::vector<real>(3));
    std::vector<std::vector<real>> RTemp(3, std::vector<real>(3));
    std::vector<std::vector<real>> R(3, std::vector<real>(3));
    std::vector<real> r_temp(3);
    std::vector<real> v_temp(3);
    std::vector<real> r_final(3);
    std::vector<real> v_final(3);

    rot_mat_z(Omega, R1);
    rot_mat_x(i, R2);
    rot_mat_z(omega, R3);
    mat_mat_mul(R1, R2, RTemp);
    mat_mat_mul(RTemp, R3, R);

    r_temp[0] = r * cos(nu);
    r_temp[1] = r * sin(nu);
    r_temp[2] = 0;

    v_temp[0] = -sqrt(GM / p) * sin(nu);
    v_temp[1] = sqrt(GM / p) * (e + cos(nu));
    v_temp[2] = 0;

    mat_vec_mul(R, r_temp, r_final);
    mat_vec_mul(R, v_temp, v_final);
    cartesianState[0] = r_final[0];
    cartesianState[1] = r_final[1];
    cartesianState[2] = r_final[2];
    cartesianState[3] = v_final[0];
    cartesianState[4] = v_final[1];
    cartesianState[5] = v_final[2];
}

void cartesian_to_keplerian(const std::vector<real> &cartesianState,
                            std::vector<real> &keplerianState, const real GM) {
    std::vector<real> rVec(3);
    std::vector<real> vVec(3);
    rVec[0] = cartesianState[0];
    rVec[1] = cartesianState[1];
    rVec[2] = cartesianState[2];
    vVec[0] = cartesianState[3];
    vVec[1] = cartesianState[4];
    vVec[2] = cartesianState[5];
    real r;
    real v;
    vnorm(rVec, r);
    vnorm(vVec, v);

    std::vector<real> hVec(3);
    vcross(rVec, vVec, hVec);
    std::vector<real> nVec(3);
    vcross({0, 0, 1}, hVec, nVec);
    std::vector<real> eVecTemp1(3);
    std::vector<real> eVecTemp2(3);
    std::vector<real> rHatVec(3);
    std::vector<real> eVec(3);
    vcross(vVec, hVec, eVecTemp1);
    vcmul(eVecTemp1, 1.0L / GM, eVecTemp2);
    vunit(rVec, rHatVec);
    vsub(eVecTemp2, rHatVec, eVec);

    real h;
    vnorm(hVec, h);
    real n;
    vnorm(nVec, n);
    real e;
    vnorm(eVec, e);
    real a = h * h / (GM * (1 - e * e));

    real i = acos(hVec[2] / h);
    if (i > M_PI / 2) {
        i = M_PI - i;
    }
    real Omega = acos(nVec[0] / n);
    if (nVec[1] < 0) {
        Omega = 2 * M_PI - Omega;
    }
    real omega = acos(
        (nVec[0] * eVec[0] + nVec[1] * eVec[1] + nVec[2] * eVec[2]) / (n * e));
    if (eVec[2] < 0) {
        omega = 2 * M_PI - omega;
    }
    real nu = acos((eVec[0] * cartesianState[0] + eVec[1] * cartesianState[1] +
                    eVec[2] * cartesianState[2]) /
                   (e * r));
    if (cartesianState[2] < 0) {
        nu = 2 * M_PI - nu;
    }
    // real E = 2*atan(tan(nu/2)*sqrt((1-e)/(1+e)));
    // real M = E-e*sin(E);

    keplerianState[0] = a;
    keplerianState[1] = e;
    keplerianState[2] = i;
    keplerianState[3] = Omega;
    keplerianState[4] = omega;
    keplerianState[5] = nu;
}

void cometary_to_cartesian(const real &epochMjd,
                           const std::vector<real> &cometaryState,
                           std::vector<real> &cartesianState, const real GM) {
    std::vector<real> keplerianState(6);
    cometary_to_keplerian(epochMjd, cometaryState, keplerianState, GM);
    keplerian_to_cartesian(keplerianState, cartesianState, GM);
}

void cartesian_to_cometary(const real &epochMjd,
                           const std::vector<real> &cartesianState,
                           std::vector<real> &cometaryState, const real GM) {
    std::vector<real> keplerianState(6);
    cartesian_to_keplerian(cartesianState, keplerianState, GM);
    keplerian_to_cometary(epochMjd, keplerianState, cometaryState, GM);
}

void cartesian_cometary_partials(const real &epochMjd,
                                 const std::vector<real> &cometaryState,
                                 std::vector<std::vector<real>> &partials,
                                 const real GM) {
    std::vector<real> keplerianState(6);
    partials = std::vector<std::vector<real>>(6, std::vector<real>(6, 0.0L));
    cometary_to_keplerian(epochMjd, cometaryState, keplerianState, GM);
    const real tp = cometaryState[2];
    const real a = keplerianState[0];
    const real e = keplerianState[1];
    const real inc = keplerianState[2];
    const real Omega = keplerianState[3];
    const real omega = keplerianState[4];
    const real nu = keplerianState[5];
    const real E = 2 * atan2(tan(nu / 2) * sqrt(1 - e), sqrt(1 + e));
    real *part = new real[6];
    // dCartde(GM, a, e, inc, Omega, omega, E, part);
    dCartdeNum(epochMjd, GM, cometaryState, part);
    for (size_t i = 0; i < 6; i++) {
        partials[i][0] = part[i];
    }
    dCartdq(epochMjd, tp, GM, a, e, inc, Omega, omega, E, part);
    for (size_t i = 0; i < 6; i++) {
        partials[i][1] = part[i];
    }
    dCartdTp(epochMjd, tp, GM, a, e, inc, Omega, omega, E, part);
    for (size_t i = 0; i < 6; i++) {
        partials[i][2] = part[i];
    }
    dCartdOmega(GM, a, e, inc, Omega, omega, E, part);
    for (size_t i = 0; i < 6; i++) {
        partials[i][3] = part[i];
    }
    dCartdomega(GM, a, e, inc, Omega, omega, E, part);
    for (size_t i = 0; i < 6; i++) {
        partials[i][4] = part[i];
    }
    dCartdinc(GM, a, e, inc, Omega, omega, E, part);
    for (size_t i = 0; i < 6; i++) {
        partials[i][5] = part[i];
    }
    delete[] part;
}

// partials from 2 sources:
// 1: Murison, https://www.researchgate.net/publication/271214791_Partial_Derivatives_of_Observables_with_Respect_to_Two-Body_Orbital_Elements
// 2: Farnocchia et al, 10.1007/s10569-013-9476-9
void dCartdeNum(const real &t, const real &GM,
                const std::vector<real> &cometaryState, real *partial){
    const real e = cometaryState[0];
    const real delta = 1e-8;
    const real pert = delta*e;
    std::vector<real> comStatePlus = cometaryState;
    std::vector<real> comStateMinus = cometaryState;
    comStatePlus[0] += pert;
    comStateMinus[0] -= pert;
    std::vector<real> cartStatePlus(6);
    std::vector<real> cartStateMinus(6);
    cometary_to_cartesian(t, comStatePlus, cartStatePlus, GM);
    cometary_to_cartesian(t, comStateMinus, cartStateMinus, GM);
    for (size_t i = 0; i < 6; i++){
        partial[i] = (cartStatePlus[i] - cartStateMinus[i])/(2.0L*pert);
    }
}

void dCartde(const real &GM, const real &a, const real &e, const real &inc,
             const real &Omega, const real &omega, const real &E,
             real *partial) {
    const real n = sqrt(GM / a / a / a);
    const real cO = cos(Omega);
    const real sO = sin(Omega);
    const real co = cos(omega);
    const real so = sin(omega);
    const real ci = cos(inc);
    const real si = sin(inc);
    const real cE = cos(E);
    const real sE = sin(E);
    const real posFac = e / sqrt(1 - e * e);
    partial[0] = a*(-cO*co + sO*so*ci + posFac*(cO*so+sO*co*ci)*sE);
    partial[1] = a*(-sO*co - cO*so*ci + posFac*(sO*so-cO*co*ci)*sE);
    partial[2] = -a*(so + posFac*co*sE)*si;
    const real velFac1 = n*a*cE/(1-e*cos(E))/(1-e*cos(E));
    const real velFac2 = (cE-e)/sqrt(1-e*e);
    partial[3] = velFac1*((-cO*co + sO*so*ci)*sE - velFac2*(cO*so+sO*co*ci));
    partial[4] = velFac1*(-(sO*co + cO*so*ci)*sE - velFac2*(sO*so-cO*co*ci));
    partial[5] = velFac1*(-so*sE - velFac2*co)*si;
}

void dCartde2(const real &GM, const real &a, const real &e, const real &inc,
             const real &Omega, const real &omega, const real &E,
             real *partial) {
    const real t2 = cos(Omega);
    const real t3 = cos(omega);
    const real t5 = sin(Omega);
    const real t6 = cos(inc);
    const real t7 = t5 * t6;
    const real t8 = sin(omega);
    const real t10 = -t2 * t3 + t7 * t8;
    const real t11 = e * e;
    const real t13 = sqrt(0.1e1 - t11);
    const real t15 = sin(E);
    const real t16 = t15 * e;
    const real t19 = t2 * t8 + t3 * t7;
    const real t22 = 0.1e1 / t13;
    const real t25 = t2 * t6;
    const real t27 = -t5 * t3 - t25 * t8;
    const real t31 = t5 * t8 - t25 * t3;
    const real t35 = sin(inc);
    const real t37 = t8 * t13;
    const real t43 = cos(E);
    const real t45 = sqrt(a);
    const real t47 = t22 * t43 / t45;
    const real t50 = -t43 + e;
    const real t53 = sqrt(GM);
    const real t57 = pow(t43 * e - 0.1e1, 0.2e1);
    const real t58 = 0.1e1 / t57;
    partial[0] = (t10 * t13 + t16 * t19) * t22 * a;
    partial[1] = (t27 * t13 + t16 * t31) * t22 * a;
    partial[2] = -a * t35 * (t37 + t3 * t15 * e) * t22;
    partial[3] = t47 * (t15 * t10 * t13 + t50 * t19) * t53 * t58;
    partial[4] = t47 * t53 * (t15 * t27 * t13 + t31 * t50) * t58;
    partial[5] = t47 * (-t37 * t15 - t3 * t50) * t53 * t35 * t58;
}

void dCartde3(const real &a, const real &e, const real &inc,
             const real &Omega, const real &omega, const real &E,
             real *partial) {
    std::vector<std::vector<real>> R1(3, std::vector<real>(3));
    std::vector<std::vector<real>> R2(3, std::vector<real>(3));
    std::vector<std::vector<real>> R3(3, std::vector<real>(3));
    std::vector<std::vector<real>> RTemp(3, std::vector<real>(3));
    std::vector<std::vector<real>> Q(3, std::vector<real>(3));

    rot_mat_z(Omega, R1);
    rot_mat_x(inc, R2);
    rot_mat_z(omega, R3);
    mat_mat_mul(R1, R2, RTemp);
    mat_mat_mul(RTemp, R3, Q);

    std::vector<real> drdeFrame(3);
    drdeFrame[0] = (a*cos(2*atan((tan(E/2)*sqrt(e+1))/sqrt(1-e)))*(e*e - 1)*(cos(2*atan((tan(E/2)*sqrt(e+1))/sqrt(1-e))) + (2*e*sin(2*atan((tan(E/2)*sqrt(e+1))/sqrt(1-e)))*(tan(E/2)/(2*sqrt(1-e)*sqrt(e+1)) + (tan(E/2)*sqrt(e+1))/(2*pow(1-e, 1.5L))))/((tan(E/2)*tan(E/2)*(e + 1))/(e - 1) - 1)))/pow(e*cos(2*atan((tan(E/2)*sqrt(e+1))/sqrt(1-e))) + 1, 2.0L) - (2*a*e*cos(2*atan((tan(E/2)*sqrt(e+1))/sqrt(1-e))))/(e*cos(2*atan((tan(E/2)*sqrt(e+1))/sqrt(1-e))) + 1) - (2*a*sin(2*atan((tan(E/2)*sqrt(e+1))/sqrt(1-e)))*(e*e - 1)*(tan(E/2)/(2*sqrt(1-e)*sqrt(e+1)) + (tan(E/2)*sqrt(e+1))/(2*pow(1-e, 1.5L))))/(((tan(E/2)*tan(E/2)*(e + 1))/(e - 1) - 1)*(e*cos(2*atan((tan(E/2)*sqrt(e+1))/sqrt(1-e))) + 1));
    drdeFrame[1] = (a*sin(2*atan((tan(E/2)*sqrt(e+1))/sqrt(1-e)))*(e*e - 1)*(cos(2*atan((tan(E/2)*sqrt(e+1))/sqrt(1-e))) + (2*e*sin(2*atan((tan(E/2)*sqrt(e+1))/sqrt(1-e)))*(tan(E/2)/(2*sqrt(1-e)*sqrt(e+1)) + (tan(E/2)*sqrt(e+1))/(2*pow(1-e, 1.5L))))/((tan(E/2)*tan(E/2)*(e + 1))/(e - 1) - 1)))/pow(e*cos(2*atan((tan(E/2)*sqrt(e+1))/sqrt(1-e))) + 1, 2.0L) - (2*a*e*sin(2*atan((tan(E/2)*sqrt(e+1))/sqrt(1-e))))/(e*cos(2*atan((tan(E/2)*sqrt(e+1))/sqrt(1-e))) + 1) + (2*a*cos(2*atan((tan(E/2)*sqrt(e+1))/sqrt(1-e)))*(e*e - 1)*(tan(E/2)/(2*sqrt(1-e)*sqrt(e+1)) + (tan(E/2)*sqrt(e+1))/(2*pow(1-e, 1.5L))))/(((tan(E/2)*tan(E/2)*(e + 1))/(e - 1) - 1)*(e*cos(2*atan((tan(E/2)*sqrt(e+1))/sqrt(1-e))) + 1));
    drdeFrame[2] = 0;

    std::vector<real> drde(3);
    mat_vec_mul(Q, drdeFrame, drde);
    partial[0] = drde[0];
    partial[1] = drde[1];
    partial[2] = drde[2];
    partial[3] = 0;
    partial[4] = 0;
    partial[5] = 0;
}

void dCartdq(const real &t, const real &tp, const real &GM, const real &a,
             const real &e, const real &inc, const real &Omega,
             const real &omega, const real &E, real *partial) {
    const real t2 = cos(Omega);
    const real t3 = sin(omega);
    const real t5 = sin(Omega);
    const real t6 = cos(inc);
    const real t7 = t6 * t5;
    const real t8 = cos(omega);
    const real t10 = t2 * t3 + t7 * t8;
    const real t11 = sin(E);
    const real t12 = t11 * e;
    const real t13 = t - tp;
    const real t17 = cos(E);
    const real n = sqrt(GM / a / a / a);
    const real t19 = (t12 + 0.3e1 / 0.2e1 * n * t13) * t17 - t11;
    const real t21 = e * e;
    const real t23 = sqrt(0.1e1 - t21);
    const real t27 = -t2 * t8 + t7 * t3;
    const real t28 = t17 * t17;
    const real t29 = t28 * e;
    const real t30 = t21 + 0.1e1;
    const real t32 = n * t11;
    const real t34 = 0.3e1 / 0.2e1 * t32 * t13;
    const real t35 = -t29 + t30 * t17 + t34 - e;
    const real t40 = t17 * e - 0.1e1;
    const real t42 = 0.1e1 / t40 / 0.2e1;
    const real t44 = t2 * t6;
    const real t46 = t5 * t3 - t44 * t8;
    const real t51 = t5 * t8 + t44 * t3;
    const real t62 = sin(inc);
    const real t68 = -t13;
    const real t71 = t28 * t17 * t21 - 0.2e1 * t29 + t17 + 0.3e1 * t32 * t68;
    const real t84 = t11 * t21 * t28 + (-0.2e1 * t12 - 0.3e1 * n * t68) * t17 + t11 +
            0.3e1 * n * t68 * e;
    const real t87 = sqrt(GM);
    const real t89 = sqrt(a);
    const real t91 = 0.1e1 / t89 / a;
    const real t92 = t40 * t40;
    const real t94 = 0.1e1 / t92 / t40;
    const real dadq = 1/(1-e);
    partial[0] = (-0.2e1 * t10 * t19 * t23 + 0.2e1 * t27 * t35) * t42 * dadq;
    partial[1] = (-0.2e1 * t46 * t19 * t23 - 0.2e1 * t35 * t51) * t42 * dadq;
    partial[2] = 0.2e1 * (t8 * t19 * t23 + (t29 - t30 * t17 - t34 + e) * t3) * t62 * t42 * dadq;
    partial[3] = (-t71 * t10 * t23 + t27 * t84) * t87 * t91 * t94 / 0.2e1 * dadq;
    partial[4] = -t91 * t87 * (t71 * t46 * t23 + t84 * t51) * t94 / 0.2e1 * dadq;
    partial[5] = -(-t71 * t8 * t23 + t3 * t84) * t91 * t87 * t62 * t94 / 0.2e1 * dadq;
}

void dCartdTp(const real &t, const real &tp, const real &GM, const real &a,
              const real &e, const real &inc, const real &Omega,
              const real &omega, const real &E, real *partial) {
    const real t2 = cos(E);
    const real t3 = cos(Omega);
    const real t4 = sin(omega);
    const real t6 = sin(Omega);
    const real t7 = cos(inc);
    const real t8 = t7 * t6;
    const real t9 = cos(omega);
    const real t11 = t3 * t4 + t8 * t9;
    const real t13 = e * e;
    const real t15 = sqrt(0.1e1 - t13);
    const real t19 = t3 * t9 - t8 * t4;
    const real t20 = sin(E);
    const real t24 = t3 * t7;
    const real t26 = t4 * t6 - t24 * t9;
    const real t31 = t6 * t9 + t24 * t4;
    const real t34 = sin(inc);
    const real t37 = t15 * t9;
    const real t40 = sqrt(a);
    const real t41 = 0.1e1 / t40;
    const real t45 = t2 - e;
    const real t49 = sqrt(GM);
    const real t52 = pow(t2 * e - 0.1e1, 0.2e1);
    const real t53 = 0.1e1 / t52;
    const real t54 = t49 * t53;
    const real n = sqrt(GM / a / a / a);
    const real r = a * (1 - e * t2);
    const real dEda = -3*n*(t-tp)/(2*r);
    const real dadq = 1/(1-e);
    const real dqdTp = 2*a*(1-e)/(3*(t-tp));
    const real fac = dEda*dadq*dqdTp;
    partial[0] = (-t2 * t11 * t15 - t19 * t20) * a * fac;
    partial[1] = (-t2 * t26 * t15 - t31 * t20) * a * fac;
    partial[2] = a * t34 * (-t20 * t4 + t37 * t2) * fac;
    partial[3] = t41 * (t20 * t11 * t15 - t19 * t45) * t54 * fac;
    partial[4] = -t41 * (t31 * t45 - t26 * t15 * t20) * t54 * fac;
    partial[5] = -t41 * t49 * t34 * (t37 * t20 + t4 * t45) * t53 * fac;
}

void dCartdOmega(const real &GM, const real &a, const real &e, const real &inc,
                 const real &Omega, const real &omega, const real &E,
                 real *partial) {
    const real t2 = sin(Omega);
    const real t3 = sin(omega);
    const real t5 = cos(Omega);
    const real t6 = cos(inc);
    const real t7 = t5 * t6;
    const real t8 = cos(omega);
    const real t11 = e * e;
    const real t13 = sqrt(0.1e1 - t11);
    const real t14 = (-t3 * t2 + t8 * t7) * t13;
    const real t15 = sin(E);
    const real t17 = cos(E);
    const real t18 = -t17 + e;
    const real t21 = t2 * t8 + t7 * t3;
    const real t25 = t2 * t6;
    const real t27 = t3 * t5 + t25 * t8;
    const real t32 = t5 * t8 - t25 * t3;
    const real t39 = sqrt(GM);
    const real t41 = sqrt(a);
    const real t42 = 0.1e1 / t41;
    const real t45 = 0.1e1 / (t17 * e - 0.1e1);
    partial[0] = (-t15 * t14 + t18 * t21) * a;
    partial[1] = (-t15 * t27 * t13 - t32 * t18) * a;
    partial[2] = 0.0e0;
    partial[3] = (-t21 * t15 + t14 * t17) * t39 * t42 * t45;
    partial[4] = t42 * (t17 * t27 * t13 + t32 * t15) * t39 * t45;
    partial[5] = 0.0e0;
}

void dCartdomega(const real &GM, const real &a, const real &e, const real &inc,
                 const real &Omega, const real &omega, const real &E,
                 real *partial) {
    const real t2 = sin(E);
    const real t3 = cos(Omega);
    const real t4 = cos(omega);
    const real t6 = sin(Omega);
    const real t7 = cos(inc);
    const real t8 = t6 * t7;
    const real t9 = sin(omega);
    const real t11 = -t3 * t4 + t8 * t9;
    const real t13 = e * e;
    const real t15 = sqrt(0.1e1 - t13);
    const real t17 = cos(E);
    const real t18 = -t17 + e;
    const real t21 = t3 * t9 + t8 * t4;
    const real t25 = t7 * t3;
    const real t27 = t6 * t4 + t9 * t25;
    const real t32 = t6 * t9 - t25 * t4;
    const real t35 = t15 * t9;
    const real t40 = sin(inc);
    const real t42 = sqrt(GM);
    const real t43 = sqrt(a);
    const real t45 = t42 / t43;
    const real t53 = 0.1e1 / (t17 * e - 0.1e1);
    partial[0] = a * (t2 * t11 * t15 + t18 * t21);
    partial[1] = (-t2 * t27 * t15 + t32 * t18) * a;
    partial[2] = -a * (t35 * t2 + t4 * t18) * t40;
    partial[3] = t45 * (-t17 * t11 * t15 - t2 * t21) * t53;
    partial[4] = t45 * (t17 * t27 * t15 - t32 * t2) * t53;
    partial[5] = t45 * t40 * (t4 * t2 + t35 * t17) * t53;
}

void dCartdinc(const real &GM, const real &a, const real &e, const real &inc,
               const real &Omega, const real &omega, const real &E,
               real *partial) {
    const real t2 = sin(Omega);
    const real t4 = cos(omega);
    const real t5 = e * e;
    const real t7 = sqrt(0.1e1 - t5);
    const real t8 = t4 * t7;
    const real t9 = sin(E);
    const real t10 = t8 * t9;
    const real t11 = sin(omega);
    const real t12 = cos(E);
    const real t13 = -t12 + e;
    const real t15 = -t10 + t11 * t13;
    const real t16 = sin(inc);
    const real t20 = cos(Omega);
    const real t25 = cos(inc);
    const real t27 = sqrt(a);
    const real t28 = 0.1e1 / t27;
    const real t29 = sqrt(GM);
    const real t30 = t29 * t28;
    const real t34 = t9 * t11 - t8 * t12;
    const real t38 = 0.1e1 / (t12 * e - 0.1e1);
    partial[0] = -a * t2 * t15 * t16;
    partial[1] = a * t15 * t20 * t16;
    partial[2] = (t10 - t11 * t13) * t25 * a;
    partial[3] = t30 * t2 * t16 * t34 * t38;
    partial[4] = -t29 * t20 * t16 * t34 * t28 * t38;
    partial[5] = t30 * t25 * t34 * t38;
}
