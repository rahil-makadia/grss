#include "utilities.h"

void get_spice_state_lt(const int &spiceID, const real &t0_mjd,
                        const Constants &consts, double state[6], double &lt) {
    real t0_et;
    mjd_to_et(t0_mjd, t0_et);
    SpiceInt center = 0;  // 0 = solar system barycenter
    ConstSpiceChar *frame =
        "J2000";  // Earth mean equator and equinox of J2000,
                  // states output will be ICRF-EME2000 frame
    spkgeo_c(spiceID, t0_et, frame, center, state, &lt);
    // std::cout << "state: " << state[0] << " " << state[1] << " " << state[2]
    // << " " << state[3] << " " << state[4] << " " << state[5] << std::endl;
    // std::cout << "du2m: " << consts.du2m << std::endl;
    // std::cout << "tu2sec: " << consts.tu2sec << std::endl;
    for (int i = 0; i < 6; i++) {
        state[i] *= 1000.0L / consts.du2m;
    }
    for (int i = 3; i < 6; i++) {
        state[i] *= consts.tu2sec;
    }
}

void get_observer_state(const real &tObsMjd,
                        const std::vector<real> &observerInfo,
                        const Constants &consts, const bool &tObsInUTC,
                        std::vector<real> &observerState) {
    SpiceInt baseBody = observerInfo[0];
    if (observerInfo[0] == 500) baseBody = 399;
    if (baseBody == 0) {
        observerState[0] = 0.0L;
        observerState[1] = 0.0L;
        observerState[2] = 0.0L;
        observerState[3] = 0.0L;
        observerState[4] = 0.0L;
        observerState[5] = 0.0L;
        return;
    }
    real t_obs_et;
    real tObsMjdTDB;
    mjd_to_et(tObsMjd, t_obs_et);
    if (tObsInUTC) {
        // std::cout << tObsMjd << " MJD UTC" << t_obs_et << " s UTC -> ";
        SpiceDouble et_minus_utc;
        real sec_past_j2000_utc = t_obs_et;
        deltet_c(sec_past_j2000_utc, "UTC", &et_minus_utc);
        t_obs_et += et_minus_utc;
        et_to_mjd(t_obs_et, tObsMjdTDB);
        // std::cout << t_obs_et << " s ET " << tObsMjdTDB << " MJD TDB" <<
        // std::endl;
    } else {
        tObsMjdTDB = tObsMjd;
    }
    double baseBodyState[6];
    double lt;
    get_spice_state_lt(baseBody, tObsMjdTDB, consts, baseBodyState, lt);

    // real a, f;
    ConstSpiceChar *baseBodyFrame;
    switch ((int)observerInfo[0]) {
        case 399:
            baseBodyFrame = "ITRF93";
            break;
        case 500:
            observerState[0] = (real) baseBodyState[0] + observerInfo[1]/consts.du2m;
            observerState[1] = (real) baseBodyState[1] + observerInfo[2]/consts.du2m;
            observerState[2] = (real) baseBodyState[2] + observerInfo[3]/consts.du2m;
            observerState[3] = (real) baseBodyState[3] + observerInfo[4]/consts.du2m*consts.tu2sec;
            observerState[4] = (real) baseBodyState[4] + observerInfo[5]/consts.du2m*consts.tu2sec;
            observerState[5] = (real) baseBodyState[5] + observerInfo[6]/consts.du2m*consts.tu2sec;
            return;
            break;
        default:
            std::cout << "Given base body: " << baseBody << std::endl;
            throw std::invalid_argument("Given base body not supported");
            break;
    }
    real lon = observerInfo[1];
    real lat = observerInfo[2];
    real rho = observerInfo[3];
    ConstSpiceDouble bodyFixedX = rho * cos(lat) * cos(lon) / 1.0e3L;
    ConstSpiceDouble bodyFixedY = rho * cos(lat) * sin(lon) / 1.0e3L;
    ConstSpiceDouble bodyFixedZ = rho * sin(lat) / 1.0e3L;
    ConstSpiceDouble bodyFixedState[6] = {bodyFixedX, bodyFixedY, bodyFixedZ,
                                          0.0,        0.0,        0.0};
    SpiceDouble observerStateInertial[6];
    SpiceDouble rotMat[6][6];
    sxform_c(baseBodyFrame, "J2000", t_obs_et, rotMat);
    mxvg_c(rotMat, bodyFixedState, 6, 6, observerStateInertial);
    observerStateInertial[0] *= (real)1.0e3L / consts.du2m;
    observerStateInertial[1] *= (real)1.0e3L / consts.du2m;
    observerStateInertial[2] *= (real)1.0e3L / consts.du2m;
    observerStateInertial[3] *= (real)1.0e3L / consts.du2m * consts.tu2sec;
    observerStateInertial[4] *= (real)1.0e3L / consts.du2m * consts.tu2sec;
    observerStateInertial[5] *= (real)1.0e3L / consts.du2m * consts.tu2sec;
    observerState[0] = baseBodyState[0] + observerStateInertial[0];
    observerState[1] = baseBodyState[1] + observerStateInertial[1];
    observerState[2] = baseBodyState[2] + observerStateInertial[2];
    observerState[3] = baseBodyState[3] + observerStateInertial[3];
    observerState[4] = baseBodyState[4] + observerStateInertial[4];
    observerState[5] = baseBodyState[5] + observerStateInertial[5];
}

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

void vnorm(const std::vector<real> &v, real &norm) {
    norm = 0;
    for (size_t i = 0; i < v.size(); i++) {
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

void vcross(const std::vector<real> &v1, const std::vector<real> &v2,
            std::vector<real> &v3) {
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
        std::cout << "utilities.cpp: ERROR: cometary_to_keplerian: Cannot "
                     "handle e = 1 right now!!!"
                  << std::endl;
        exit(1);
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
