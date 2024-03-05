#include "elements.h"

void kepler_solve_elliptic(const real &M, const real &e, real &E, const real &tol,
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
}

void kepler_solve_hyperbolic(const real &M, const real &e, real &EHyp,
                             const real &tol, const int &max_iter) {
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
}

void kepler_solve(const real &epochMjD, const std::vector<real> &cometaryState,
                  const real &GM, real &M, real &E, real &nu,
                  const real &tol, const int &max_iter) {
    const real a = cometaryState[1] / (1 - cometaryState[0]);
    const real a3 = a*a*a;
    const real e = cometaryState[0];
    if (e < 1) {
        const real n = sqrt(GM / a3);
        M = n * (epochMjD - cometaryState[2]);
        kepler_solve_elliptic(M, cometaryState[0], E, tol, max_iter);
        nu = 2 * atan2(tan(E / 2) * sqrt(1 + e), sqrt(1 - e));
    } else if (e > 1) {
        const real n = sqrt(-GM / a3);
        M = n * (epochMjD - cometaryState[2]);
        kepler_solve_hyperbolic(M, cometaryState[0], E, tol, max_iter);
        nu = 2 * atan2(tanh(E / 2) * sqrt(e + 1), sqrt(e - 1));
    } else {
        throw std::runtime_error(
            "utilities.cpp: kepler_solve: Cannot handle e = 1 right "
            "now!!!");
    }
}

void cometary_to_keplerian(const real &epochMjD,
                           const std::vector<real> &cometaryState,
                           std::vector<real> &keplerianState, const real GM) {
    real a = cometaryState[1] / (1 - cometaryState[0]);
    real M, E, nu;
    kepler_solve(epochMjD, cometaryState, GM, M, E, nu);
    keplerianState[0] = a;
    keplerianState[1] = cometaryState[0];
    if (keplerianState[1] < 0) {
        throw std::runtime_error(
            "cometary_to_keplerian: e cannot be negative");
    }
    keplerianState[2] = cometaryState[5];
    keplerianState[3] = cometaryState[3];
    keplerianState[4] = cometaryState[4];
    keplerianState[5] = nu;
    // check that there are no NaNs
    for (size_t i = 0; i < 6; i++) {
        if (std::isnan(keplerianState[i])) {
            std::cout << "cometary_to_keplerian: cometaryState: ";
            for (size_t j = 0; j < 6; j++) {
                std::cout << cometaryState[j] << " ";
            }
            std::cout << std::endl;
            std::cout << "cometary_to_keplerian: keplerianState: ";
            for (size_t j = 0; j < 6; j++) {
                std::cout << keplerianState[j] << " ";
            }
            std::cout << std::endl;
            throw std::runtime_error(
                "cometary_to_keplerian: NaN in keplerian state");
        }
    }
}

void keplerian_to_cometary(const real &epochMjD,
                           const std::vector<real> &keplerianState,
                           std::vector<real> &cometaryState, const real GM) {
    real a = keplerianState[0];
    real e = keplerianState[1];
    if (e < 0) {
        throw std::runtime_error(
            "keplerian_to_cometary: e cannot be negative");
    }
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
    // check that there are no NaNs
    for (size_t i = 0; i < 6; i++) {
        if (std::isnan(cometaryState[i])) {
            std::cout << "keplerian_to_cometary: keplerianState: ";
            for (size_t j = 0; j < 6; j++) {
                std::cout << keplerianState[j] << " ";
            }
            std::cout << std::endl;
            std::cout << "keplerian_to_cometary: cometaryState: ";
            for (size_t j = 0; j < 6; j++) {
                std::cout << cometaryState[j] << " ";
            }
            std::cout << std::endl;
            throw std::runtime_error(
                "keplerian_to_cometary: NaN in cometary state");
        }
    }
}

void keplerian_to_cartesian(const std::vector<real> &keplerianState,
                            std::vector<real> &cartesianState, const real GM) {
    real a = keplerianState[0];
    real e = keplerianState[1];
    if (e < 0) {
        throw std::runtime_error(
            "keplerian_to_cartesian: e cannot be negative");
    }
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
    // check that there are no NaNs
    for (size_t i = 0; i < 6; i++) {
        if (std::isnan(cartesianState[i])) {
            std::cout << "keplerian_to_cartesian: keplerianState: ";
            for (size_t j = 0; j < 6; j++) {
                std::cout << keplerianState[j] << " ";
            }
            std::cout << std::endl;
            std::cout << "keplerian_to_cartesian: cartesianState: ";
            for (size_t j = 0; j < 6; j++) {
                std::cout << cartesianState[j] << " ";
            }
            std::cout << std::endl;
            throw std::runtime_error(
                "keplerian_to_cartesian: NaN in cartesian state");
        }
    }
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
    // check that there are no NaNs
    for (size_t i = 0; i < 6; i++) {
        if (std::isnan(keplerianState[i])) {
            std::cout << "cartesian_to_keplerian: cartesianState: ";
            for (size_t j = 0; j < 6; j++) {
                std::cout << cartesianState[j] << " ";
            }
            std::cout << std::endl;
            std::cout << "cartesian_to_keplerian: keplerianState: ";
            for (size_t j = 0; j < 6; j++) {
                std::cout << keplerianState[j] << " ";
            }
            std::cout << std::endl;
            throw std::runtime_error(
                "cartesian_to_keplerian: NaN in keplerian state");
        }
    }
}

void cometary_to_cartesian(const real &epochMjd,
                           const std::vector<real> &cometaryState,
                           std::vector<real> &cartesianState, const real GM) {
    std::vector<real> keplerianState(6);
    if (cometaryState[0] < 0) {
        throw std::runtime_error(
            "cometary_to_cartesian: e cannot be negative");
    }
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

void get_elements_partials(const real &epochMjd, const std::vector<real> &elems,
                           const std::string conversion,
                           std::vector<std::vector<real>> &partials,
                           const real GM) {
    // real elemTol = 1e-14;
    // only valid conversion is com2cart or kep2cart
    real e, q, om, w, i, a, nu, E, M;
    if (conversion == "com2cart") {
        e = elems[0];
        q = elems[1];
        om = elems[3];
        w = elems[4];
        i = elems[5];
        a = q / (1 - e);
        kepler_solve(epochMjd, elems, GM, M, E, nu);
    } else if (conversion == "kep2cart") {
        a = elems[0];
        e = elems[1];
        q = a * (1 - e);
        i = elems[2];
        om = elems[3];
        w = elems[4];
        nu = elems[5];
        E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(nu / 2));
        M = E - e * sin(E);
    } else {
        throw std::invalid_argument("get_cartesian_partials: invalid conversion "
                                    "type, must be com2cart or kep2cart");
    }
    const std::vector<real> kep = {a, e, i, om, w, nu};
    std::vector<real> cart(6);
    keplerian_to_cartesian(kep, cart, GM);

    real fun, fun1, fun2, den;
    real *partial_fun, *partial_fun1, *partial_fun2, *partial_den;

    real *pos = new real[3];
    pos[0] = cart[0];
    pos[1] = cart[1];
    pos[2] = cart[2];
    real r;
    vnorm(pos, 3, r);
    real *vel = new real[3];
    vel[0] = cart[3];
    vel[1] = cart[4];
    vel[2] = cart[5];
    real v;
    vnorm(vel, 3, v);

    real *h_vec = new real[3];
    vcross(pos, vel, h_vec);
    real h_mag;
    vnorm(h_vec, 3, h_mag);

    real *n_vec = new real[3];
    n_vec[0] = -h_vec[1];
    n_vec[1] = h_vec[0];
    n_vec[2] = 0;
    real n_mag;
    vnorm(n_vec, 3, n_mag);

    real *e_vec = new real[3];
    vcross(vel, h_vec, e_vec);
    for (size_t i = 0; i < 3; i++) {
        e_vec[i] /= GM;
        e_vec[i] -= pos[i] / r;
    }
    real e_mag;
    vnorm(e_vec, 3, e_mag);
    // if (fabs(e_mag - e) > elemTol) {
    //     std::cout << "get_elements_partials: WARNING: e_mag - e = " << e_mag - e
    //               << std::endl;
    // }
    e = e_mag;

    real **partial_r_vec = new real*[6];
    for (size_t i = 0; i < 6; i++) {
        partial_r_vec[i] = new real[3];
        partial_r_vec[i][0] = 0;
        partial_r_vec[i][1] = 0;
        partial_r_vec[i][2] = 0;
    }
    partial_r_vec[0][0] = 1;
    partial_r_vec[1][1] = 1;
    partial_r_vec[2][2] = 1;
    real *partial_r_mag = new real[6];
    for (size_t i = 0; i < 6; i++) {
        if (i > 2) {
            partial_r_mag[i] = 0;
            continue;
        }
        partial_r_mag[i] = pos[i] / r;
    }

    real **partial_v_vec = new real*[6];
    for (size_t i = 0; i < 6; i++) {
        partial_v_vec[i] = new real[3];
        partial_v_vec[i][0] = 0;
        partial_v_vec[i][1] = 0;
        partial_v_vec[i][2] = 0;
    }
    partial_v_vec[3][0] = 1;
    partial_v_vec[4][1] = 1;
    partial_v_vec[5][2] = 1;

    real **partial_h_vec = new real*[6];
    real *partial_h_mag = new real[6];
    for (size_t i = 0; i < 6; i++) {
        partial_h_vec[i] = new real[3];
    }
    for (size_t i = 0; i < 6; i++) {
        vcross(partial_r_vec[i], vel, partial_h_vec[i]);
        real temp[3];
        vcross(pos, partial_v_vec[i], temp);
        for (size_t j = 0; j < 3; j++) {
            partial_h_vec[i][j] += temp[j];
        }
        vdot(partial_h_vec[i], h_vec, 3, partial_h_mag[i]);
        partial_h_mag[i] /= h_mag;
    }

    real **partial_n_vec = new real*[6];
    real *partial_n_mag = new real[6];
    for (size_t i = 0; i < 6; i++) {
        partial_n_vec[i] = new real[3];
    }
    for (size_t i = 0; i < 6; i++) {
        partial_n_vec[i][0] = -partial_h_vec[i][1];
        partial_n_vec[i][1] = partial_h_vec[i][0];
        partial_n_vec[i][2] = 0;
        vdot(partial_n_vec[i], n_vec, 3, partial_n_mag[i]);
        partial_n_mag[i] /= n_mag;
    }

    // eccentricity vector & magnitude
    real **partial_e_vec = new real*[6];
    real *partial_e_mag = new real[6];
    for (size_t i = 0; i < 6; i++) {
        partial_e_vec[i] = new real[3];
    }
    for (size_t i = 0; i < 6; i++) {
        real temp[3];
        vcross(partial_v_vec[i], h_vec, temp);
        for (size_t j = 0; j < 3; j++) {
            partial_e_vec[i][j] = temp[j] / GM;
        }
        vcross(vel, partial_h_vec[i], temp);
        for (size_t j = 0; j < 3; j++) {
            partial_e_vec[i][j] += temp[j] / GM;
        }
        for (size_t j = 0; j < 3; j++) {
            partial_e_vec[i][j] -= partial_r_vec[i][j] / r;
        }
        vdot(pos, partial_r_vec[i], 3, temp[0]);
        for (size_t j = 0; j < 3; j++) {
            partial_e_vec[i][j] += temp[0] * pos[j] / (r * r * r);
        }
        vdot(partial_e_vec[i], e_vec, 3, partial_e_mag[i]);
        partial_e_mag[i] /= e_mag;
    }

    // semi-major axis
    den = 2 / r - v * v / GM;
    // if (fabs(1 / den - a) > elemTol) {
    //     std::cout << "get_elements_partials: WARNING: 1 / den - a = " << 1 / den - a
    //               << std::endl;
    // }
    a = 1 / den;
    partial_den = new real[6];
    partial_den[0] = -2 * pos[0] / (r * r * r);
    partial_den[1] = -2 * pos[1] / (r * r * r);
    partial_den[2] = -2 * pos[2] / (r * r * r);
    partial_den[3] = -2 * vel[0] / GM;
    partial_den[4] = -2 * vel[1] / GM;
    partial_den[5] = -2 * vel[2] / GM;
    real *partial_a = new real[6];
    for (size_t i = 0; i < 6; i++) {
        partial_a[i] = -partial_den[i] / (den * den);
    }

    // perihelion distance
    real *partial_q = new real[6];
    for (size_t i = 0; i < 6; i++) {
        partial_q[i] = partial_a[i] - e * partial_a[i] - a * partial_e_mag[i];
    }
    den = 0;
    delete[] partial_den;

    // true anomaly
    vdot(e_vec, pos, 3, fun);
    fun /= (e_mag * r);
    partial_fun = new real[6];
    const real ecc_times_r = e_mag * r;
    real ecc_dot_r;
    vdot(e_vec, pos, 3, ecc_dot_r);
    real temp1, temp2;
    for (size_t i = 0; i < 6; i++) {
        vdot(partial_e_vec[i], pos, 3, temp1);
        vdot(e_vec, partial_r_vec[i], 3, temp2);
        partial_fun[i] = (ecc_times_r * (temp1 + temp2) -
                          ecc_dot_r * (partial_e_mag[i] * r + e_mag * partial_r_mag[i])) /
                         (ecc_times_r * ecc_times_r);
    }
    real *partial_nu = new real[6];
    for (size_t i = 0; i < 6; i++) {
        partial_nu[i] = -partial_fun[i] / sqrt(1 - fun * fun);
    }
    real posDotVel;
    vdot(pos, vel, 3, posDotVel);
    if (posDotVel < 0) {
        for (size_t i = 0; i < 6; i++) {
            partial_nu[i] *= -1;
        }
    }
    fun = 0;
    delete[] partial_fun;

    // eccentric anomaly
    fun1 = tan(nu / 2);
    fun2 = pow((1 + e) / (1 - e), 0.5);
    partial_fun1 = new real[6];
    partial_fun2 = new real[6];
    for (size_t i = 0; i < 6; i++) {
        partial_fun1[i] = partial_nu[i] / (2 * cos(nu / 2) * cos(nu / 2));
        partial_fun2[i] = (sqrt(1 - e) * partial_e_mag[i] / sqrt(1 + e) / 2 -
                           sqrt(1 + e) * -partial_e_mag[i] / sqrt(1 - e) / 2) /
                          (1 - e);
    }
    fun = fun1 / fun2;
    partial_fun = new real[6];
    for (size_t i = 0; i < 6; i++) {
        partial_fun[i] = (fun2 * partial_fun1[i] - fun1 * partial_fun2[i]) / (fun2 * fun2);
    }
    real *partial_E = new real[6];
    for (size_t i = 0; i < 6; i++) {
        partial_E[i] = 2 * partial_fun[i] / (1 + fun * fun);
    }
    fun = 0;
    fun1 = 0;
    fun2 = 0;
    delete[] partial_fun1;
    delete[] partial_fun2;
    delete[] partial_fun;
    /*
    alternative for eccentric anomaly
    fun = 1 / e - r / (a * e);
    partial_fun = new real[6];
    for (size_t i = 0; i < 6; i++) {
        partial_fun[i] = -partial_e_mag[i] / (e * e) -
                         ((a * e * partial_r_mag[i] - r * (partial_a[i] * e + a * partial_e_mag[i])) /
                          (a * a * e * e));
    }
    partial_E = new real[6];
    for (size_t i = 0; i < 6; i++) {
        partial_E[i] = -partial_fun[i] / sqrt(1 - fun * fun);
    }
    fun = 0;
    delete[] partial_fun;
    */ 

    // mean anomaly
    // if (fabs(E - e * sin(E) - M) > elemTol) {
    //     std::cout << "get_elements_partials: WARNING: E - e*sin(E) - M = " << E - e * sin(E) - M
    //               << std::endl;
    // }
    M = E - e * sin(E);
    real *partial_M = new real[6];
    for (size_t i = 0; i < 6; i++) {
        partial_M[i] = partial_E[i] - e * cos(E) * partial_E[i] - partial_e_mag[i] * sin(E);
    }

    // time of periapsis passage
    real *partial_tp = new real[6];
    for (size_t i = 0; i < 6; i++) {
        partial_tp[i] = -M / sqrt(GM) * 3 * sqrt(a) * partial_a[i] / 2 - partial_M[i] * sqrt(a * a * a / GM);
    }

    // longitude of ascending node
    fun = -h_vec[1] / n_mag;
    partial_fun = new real[6];
    for (size_t i = 0; i < 6; i++) {
        partial_fun[i] = (n_mag * -partial_h_vec[i][1] - -h_vec[1] * partial_n_mag[i]) / (n_mag * n_mag);
    }
    real *partial_Omega = new real[6];
    for (size_t i = 0; i < 6; i++) {
        partial_Omega[i] = -partial_fun[i] / sqrt(1 - fun * fun);
    }
    if (n_vec[1] < 0) {
        for (size_t i = 0; i < 6; i++) {
            partial_Omega[i] *= -1;
        }
    }
    fun = 0;
    delete[] partial_fun;

    // argument of periapsis
    vdot(n_vec, e_vec, 3, fun);
    fun /= (n_mag * e_mag);
    partial_fun = new real[6];
    real n_times_e = n_mag * e_mag;
    real n_dot_e;
    vdot(n_vec, e_vec, 3, n_dot_e);
    for (size_t i = 0; i < 6; i++) {
        real temp1, temp2;
        vdot(partial_n_vec[i], e_vec, 3, temp1);
        vdot(n_vec, partial_e_vec[i], 3, temp2);
        partial_fun[i] = (n_times_e * (temp1 + temp2) -
                          n_dot_e * (partial_n_mag[i] * e_mag + n_mag * partial_e_mag[i])) /
                         (n_times_e * n_times_e);
    }
    real *partial_w = new real[6];
    for (size_t i = 0; i < 6; i++) {
        partial_w[i] = -partial_fun[i] / sqrt(1 - fun * fun);
    }
    if (e_vec[2] < 0) {
        for (size_t i = 0; i < 6; i++) {
            partial_w[i] *= -1;
        }
    }
    fun = 0;
    delete[] partial_fun;

    // inclination
    fun = h_vec[2] / h_mag;
    partial_fun = new real[6];
    for (size_t i = 0; i < 6; i++) {
        partial_fun[i] = (h_mag * partial_h_vec[i][2] - h_vec[2] * partial_h_mag[i]) / (h_mag * h_mag);
    }
    real *partial_i = new real[6];
    for (size_t i = 0; i < 6; i++) {
        partial_i[i] = -partial_fun[i] / sqrt(1 - fun * fun);
    }

    // assign partials
    if (conversion == "com2cart") {
        for (size_t i = 0; i < 6; i++) {
            partials[0][i] = partial_e_mag[i];
            partials[1][i] = partial_q[i];
            partials[2][i] = partial_tp[i];
            partials[3][i] = partial_Omega[i];
            partials[4][i] = partial_w[i];
            partials[5][i] = partial_i[i];
        }
    } else {
        for (size_t i = 0; i < 6; i++) {
            partials[0][i] = partial_a[i];
            partials[1][i] = partial_e_mag[i];
            partials[2][i] = partial_i[i];
            partials[3][i] = partial_Omega[i];
            partials[4][i] = partial_w[i];
            partials[5][i] = partial_nu[i];
        }
    }

    // clean up
    delete[] pos;
    delete[] vel;
    delete[] h_vec;
    delete[] n_vec;
    delete[] e_vec;
    for (size_t i = 0; i < 6; i++) {
        delete[] partial_r_vec[i];
        delete[] partial_v_vec[i];
        delete[] partial_h_vec[i];
        delete[] partial_n_vec[i];
        delete[] partial_e_vec[i];
    }
    delete[] partial_r_vec;
    delete[] partial_v_vec;
    delete[] partial_h_vec;
    delete[] partial_n_vec;
    delete[] partial_e_vec;
    delete[] partial_r_mag;
    delete[] partial_h_mag;
    delete[] partial_n_mag;
    delete[] partial_e_mag;
    delete[] partial_a;
    delete[] partial_q;
    delete[] partial_nu;
    delete[] partial_E;
    delete[] partial_M;
    delete[] partial_tp;
    delete[] partial_Omega;
    delete[] partial_w;
    delete[] partial_i;
}

void get_cartesian_partials(const real &epochMjd,
                            const std::vector<real> &elems,
                            const std::string &conversion,
                            std::vector<std::vector<real>> &partials,
                            const real GM) {
    std::vector<std::vector<real> > elemPartials(6, std::vector<real>(6));
    get_elements_partials(epochMjd, elems, conversion, elemPartials, GM);
    mat_inv(elemPartials, partials);
}
