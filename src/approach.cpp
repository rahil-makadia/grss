/**
 * @file    approach.cpp
 * @brief   Source file for close approach and impact detection.
 * @author  Rahil Makadia <makadia2@illinois.edu>
 *
 * @section     LICENSE
 * Copyright (C) 2022-2025 Rahil Makadia
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, see <https://www.gnu.org/licenses>.
 */

#include "approach.h"

/**
 * @brief Get the atmosphere thickness offset for the central body.
 * 
 * @param[in] centralBodySpiceId SPICE ID of the central body.
 * @return atm_offset Offset for the atmosphere of the central body.
 */
static real get_atm_offset(const int &centralBodySpiceId){
    real atm_offset = 0;
    real au2km = 149597870.7;
    switch (centralBodySpiceId)
    {
    case 399:
        atm_offset = 100.0/au2km;
        break;
    default:
        atm_offset = 0;
        break;
    }
    return atm_offset;
}

/**
 * @brief Convert body-fixed rectangular coordinates to geodetic coordinates for a spheroid.
 * 
 * @param[in] x X-coordinate in body-fixed frame.
 * @param[in] y Y-coordinate in body-fixed frame.
 * @param[in] z Z-coordinate in body-fixed frame.
 * @param[out] lon Longitude in geodetic frame.
 * @param[out] lat Latitude in geodetic frame.
 * @param[out] h Height above spheroid in geodetic frame.
 */
static void rec_to_geodetic(const real &x, const real &y, const real &z,
                            real &lon, real &lat, real &h) {
    const real a = EARTH_RAD_WGS84;
    const real f = EARTH_FLAT_WGS84;
    if (x == 0.0 && y == 0.0) {
        lon = 0.0;
    } else {
        lon = atan2(y, x);
        if (lon < 0.0) {
            lon += 2 * PI;
        }
    }
    const real xy = sqrt(x*x + y*y);
    const real e2 = 2*f - f*f;
    if (xy == 0.0) {
        if (z >= 0.0) {
            lat = PI/2;
        } else {
            lat = -PI/2;
        }
        h = fabs(z) - a*sqrt(1-e2);
        return;
    }
    lat = atan(z/(xy*(1-e2)));
    real N;
    for (size_t k = 0; k < 5; k++) {
        N = a/sqrt(1-e2*sin(lat)*sin(lat));
        h = xy/cos(lat) - N;
        lat = atan(z*(N+h)/(xy*(N*(1-e2)+h)));
    }
}

/**
 * @param[inout] propSim PropSimulation object for the integration.
 * @param[in] tOld Time at the previous integrator epoch.
 * @param[in] xIntegOld State at the previous integrator epoch.
 * @param[in] t Time at the current integrator epoch.
 * @param[in] xInteg State at the current integrator epoch.
 */
void check_ca_or_impact(PropSimulation *propSim, const real &tOld,
                        const std::vector<real> xIntegOld, const real &t,
                        const std::vector<real> xInteg) {
    // Check for close approach or impact
    // If a close approach is detected, the time of closest approach is
    // determined by root finding using Brent's method.
    // If an impact is detected, the simulation is terminated.
    const bool forwardProp = propSim->integParams.t0 < propSim->integParams.tf;
    const bool backwardProp = propSim->integParams.t0 > propSim->integParams.tf;
    size_t starti = 0;
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        if (!propSim->integBodies[i].logCA) {
            continue;
        }
        size_t startj = 0;
        for (size_t j = 0; j < propSim->integParams.nTotal; j++) {
            if (i != j) {
                real caTol, relPosOld[3], relPos[3], relVelOld[3], relVel[3];
                if (j < propSim->integParams.nInteg) {
                    caTol = propSim->integBodies[j].caTol;
                    for (size_t k = 0; k < 3; k++) {
                        relPosOld[k] = xIntegOld[starti + k] - xIntegOld[startj + k];
                        relPos[k] = xInteg[starti + k] - xInteg[startj + k];
                        relVelOld[k] = xIntegOld[starti + 3 + k] - xIntegOld[startj + 3 + k];
                        relVel[k] = xInteg[starti + 3 + k] - xInteg[startj + 3 + k];
                    }
                } else {
                    int spiceId = propSim->spiceBodies[j-propSim->integParams.nInteg].spiceId;
                    caTol = propSim->spiceBodies[j-propSim->integParams.nInteg].caTol;
                    double xSpice[9], xSpiceOld[9];
                    get_spk_state(spiceId, t, propSim->spkEphem, xSpice);
                    get_spk_state(spiceId, tOld, propSim->spkEphem, xSpiceOld);
                    for (size_t k = 0; k < 3; k++) {
                        relPosOld[k] = xIntegOld[starti + k] - xSpiceOld[k];
                        relPos[k] = xInteg[starti + k] - xSpice[k];
                        relVelOld[k] = xIntegOld[starti + 3 + k] - xSpiceOld[3 + k];
                        relVel[k] = xInteg[starti + 3 + k] - xSpice[3 + k];
                    }
                }
                const real relDistOld = sqrt(relPosOld[0] * relPosOld[0] +
                                             relPosOld[1] * relPosOld[1] +
                                             relPosOld[2] * relPosOld[2]);
                const real relDist =
                    sqrt(relPos[0] * relPos[0] + relPos[1] * relPos[1] +
                         relPos[2] * relPos[2]);
                // check close approach
                const real radialVelOld =
                    (relPosOld[0] * relVelOld[0] + relPosOld[1] * relVelOld[1] +
                     relPosOld[2] * relVelOld[2]) /
                    relDistOld;
                const real radialVel =
                    (relPos[0] * relVel[0] + relPos[1] * relVel[1] +
                     relPos[2] * relVel[2]) /
                    relDist;
                const bool relDistMinimum = (forwardProp && radialVelOld < 0.0 && radialVel >= 0.0) ||
                                            (backwardProp && radialVelOld > 0.0 && radialVel <= 0.0);
                const bool relDistWithinTol = relDist <= caTol || relDistOld <= caTol;
                if (relDistMinimum && relDistWithinTol) {
                    CloseApproachParameters ca;
                    get_ca_or_impact_time(propSim, i, j, tOld, t, ca.t, ca_rdot_calc);
                    ca.xRel = get_rel_state(propSim, i, j, ca.t);
                    ca.tCA = ca.t;
                    ca.xRelCA = ca.xRel;
                    ca.flybyBodyIdx = i;
                    ca.centralBodyIdx = j;
                    ca.get_ca_parameters(propSim, ca.t);
                    propSim->caParams.push_back(ca);
                    if (ca.impact){
                        real rCA;
                        impact_r_calc(propSim, i, j, ca.t, rCA);
                        // near miss check
                        if (rCA > 0){
                            ca.impact = false;
                            break;
                        }
                        ImpactParameters impact;
                        const real maxImpactCaTDiff = 600.0/86400.0; // 10 minutes
                        real tImpStart;
                        if (forwardProp){
                            tImpStart = fmax(tOld-maxImpactCaTDiff, propSim->integParams.t0-propSim->tEvalMargin);
                        } else {
                            tImpStart = fmin(tOld+maxImpactCaTDiff, propSim->integParams.t0+propSim->tEvalMargin);
                        }
                        get_ca_or_impact_time(propSim, i, j, tImpStart, ca.t, impact.t, impact_r_calc);
                        impact.xRel = get_rel_state(propSim, i, j, impact.t);
                        impact.tCA = ca.t;
                        impact.xRelCA = ca.xRel;
                        impact.flybyBodyIdx = i;
                        impact.centralBodyIdx = j;
                        impact.get_ca_parameters(propSim, impact.t);
                        impact.get_impact_parameters(propSim);
                        propSim->impactParams.push_back(impact);
                        propSim->integBodies[i].logCA = false;
                    }
                }
            }
            if (j < propSim->integParams.nInteg) {
                startj += 2*propSim->integBodies[j].n2Derivs;
            }
        }
        starti += 2*propSim->integBodies[i].n2Derivs;
    }
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] i Index of the first body.
 * @param[in] j Index of the second body.
 * @param[in] t Time at the current integrator epoch.
 * @param[out] rDot Relative radial velocity.
 */
void ca_rdot_calc(PropSimulation *propSim, const size_t &i, const size_t &j,
                  const real &t, real &rDot) {
    // Calculate the radial velocity between two bodies at a given time
    // This is used in get_ca_or_impact_time to find the time of closest approach
    std::vector<real> xRel = get_rel_state(propSim, i, j, t);
    real relDist =
        sqrt(xRel[0] * xRel[0] + xRel[1] * xRel[1] + xRel[2] * xRel[2]);
    rDot =
        (xRel[0] * xRel[3] + xRel[1] * xRel[4] + xRel[2] * xRel[5]) / relDist;
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] i Index of the first body.
 * @param[in] j Index of the second body.
 * @param[in] t Time at the current integrator epoch.
 * @param[out] r Relative distance.
 */
void impact_r_calc(PropSimulation *propSim, const size_t &i, const size_t &j,
                  const real &t, real &r) {
    // Calculate the distance between two bodies at a given time, accounting for
    // the radius of the bodies
    std::vector<real> xRel = get_rel_state(propSim, i, j, t);
    if (j < propSim->integParams.nInteg) {
        const real relDist =
            sqrt(xRel[0] * xRel[0] + xRel[1] * xRel[1] + xRel[2] * xRel[2]);
        r = relDist - propSim->integBodies[i].radius - propSim->integBodies[j].radius;
        return;
    } else {
        std::string baseBodyFrame;
        const int spiceId = propSim->spiceBodies[j - propSim->integParams.nInteg].spiceId;
        const real atmOffset = get_atm_offset(spiceId);
        if (spiceId == 399){
            get_baseBodyFrame(spiceId, t, baseBodyFrame);
            std::vector<std::vector<real>> rotMat(6, std::vector<real>(6));
            get_pck_rotMat("J2000", baseBodyFrame, t, propSim->pckEphem, rotMat);
            std::vector<real> xBody(6);
            mat_vec_mul(rotMat, xRel, xBody);
            const real x = xBody[0];
            const real y = xBody[1];
            const real z = xBody[2];
            real lon, lat, h;
            rec_to_geodetic(x, y, z, lon, lat, h);
            r = h - atmOffset;
            return;
        } else {
            const real relDist =
                sqrt(xRel[0] * xRel[0] + xRel[1] * xRel[1] + xRel[2] * xRel[2]);
            r = relDist - propSim->integBodies[i].radius - propSim->spiceBodies[j - propSim->integParams.nInteg].radius - atmOffset;
            return;
        }
    }
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] i Index of the first body.
 * @param[in] j Index of the second body.
 * @param[in] t Time to compute the relative state.
 * @return xRel Relative state of the body.
 */
std::vector<real> get_rel_state(PropSimulation *propSim, const size_t &i,
                                const size_t &j, const real &t) {
    std::vector<real> xInterp = propSim->interpolate(t);
    std::vector<real> xRel(2*propSim->integBodies[i].n2Derivs, std::numeric_limits<real>::quiet_NaN());
    size_t starti = 0;
    for (size_t k = 0; k < i; k++) {
        starti += 2*propSim->integBodies[k].n2Derivs;
    }
    for (size_t k = 0; k < 2*propSim->integBodies[i].n2Derivs; k++) {
        xRel[k] = xInterp[starti + k];
    }
    if (j < propSim->integParams.nInteg) {
        size_t startj = 0;
        for (size_t k = 0; k < j; k++) {
            startj += 2*propSim->integBodies[k].n2Derivs;
        }
        for (size_t k = 0; k < 6; k++) {
            xRel[k] -= xInterp[startj + k];
        }
    } else {
        double xSpice[9];
        get_spk_state(
            propSim->spiceBodies[j - propSim->integParams.nInteg].spiceId, t,
            propSim->spkEphem, xSpice);
        for (size_t k = 0; k < 6; k++) {
            xRel[k] -= xSpice[k];
        }
    }
    return xRel;
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] i Index of the first body.
 * @param[in] j Index of the second body.
 * @param[in] x1 Initial time of the bracketed interval.
 * @param[in] x2 Final time of the bracketed interval.
 * @param[out] tCA Time of close approach or impact.
 * @param[in] zero_func Function to compute the zero of (rDot for CA or r for impacts).
 */
void get_ca_or_impact_time(PropSimulation *propSim, const size_t &i,
                           const size_t &j, const real &x1, const real &x2,
                           real &tCA,
                           void (*zero_func)(PropSimulation *, const size_t &,
                                             const size_t &, const real &,
                                             real &)) {
    // Brent's method for root finding, from Numerical Recipes in C/C++, 3rd
    // edition, p. 454
    // make sure x1 and x2 are not nan
    if (std::isnan(x1) || std::isnan(x2)) {
        std::cout << "x1: " << x1 << " x2: " << x2 << std::endl;
        throw std::runtime_error("get_ca_or_impact_time: One of the interval points is nan.");
    }
    const real tol = 1.0e-6 / 86400.0;  // 1 microsecond
    const size_t maxIter = 100;
    const real eps = std::numeric_limits<real>::epsilon();
    real a = x1;
    real b = x2;
    real c = x2;
    real d = 0.0;
    real e = 0.0;
    real min1, min2;
    real fa, fb, fc, p, q, r, s, tol1, xm;
    zero_func(propSim, i, j, a, fa);
    zero_func(propSim, i, j, b, fb);
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
        std::cout << "t1: " << x1 << " t2: " << x2 << " f1: " << fa << " f2: " << fb << std::endl;
        throw std::runtime_error("get_ca_or_impact_time: Root must be bracketed.");
    }
    fc = fb;
    size_t iter;
    for (iter = 0; iter < maxIter; ++iter) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a;
            fc = fa;
            e = d = b - a;
        }
        if (fabs(fc) < fabs(fb)) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        tol1 = 2.0 * eps * fabs(b) + 0.5 * tol;
        xm = 0.5 * (c - b);
        if (fabs(xm) <= tol1 || fb == 0.0) {
            tCA = b;
            return;
        }
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            s = fb / fa;
            if (a == c) {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0) {
                q = -q;
            }
            p = fabs(p);
            min1 = 3.0 * xm * q - fabs(tol1 * q);
            min2 = fabs(e * q);
            if (2.0 * p < (min1 < min2 ? min1 : min2)) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        if (fabs(d) > tol1) {
            b += d;
        } else {
            b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
        }
        zero_func(propSim, i, j, b, fb);
    }
    std::cout << "WARNING: Maximum number of iterations exceeded in "
                 "get_ca_or_impact_time!!! Impact/CA time may not be accurate."
                << std::endl;
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[inout] ca CloseApproachParameters object for the close approach.
 * @param[in] mu Gravitational parameter of the central body.
 * @param[in] radius Radius of the central body.
 */
void get_bplane_partials(PropSimulation *propSim, CloseApproachParameters *ca,
                         const real &mu, const real &radius) {
    const size_t i = ca->flybyBodyIdx;
    const size_t j = ca->centralBodyIdx;
    real h, pos[3], vel[3], hVec[3];
    for (size_t k = 0; k < 3; k++) {
        pos[k] = ca->xRelCA[k];
        vel[k] = ca->xRelCA[k + 3];
    }
    real r, v, posDotVel;
    vnorm(pos, 3, r);
    vnorm(vel, 3, v);
    vdot(pos, vel, 3, posDotVel);
    vcross(pos, vel, hVec);
    vnorm(hVec, 3, h);
    real vCrossH[3], eVec[3], e;
    vcross(vel, hVec, vCrossH);
    for (size_t k = 0; k < 3; k++) {
        eVec[k] = vCrossH[k] / mu - pos[k] / r;
    }
    vnorm(eVec, 3, e);
    const real a = (mu * r) / (2 * mu - r * v * v);
    const real vInf = sqrt(-mu / a);
    const real n = sqrt(-mu / a / a / a);
    const real sinhF = vInf*posDotVel/mu/e;
    const real F = asinh(sinhF); // or F = -log(2*r/a/e);

    real **partial_r_vec = new real*[6];
    real **partial_v_vec = new real*[6];
    for (size_t k = 0; k < 6; k++) {
        partial_r_vec[k] = new real[3];
        partial_r_vec[k][0] = 0;
        partial_r_vec[k][1] = 0;
        partial_r_vec[k][2] = 0;
        partial_v_vec[k] = new real[3];
        partial_v_vec[k][0] = 0;
        partial_v_vec[k][1] = 0;
        partial_v_vec[k][2] = 0;
    }
    partial_r_vec[0][0] = 1;
    partial_r_vec[1][1] = 1;
    partial_r_vec[2][2] = 1;
    partial_v_vec[3][0] = 1;
    partial_v_vec[4][1] = 1;
    partial_v_vec[5][2] = 1;
    real *partial_alpha = new real[6];
    real *partial_vInf = new real[6];
    real *partial_a = new real[6];
    real *partial_n = new real[6];
    real temp1, temp2;
    for (size_t k = 0; k < 6; k++) {
        vdot(vel, partial_v_vec[k], 3, temp1);
        vdot(pos, partial_r_vec[k], 3, temp2);
        partial_alpha[k] = 2*(temp1 + mu*temp2/r/r/r);
        partial_vInf[k] = partial_alpha[k]/2/vInf;
        partial_a[k] = -a*partial_alpha[k]/(v*v - 2*mu/r);
        partial_n[k] = -3*n*partial_a[k]/2/a;
    }
    real **partial_hVec = new real*[6];
    real *partial_h = new real[6];
    for (size_t k = 0; k < 6; k++) {
        partial_hVec[k] = new real[3];
    }
    real temp1Vec[3];
    for (size_t k = 0; k < 6; k++) {
        vcross(partial_r_vec[k], vel, partial_hVec[k]);
        vcross(pos, partial_v_vec[k], temp1Vec);
        for (size_t k2 = 0; k2 < 3; k2++) {
            partial_hVec[k][k2] += temp1Vec[k2];
        }
        vdot(partial_hVec[k], hVec, 3, partial_h[k]);
        partial_h[k] /= h;
    }
    real **partial_eVec = new real*[6];
    real *partial_e = new real[6];
    for (size_t k = 0; k < 6; k++) {
        partial_eVec[k] = new real[3];
        vcross(partial_v_vec[k], hVec, temp1Vec);
        for (size_t k2 = 0; k2 < 3; k2++) {
            partial_eVec[k][k2] = temp1Vec[k2] / mu;
        }
        vcross(vel, partial_hVec[k], temp1Vec);
        for (size_t k2 = 0; k2 < 3; k2++) {
            partial_eVec[k][k2] += temp1Vec[k2] / mu;
        }
        for (size_t k2 = 0; k2 < 3; k2++) {
            partial_eVec[k][k2] -= partial_r_vec[k][k2] / r;
        }
        vdot(pos, partial_r_vec[k], 3, temp1Vec[0]);
        for (size_t k2 = 0; k2 < 3; k2++) {
            partial_eVec[k][k2] += temp1Vec[0] * pos[k2] / r/r/r;
        }
        vdot(partial_eVec[k], eVec, 3, partial_e[k]);
        partial_e[k] /= e;
    }
    real *partial_F = new real[6];
    for (size_t k = 0; k < 6; k++) {
        vdot(pos, partial_v_vec[k], 3, temp1);
        vdot(vel, partial_r_vec[k], 3, temp2);
        partial_F[k] = (posDotVel*partial_vInf[k]/mu + vInf*(temp1+temp2)/mu - partial_e[k]*sinhF)/e/cosh(F);
    }
    std::vector<real> accMap(propSim->integParams.n2Derivs, 0.0);
    std::vector<real> xMap = propSim->interpolate(ca->tMap);
    get_state_der(propSim, ca->tMap, xMap, accMap);
    real accRelMap[3], accPlanet[3];
    if (j < propSim->integParams.nInteg) {
        accPlanet[0] = propSim->integBodies[j].acc[0];
        accPlanet[1] = propSim->integBodies[j].acc[1];
        accPlanet[2] = propSim->integBodies[j].acc[2];
    } else {
        accPlanet[0] = propSim->spiceBodies[j - propSim->integParams.nInteg].acc[0];
        accPlanet[1] = propSim->spiceBodies[j - propSim->integParams.nInteg].acc[1];
        accPlanet[2] = propSim->spiceBodies[j - propSim->integParams.nInteg].acc[2];
    }
    for (size_t k = 0; k < 3; k++) {
        accRelMap[k] = propSim->integBodies[i].acc[k] - accPlanet[k];
    }
    real posMapDotAccMap, vMap;
    real posMap[3] = {ca->xRelMap[0], ca->xRelMap[1], ca->xRelMap[2]};
    real velMap[3] = {ca->xRelMap[3], ca->xRelMap[4], ca->xRelMap[5]};
    vnorm(velMap, 3, vMap);
    vdot(posMap, accRelMap, 3, posMapDotAccMap);
    std::vector<real> dtLinMinustCA(6, 0.0);
    const real tLin = ca->tMap + (F - e*sinhF)/n - log(e)/n;
    for (size_t k = 0; k < 6; k++) {
        dtLinMinustCA[k] = (r*partial_F[k]/a - (e*sinhF + 1)*partial_e[k]/e - (tLin-ca->tCA)*partial_n[k])/n;
        vdot(velMap, partial_r_vec[k], 3, temp2);
        vdot(posMap, partial_v_vec[k], 3, temp1);
        ca->dt[k] = -(temp1+temp2)/(posMapDotAccMap + vMap*vMap);
    }
    real **partial_xCA = new real*[6];
    for (size_t k = 0; k < 6; k++) {
        partial_xCA[k] = new real[6];
        for (size_t k2 = 0; k2 < 6; k2++) {
            partial_xCA[k][k2] = 0;
        }
    }
    for (size_t k = 0; k < 3; k++) {
        for (size_t k2 = 0; k2 < 6; k2++) {
            partial_xCA[k][k2] = partial_r_vec[k2][k] + vel[k]*ca->dt[k2];
            partial_xCA[k+3][k2] = partial_v_vec[k2][k] + accRelMap[k]*ca->dt[k2];
        }
    }
    real pHat[3], qHat[3], hCrossPHat[3];
    vunit(eVec, 3, pHat);
    vcross(hVec, pHat, hCrossPHat);
    for (size_t k = 0; k < 3; k++) {
        qHat[k] = hCrossPHat[k] / h;
    }
    real rHat[3], sHat[3], tHat[3], vCrossSHatVec[3];
    for (size_t k = 0; k < 3; k++) {
        sHat[k] = pHat[k]/e + sqrt(e * e - 1) * qHat[k]/e;
    }
    real vHat[3] = {0.0, 0.0, -1.0};
    vcross(vHat, sHat, vCrossSHatVec);
    vunit(vCrossSHatVec, 3, tHat);
    vcross(sHat, tHat, rHat);
    real **partial_pHat = new real*[6];
    real **partial_qHat = new real*[6];
    real **partial_rHat = new real*[6];
    real **partial_sHat = new real*[6];
    real **partial_tHat = new real*[6];
    real temp2Vec[3];
    for (size_t k = 0; k < 6; k++) {
        partial_pHat[k] = new real[3];
        partial_qHat[k] = new real[3];
        partial_rHat[k] = new real[3];
        partial_sHat[k] = new real[3];
        partial_tHat[k] = new real[3];
        for (size_t k2 = 0; k2 < 3; k2++) {
            partial_pHat[k][k2] = partial_eVec[k][k2]/e - partial_e[k]*eVec[k2]/e/e;
        }
        vcross(partial_hVec[k], pHat, temp1Vec);
        vcross(hVec, partial_pHat[k], temp2Vec);
        for (size_t k2 = 0; k2 < 3; k2++) {
            partial_qHat[k][k2] = (temp1Vec[k2] + temp2Vec[k2] - partial_h[k]*qHat[k2])/h;
        }
        for (size_t k2 = 0; k2 < 3; k2++) {
            partial_sHat[k][k2] = -partial_e[k] * pHat[k2] / e / e +
                                    partial_pHat[k][k2] / e +
                                    partial_e[k] * qHat[k2] / e / e / sqrt(e * e - 1) +
                                    sqrt(e * e - 1) * partial_qHat[k][k2] / e;
        }
        vcross(vHat, partial_sHat[k], temp1Vec);
        vdot(tHat, temp1Vec, 3, temp1);
        vcross(vHat, sHat, temp2Vec);
        vnorm(temp2Vec, 3, temp2);
        for (size_t k2 = 0; k2 < 3; k2++) {
            partial_tHat[k][k2] = (temp1Vec[k2] - temp1*tHat[k2])/temp2;
        }
        vcross(partial_sHat[k], tHat, temp1Vec);
        vcross(sHat, partial_tHat[k], temp2Vec);
        for (size_t k2 = 0; k2 < 3; k2++) {
            partial_rHat[k][k2] = temp1Vec[k2] + temp2Vec[k2];
        }
    }
    real lambda = sqrt(1.0 + 2 * mu / radius / vInf / vInf);
    real *partial_lambda = new real[6];
    std::vector<real> k_dx(6, 0.0);
    std::vector<real> k_dy(6, 0.0);
    std::vector<real> s_dx(6, 0.0);
    std::vector<real> s_dy(6, 0.0);
    for (size_t k = 0; k < 6; k++) {
        vdot(tHat, partial_hVec[k], 3, temp1);
        vdot(hVec, partial_tHat[k], 3, temp2);
        k_dx[k] = (temp1 + temp2 - ca->kizner.x*partial_vInf[k])/vInf;
        vdot(rHat, partial_hVec[k], 3, temp1);
        vdot(hVec, partial_rHat[k], 3, temp2);
        k_dy[k] = -(temp1 + temp2 + ca->kizner.y*partial_vInf[k])/vInf;
        partial_lambda[k] = -2*mu*partial_vInf[k]/lambda/radius/vInf/vInf/vInf;
        s_dx[k] = (k_dx[k] - ca->scaled.x*partial_lambda[k])/lambda;
        s_dy[k] = (k_dy[k] - ca->scaled.y*partial_lambda[k])/lambda;
    }
    // map off-nominal trajectory for full partials
    vec_mat_mul(k_dx, partial_xCA, 6, ca->kizner.dx);
    vec_mat_mul(k_dy, partial_xCA, 6, ca->kizner.dy);
    vec_mat_mul(s_dx, partial_xCA, 6, ca->scaled.dx);
    vec_mat_mul(s_dy, partial_xCA, 6, ca->scaled.dy);
    vec_mat_mul(dtLinMinustCA, partial_xCA, 6, ca->dtLin);
    for (size_t k = 0; k < 6; k++) {
        ca->dtLin[k] += ca->dt[k];
    }
    // Acceleration of the Sun for opik formulation
    real accSun[3];
    for (size_t k = 0; k < propSim->integParams.nSpice; k++) {
        if (propSim->spiceBodies[k].spiceId == 10) {
            accSun[0] = propSim->spiceBodies[k].acc[0];
            accSun[1] = propSim->spiceBodies[k].acc[1];
            accSun[2] = propSim->spiceBodies[k].acc[2];
            break;
        }
    }
    real **partial_vel_planet = new real*[6];
    for (size_t k = 0; k < 6; k++) {
        partial_vel_planet[k] = new real[3];
        for (size_t k2 = 0; k2 < 3; k2++) {
            partial_vel_planet[k][k2] = ca->dt[k]*(accPlanet[k2] - accSun[k2]);
        }
    }
    // partials of xi and zeta w.r.t planet velocity are needed for total derivative
    real **partial_vpl_vpl = new real*[3];
    for (size_t k = 0; k < 3; k++) {
        partial_vpl_vpl[k] = new real[3];
        for (size_t k2 = 0; k2 < 3; k2++) {
            partial_vpl_vpl[k][k2] = 0;
        }
    }
    partial_vpl_vpl[0][0] = 1;
    partial_vpl_vpl[1][1] = 1;
    partial_vpl_vpl[2][2] = 1;

    double xCentralBody[9];
    get_spk_state(ca->centralBodySpiceId, ca->tCA, propSim->spkEphem, xCentralBody);
    double xSun[9];
    get_spk_state(10, ca->tCA, propSim->spkEphem, xSun);
    real vCentralBodyHelio[3];
    for (size_t k = 0; k < 3; k++) {
        vCentralBodyHelio[k] = xCentralBody[3+k]-xSun[3+k];
    }
    real xiHat[3], zetaHat[3], vPlanetCrossSHatVec[3];
    vcross(vCentralBodyHelio, sHat, vPlanetCrossSHatVec);
    vunit(vPlanetCrossSHatVec, 3, xiHat);
    vcross(sHat, xiHat, zetaHat);
    for (size_t k = 0; k < 3; k++) {
        zetaHat[k] *= -1;
    }
    real vPlanetCrossSHat;
    vnorm(vPlanetCrossSHatVec, 3, vPlanetCrossSHat);
    real **partial_vPlanetCrossSHatVec_vpl = new real*[3];
    for (size_t k = 0; k < 3; k++) {
        partial_vPlanetCrossSHatVec_vpl[k] = new real[3];
        vcross(partial_vpl_vpl[k], sHat, partial_vPlanetCrossSHatVec_vpl[k]);
    }
    real *partial_vPlanetCrossSHat_vpl = new real[3];
    vcross(xiHat, sHat, partial_vPlanetCrossSHat_vpl);
    for (size_t k = 0; k < 3; k++){
        partial_vPlanetCrossSHat_vpl[k] *= -1.0;
    }
    real **partial_xiHat_vpl = new real*[3];
    real **partial_zetaHat_vpl = new real*[3];
    for (size_t k = 0; k < 3; k++) {
        partial_xiHat_vpl[k] = new real[3];
        for (size_t k2 = 0; k2 < 3; k2++) {
            partial_xiHat_vpl[k][k2] =
                (vPlanetCrossSHat * partial_vPlanetCrossSHatVec_vpl[k][k2] -
                    vPlanetCrossSHatVec[k2] * partial_vPlanetCrossSHat_vpl[k]) /
                vPlanetCrossSHat / vPlanetCrossSHat;
        }
        partial_zetaHat_vpl[k] = new real[3];
        vcross(sHat, partial_xiHat_vpl[k], partial_zetaHat_vpl[k]);
        for (size_t k2 = 0; k2 < 3; k2++) {
            partial_zetaHat_vpl[k][k2] *= -1.0;
        }
    }
    real *partial_xi_vpl = new real[3];
    real *partial_zeta_vpl = new real[3];
    for (size_t k = 0; k < 3; k++) {
        vdot(partial_zetaHat_vpl[k], hVec, 3, partial_xi_vpl[k]);
        partial_xi_vpl[k] /= vInf;
        vdot(partial_xiHat_vpl[k], hVec, 3, partial_zeta_vpl[k]);
        partial_zeta_vpl[k] /= -vInf;
    }
    std::vector<std::vector<real>> opikTotalDerivTerm2(2, std::vector<real>(6,0.0));
    for (size_t k = 0; k < 6; k++) {
        vdot(partial_xi_vpl, partial_vel_planet[k], 3, opikTotalDerivTerm2[0][k]);
        vdot(partial_zeta_vpl, partial_vel_planet[k], 3, opikTotalDerivTerm2[1][k]);
    }
    real **partial_xiHat = new real*[6];
    for (size_t k = 0; k < 6; k++) {
        partial_xiHat[k] = new real[3];
        vcross(vCentralBodyHelio, partial_sHat[k], temp2Vec);
        temp1 = xiHat[0] * temp2Vec[0] + xiHat[1] * temp2Vec[1] +
                xiHat[2] * temp2Vec[2];
        for (size_t k2 = 0; k2 < 3; k2++) {
            partial_xiHat[k][k2] = (temp2Vec[k2] - temp1*xiHat[k2])/vPlanetCrossSHat;
        }
    }
    real **partial_zetaHat = new real*[6];
    for (size_t k = 0; k < 6; k++) {
        partial_zetaHat[k] = new real[3];
        vcross(partial_sHat[k], xiHat, temp1Vec);
        vcross(sHat, partial_xiHat[k], temp2Vec);
        for (size_t k2 = 0; k2 < 3; k2++) {
            partial_zetaHat[k][k2] = -(temp1Vec[k2] + temp2Vec[k2]);
        }
    }
    std::vector<real> o_dx(6, 0.0);
    std::vector<real> o_dy(6, 0.0);
    for (size_t k = 0; k < 6; k++) {
        vdot(zetaHat, partial_hVec[k], 3, temp1);
        vdot(hVec, partial_zetaHat[k], 3, temp2);
        o_dx[k] = (temp1 + temp2 - ca->opik.x*partial_vInf[k])/vInf;
        vdot(xiHat, partial_hVec[k], 3, temp1);
        vdot(hVec, partial_xiHat[k], 3, temp2);
        o_dy[k] = -(temp1 + temp2 + ca->opik.y*partial_vInf[k])/vInf;
    }
    // map off-nominal trajectory for full partials
    vec_mat_mul(o_dx, partial_xCA, 6, ca->opik.dx);
    vec_mat_mul(o_dy, partial_xCA, 6, ca->opik.dy);
    for (size_t k = 0; k < 6; k++) {
        ca->opik.dx[k] += opikTotalDerivTerm2[0][k];
        ca->opik.dy[k] += opikTotalDerivTerm2[1][k];
    }
    // calculate mtp derivatives
    real posCA[3], velCA[3];
    posCA[0] = ca->xRelCA[0];
    posCA[1] = ca->xRelCA[1];
    posCA[2] = ca->xRelCA[2];
    velCA[0] = ca->xRelCA[3];
    velCA[1] = ca->xRelCA[4];
    velCA[2] = ca->xRelCA[5];
    real eHatX[3], eHatY[3], eHatZ[3], vHatCrosseHatZ[3];
    vunit(velCA, 3, eHatZ);
    vcross(vHat, eHatZ, vHatCrosseHatZ);
    vunit(vHatCrosseHatZ, 3, eHatY);
    vcross(eHatY, eHatZ, eHatX);
    real **partial_eHatZ = new real*[6];
    real vCA;
    vnorm(velCA, 3, vCA);
    for (size_t k = 0; k < 6; k++) {
        partial_eHatZ[k] = new real[3];
        for (size_t k2 = 0; k2 < 3; k2++) {
            temp1Vec[k2] = partial_xCA[k2+3][k];
            partial_eHatZ[k][k2] = 0;
        }
        vdot(eHatZ, temp1Vec, 3, temp1);
        for (size_t k2 = 0; k2 < 3; k2++) {
            partial_eHatZ[k][k2] = (temp1Vec[k2] - temp1*eHatZ[k2])/vCA;
        }
    }
    real **partial_eHatY = new real*[6];
    vnorm(vHatCrosseHatZ, 3, temp2);
    for (size_t k = 0; k < 6; k++) {
        partial_eHatY[k] = new real[3];
        for (size_t k2 = 0; k2 < 3; k2++) {
            vcross(vHat, partial_eHatZ[k], temp1Vec);
            vdot(eHatY, temp1Vec, 3, temp1);
            partial_eHatY[k][k2] = (temp1Vec[k2] - temp1*eHatY[k2])/temp2;
        }
    }
    real **partial_eHatX = new real*[6];
    for (size_t k = 0; k < 6; k++) {
        partial_eHatX[k] = new real[3];
        vcross(partial_eHatY[k], eHatZ, temp1Vec);
        vcross(eHatY, partial_eHatZ[k], temp2Vec);
        for (size_t k2 = 0; k2 < 3; k2++) {
            partial_eHatX[k][k2] = temp1Vec[k2] + temp2Vec[k2];
        }
    }
    std::vector<real> m_dx(6, 0.0);
    std::vector<real> m_dy(6, 0.0);
    for (size_t k = 0; k < 6; k++) {
        for (size_t k2 = 0; k2 < 3; k2++) {
            temp1Vec[k2] = partial_xCA[k2][k];
        }
        vdot(temp1Vec, eHatX, 3, temp1);
        vdot(posCA, partial_eHatX[k], 3, temp2);
        m_dx[k] = temp1 + temp2;
        vdot(temp1Vec, eHatY, 3, temp1);
        vdot(posCA, partial_eHatY[k], 3, temp2);
        m_dy[k] = temp1 + temp2;
    }
    // map off-nominal trajectory for full partials
    vec_mat_mul(m_dx, partial_xCA, 6, ca->mtp.dx);
    vec_mat_mul(m_dy, partial_xCA, 6, ca->mtp.dy);
    // clean up
    for (size_t k = 0; k < 6; k++) {
        delete[] partial_r_vec[k];
        delete[] partial_v_vec[k];
        delete[] partial_hVec[k];
        delete[] partial_eVec[k];
        delete[] partial_pHat[k];
        delete[] partial_qHat[k];
        delete[] partial_rHat[k];
        delete[] partial_sHat[k];
        delete[] partial_tHat[k];
        delete[] partial_vel_planet[k];
        delete[] partial_xiHat[k];
        delete[] partial_zetaHat[k];
        delete[] partial_eHatZ[k];
        delete[] partial_eHatY[k];
        delete[] partial_eHatX[k];
    }
    for (size_t k = 0; k < 3; k++) {
        delete[] partial_vpl_vpl[k];
        delete[] partial_vPlanetCrossSHatVec_vpl[k];
        delete[] partial_xiHat_vpl[k];
        delete[] partial_zetaHat_vpl[k];
    }
    delete[] partial_alpha;
    delete[] partial_vInf;
    delete[] partial_a;
    delete[] partial_n;
    delete[] partial_hVec;
    delete[] partial_h;
    delete[] partial_e;
    delete[] partial_F;
    delete[] partial_pHat;
    delete[] partial_qHat;
    delete[] partial_rHat;
    delete[] partial_sHat;
    delete[] partial_tHat;
    delete[] partial_lambda;
    delete[] partial_vel_planet;
    delete[] partial_vpl_vpl;
    delete[] partial_vPlanetCrossSHatVec_vpl;
    delete[] partial_vPlanetCrossSHat_vpl;
    delete[] partial_xiHat_vpl;
    delete[] partial_zetaHat_vpl;
    delete[] partial_xi_vpl;
    delete[] partial_zeta_vpl;
    delete[] partial_xiHat;
    delete[] partial_zetaHat;
    delete[] partial_eHatZ;
    delete[] partial_eHatY;
    delete[] partial_eHatX;
}

/**
 * @param[in] propSim PropSimulation object.
 * @param[in] tMap Time of the B-plane map.
 */
void CloseApproachParameters::get_ca_parameters(PropSimulation *propSim, const real &tMap) {
    // Calculate the parameters of a close approach
    const size_t i = this->flybyBodyIdx;
    const size_t j = this->centralBodyIdx;
    this->tMap = tMap;
    this->xRelMap = get_rel_state(propSim, i, j, tMap);
    this->flybyBody = propSim->integBodies[i].name;
    real mu;
    real radius;
    if (j < propSim->integParams.nInteg) {
        this->centralBody = propSim->integBodies[j].name;
        this->centralBodySpiceId = propSim->integBodies[j].spiceId;
        mu = propSim->integBodies[j].mass * propSim->consts.G;
        radius = propSim->integBodies[j].radius;
    } else {
        this->centralBody =
            propSim->spiceBodies[j - propSim->integParams.nInteg].name;
        this->centralBodySpiceId =
            propSim->spiceBodies[j - propSim->integParams.nInteg].spiceId;
        mu = propSim->spiceBodies[j - propSim->integParams.nInteg].mass *
            propSim->consts.G;
        radius = propSim->spiceBodies[j - propSim->integParams.nInteg].radius;
        radius += get_atm_offset(this->centralBodySpiceId);
    }
    // calculate B-plane parameters (Kizner B.R, B.T formulation)
    const real r = sqrt(xRelMap[0] * xRelMap[0] + xRelMap[1] * xRelMap[1] +
                        xRelMap[2] * xRelMap[2]);
    this->dist = r;
    const real v = sqrt(xRelMap[3] * xRelMap[3] + xRelMap[4] * xRelMap[4] +
                        xRelMap[5] * xRelMap[5]);
    this->vel = v;
    const real a = (mu * r) / (2 * mu - r * v * v);
    const real vInf = sqrt(-mu / a);
    this->vInf = vInf;
    real h, pos[3], vel[3], hVec[3];
    for (size_t k = 0; k < 3; k++) {
        pos[k] = xRelMap[k];
        vel[k] = xRelMap[k + 3];
    }
    vcross(pos, vel, hVec);
    vnorm(hVec, 3, h);
    real vCrossH[3], eVec[3], e;
    vcross(vel, hVec, vCrossH);
    for (size_t k = 0; k < 3; k++) {
        eVec[k] = vCrossH[k] / mu - pos[k] / r;
    }
    vnorm(eVec, 3, e);
    real pHat[3], qHat[3], hCrossPHat[3];
    vunit(eVec, 3, pHat);
    vcross(hVec, pHat, hCrossPHat);
    for (size_t k = 0; k < 3; k++) {
        qHat[k] = hCrossPHat[k] / h;
    }
    real sHat[3];
    for (size_t k = 0; k < 3; k++) {
        sHat[k] = pHat[k]/e + sqrt(e * e - 1) * qHat[k]/e;
    }
    real sHatCrossHVec[3], bVec[3];
    vcross(sHat, hVec, sHatCrossHVec);
    for (size_t k = 0; k < 3; k++) {
        bVec[k] = sHatCrossHVec[k] / vInf;
        this->bVec[k] = bVec[k];
    }
    vnorm(bVec, 3, this->bMag);
    real vHat[3] = {0.0, 0.0, -1.0};
    real rHat[3], tHat[3], vCrossSHatVec[3];
    vcross(vHat, sHat, vCrossSHatVec);
    vunit(vCrossSHatVec, 3, tHat);
    vcross(sHat, tHat, rHat);
    vdot(bVec, rHat, 3, this->kizner.x);
    vdot(bVec, tHat, 3, this->kizner.y);
    vdot(bVec, sHat, 3, this->kizner.z);
    this->gravFocusFactor = sqrt(1.0 + 2 * mu / radius / vInf / vInf);
    this->impact = this->bMag <= radius * this->gravFocusFactor;
    // if vInf is nan, compute distance at CA and check for impact
    if (std::isnan(vInf)) {
        const real distCA = sqrt(this->xRelCA[0] * this->xRelCA[0] +
                                 this->xRelCA[1] * this->xRelCA[1] +
                                 this->xRelCA[2] * this->xRelCA[2]);
        this->impact = distCA <= radius;
    }
    this->scaled.x = this->kizner.x / this->gravFocusFactor;
    this->scaled.y = this->kizner.y / this->gravFocusFactor;
    this->scaled.z = this->kizner.z / this->gravFocusFactor;
    real posDotVel;
    vdot(pos, vel, 3, posDotVel);
    const real sinhF = vInf*posDotVel/mu/e;
    const real F = asinh(sinhF); // or F = -log(2*r/a/e);
    const real n = sqrt(-mu / a / a / a);
    this->tPeri = tMap + (F - e*sinhF)/n;
    this->tLin = this->tPeri - log(e) / n;
    // calculate B-plane parameters (Öpik xi, zeta formulation)
    double xCentralBody[9];
    get_spk_state(this->centralBodySpiceId, tMap, propSim->spkEphem, xCentralBody);
    double xSun[9];
    get_spk_state(10, tMap, propSim->spkEphem, xSun);
    real vCentralBodyHelio[3];
    for (size_t k = 0; k < 3; k++) {
        vCentralBodyHelio[k] = xCentralBody[3+k]-xSun[3+k];
    }
    real xiHat[3], zetaHat[3], vPlanetCrossSHatVec[3];
    vcross(vCentralBodyHelio, sHat, vPlanetCrossSHatVec);
    vunit(vPlanetCrossSHatVec, 3, xiHat);
    vcross(sHat, xiHat, zetaHat);
    for (size_t k = 0; k < 3; k++) {
        zetaHat[k] *= -1;
    }
    vdot(bVec, xiHat, 3, this->opik.x);
    vdot(bVec, zetaHat, 3, this->opik.y);
    vdot(bVec, sHat, 3, this->opik.z);
    real posCA[3], velCA[3];
    posCA[0] = this->xRelCA[0];
    posCA[1] = this->xRelCA[1];
    posCA[2] = this->xRelCA[2];
    velCA[0] = this->xRelCA[3];
    velCA[1] = this->xRelCA[4];
    velCA[2] = this->xRelCA[5];
    real eHatX[3], eHatY[3], eHatZ[3], vHatCrosseHatZ[3];
    vunit(velCA, 3, eHatZ);
    vcross(vHat, eHatZ, vHatCrosseHatZ);
    vunit(vHatCrosseHatZ, 3, eHatY);
    vcross(eHatY, eHatZ, eHatX);
    vdot(posCA, eHatX, 3, this->mtp.x);
    vdot(posCA, eHatY, 3, this->mtp.y);
    vdot(posCA, eHatZ, 3, this->mtp.z);
    if (propSim->integBodies[i].propStm){
        get_bplane_partials(propSim, this, mu, radius);
    }
}

/**
 * @param[in] prec Precision of the output.
 */
void CloseApproachParameters::print_summary(int prec){
    std::cout.precision(prec);
    std::cout << "MJD " << this->t << " TDB:" << std::endl;
    std::cout << "    " << this->flybyBody << " approached " << this->centralBody << " at " << this->dist << " AU." << std::endl;
    std::cout << "    Relative Velocity: " << this->vel << " AU/d. V-infinity: " << this->vInf << " AU/d." << std::endl;
    std::cout << "    Gravitational focusing factor: " << this->gravFocusFactor << ". Impact: " << std::boolalpha << this->impact << std::endl;
}

/**
 * @param[in] propSim PropSimulation object.
 */
void ImpactParameters::get_impact_parameters(PropSimulation *propSim){
    std::string baseBodyFrame;
    get_baseBodyFrame(this->centralBodySpiceId, this->t, baseBodyFrame);
    std::vector<std::vector<real>> rotMat(6, std::vector<real>(6));
    get_pck_rotMat("J2000", baseBodyFrame, this->t, propSim->pckEphem, rotMat);
    mat_vec_mul(rotMat, this->xRel, this->xRelBodyFixed);
    real x, y, z, lon, lat, alt;
    x = this->xRelBodyFixed[0];
    y = this->xRelBodyFixed[1];
    z = this->xRelBodyFixed[2];
    if (this->centralBodySpiceId == 399){
        rec_to_geodetic(x, y, z, lon, lat, alt);
    } else {
        const real dist = sqrt(x*x + y*y + z*z);
        lat = asin(z/dist);
        lon = atan2(y, x);
        if (lon < 0.0) {
            lon += 2 * PI;
        }
        real centralBodyRadius;
        if ((size_t) this->centralBodyIdx < propSim->integParams.nInteg) {
            centralBodyRadius = propSim->integBodies[this->centralBodyIdx].radius;
        } else {
            centralBodyRadius = propSim->spiceBodies[this->centralBodyIdx - propSim->integParams.nInteg].radius;
        }
        alt = dist-centralBodyRadius;
    }
    this->lon = lon;
    this->lat = lat;
    this->alt = alt;
}

/**
 * @param[in] prec Precision of the output.
 */
void ImpactParameters::print_summary(int prec){
    std::cout.precision(prec);
    std::cout << "MJD " << this->t << " TDB:" << std::endl;
    std::cout << "    " << this->flybyBody << " impacted " << this->centralBody << " with a relative velocity of " << this->vel << " AU/d." << std::endl;
    std::cout << "    Impact location: " << std::endl;
    std::cout << "        Longitude: " << this->lon*RAD2DEG << " deg" << std::endl;
    std::cout << "        Latitude: " << this->lat*RAD2DEG << " deg" << std::endl;
    std::cout << "        Altitude: " << this->alt*1.495978707e8 << " km" << std::endl;
}
