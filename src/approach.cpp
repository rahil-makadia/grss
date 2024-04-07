#include "approach.h"

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
    real relPosOld[3], relPos[3], relVelOld[3], relVel[3];
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        const real flybyBodyRadius = propSim->integBodies[i].radius;
        real centralBodyRadius, caTol;
        size_t startj = 0;
        for (size_t j = 0; j < propSim->integParams.nTotal; j++) {
            if (i != j) {
                real relDistOld, relDist;
                if (j < propSim->integParams.nInteg) {
                    centralBodyRadius = propSim->integBodies[j].radius;
                    caTol = propSim->integBodies[j].caTol;
                    for (size_t k = 0; k < 3; k++) {
                        relPosOld[k] = xIntegOld[starti + k] - xIntegOld[startj + k];
                        relPos[k] = xInteg[starti + k] - xInteg[startj + k];
                        relVelOld[k] = xIntegOld[starti + 3 + k] - xIntegOld[startj + 3 + k];
                        relVel[k] = xInteg[starti + 3 + k] - xInteg[startj + 3 + k];
                    }
                } else {
                    SpiceBody bodyj = propSim->spiceBodies[j - propSim->integParams.nInteg];
                    centralBodyRadius = bodyj.radius;
                    caTol = bodyj.caTol;
                    const real atm_offset = 100.0e3/propSim->consts.du2m;
                    centralBodyRadius += atm_offset;
                    double xSpice[9], xSpiceOld[9];
                    get_spk_state(bodyj.spiceId, t, propSim->spkEphem, xSpice);
                    get_spk_state(bodyj.spiceId, tOld, propSim->spkEphem, xSpiceOld);
                    for (size_t k = 0; k < 3; k++) {
                        relPosOld[k] = xIntegOld[starti + k] - xSpiceOld[k];
                        relPos[k] = xInteg[starti + k] - xSpice[k];
                        relVelOld[k] = xIntegOld[starti + 3 + k] - xSpiceOld[3 + k];
                        relVel[k] = xInteg[starti + 3 + k] - xSpice[3 + k];
                    }
                }
                relDistOld = sqrt(relPosOld[0] * relPosOld[0] +
                                  relPosOld[1] * relPosOld[1] +
                                  relPosOld[2] * relPosOld[2]);
                relDist = sqrt(relPos[0] * relPos[0] + relPos[1] * relPos[1] +
                               relPos[2] * relPos[2]);
                // check impacts with spiceBodies first
                if (j > propSim->integParams.nInteg &&
                    relDist <= flybyBodyRadius + centralBodyRadius &&
                    relDistOld > flybyBodyRadius + centralBodyRadius) {
                    ImpactParameters impact;
                    get_ca_or_impact_time(propSim, i, j, tOld, t, impact.t, impact_r_calc);
                    impact.xRel = get_rel_state(propSim, i, j, impact.t);
                    impact.flybyBodyIdx = i;
                    impact.centralBodyIdx = j;
                    impact.get_ca_parameters(propSim, impact.t);
                    impact.get_impact_parameters(propSim);
                    impact.impact = true;
                    propSim->impactParams.push_back(impact);
                    // std::cout << "Impact detected at MJD " << impact.t
                    //             << " TDB. " << impact.flybyBody
                    //             << " collided with " << impact.centralBody
                    //             << "!" << std::endl;
                }
                // check close approach
                real radialVel, radialVelOld;
                radialVelOld =
                    (relPosOld[0] * relVelOld[0] + relPosOld[1] * relVelOld[1] +
                     relPosOld[2] * relVelOld[2]) /
                    relDistOld;
                radialVel = (relPos[0] * relVel[0] + relPos[1] * relVel[1] +
                             relPos[2] * relVel[2]) /
                    relDist;
                const bool relDistMinimum = (forwardProp && radialVelOld < 0.0 && radialVel >= 0.0) ||
                                            (backwardProp && radialVelOld > 0.0 && radialVel <= 0.0);
                if (j < propSim->integParams.nInteg) {
                    caTol = propSim->integBodies[j].caTol;
                } else {
                    caTol = propSim->spiceBodies[j - propSim->integParams.nInteg].caTol;
                }
                const bool relDistWithinTol = relDist <= caTol || relDistOld <= caTol;
                if (relDistMinimum && relDistWithinTol) {
                    CloseApproachParameters ca;
                    get_ca_or_impact_time(propSim, i, j, tOld, t, ca.t, ca_rdot_calc);
                    ca.xRel = get_rel_state(propSim, i, j, ca.t);
                    ca.flybyBodyIdx = i;
                    ca.centralBodyIdx = j;
                    ca.get_ca_parameters(propSim, ca.t);
                    propSim->caParams.push_back(ca);
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
    real relDist =
        sqrt(xRel[0] * xRel[0] + xRel[1] * xRel[1] + xRel[2] * xRel[2]);
    real flybyBodyRadius, centralBodyRadius;
    flybyBodyRadius = propSim->integBodies[i].radius;
    if (j < propSim->integParams.nInteg) {
        centralBodyRadius = propSim->integBodies[j].radius;
    } else {
        centralBodyRadius = propSim->spiceBodies[j - propSim->integParams.nInteg].radius;
        const real atm_offset = 100.0e3/propSim->consts.du2m;
        centralBodyRadius += atm_offset;
    }
    r = relDist - flybyBodyRadius - centralBodyRadius;
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] i Index of the first body.
 * @param[in] j Index of the second body.
 * @param[in] t Time to compute the relative state.
 * @return xRel Relative state of the body.
 */
static std::vector<real> get_rel_state(PropSimulation *propSim, const size_t &i,
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
        throw std::runtime_error("Root must be bracketed in get_ca_or_impact_time");
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
                 "get_ca_or_impact_time!!! Impact/CA time may not be accurate.";
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
    const real fac1 = 1 / e;
    const real fac2 = sqrt(e * e - 1);
    real sHat[3];
    for (size_t k = 0; k < 3; k++) {
        sHat[k] = fac1 * pHat[k] + fac1 * fac2 * qHat[k];
    }
    real sHatCrossHVec[3], bVec[3], b;
    vcross(sHat, hVec, sHatCrossHVec);
    for (size_t k = 0; k < 3; k++) {
        bVec[k] = sHatCrossHVec[k] / vInf;
        this->bVec[k] = bVec[k];
    }
    vnorm(bVec, 3, b);
    this->bMag = b;
    real vVec[3];
    vVec[0] = 0;
    vVec[1] = 0;
    vVec[2] = -1;
    real rHat[3], tHat[3], vCrossSHatVec[3];
    vcross(vVec, sHat, vCrossSHatVec);
    vunit(vCrossSHatVec, 3, tHat);
    vcross(sHat, tHat, rHat);
    vdot(bVec, rHat, 3, this->kizner.x);
    vdot(bVec, tHat, 3, this->kizner.y);
    vdot(bVec, sHat, 3, this->kizner.z);
    this->gravFocusFactor = sqrt(1.0 + 2 * mu / radius / vInf / vInf);
    this->impact = b <= radius * this->gravFocusFactor;
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
    real xiHat[3], zetaHat[3], vCentralBodyHelioCrossSHat[3];
    vcross(vCentralBodyHelio, sHat, vCentralBodyHelioCrossSHat);
    vunit(vCentralBodyHelioCrossSHat, 3, xiHat);
    vcross(sHat, xiHat, zetaHat);
    for (size_t k = 0; k < 3; k++) {
        zetaHat[k] *= -1;
    }
    vdot(bVec, xiHat, 3, this->opik.x);
    vdot(bVec, zetaHat, 3, this->opik.y);
    vdot(bVec, sHat, 3, this->opik.z);
    real posCA[3], velCA[3];
    posCA[0] = this->xRel[0];
    posCA[1] = this->xRel[1];
    posCA[2] = this->xRel[2];
    velCA[0] = this->xRel[3];
    velCA[1] = this->xRel[4];
    velCA[2] = this->xRel[5];
    real eHatX[3], eHatY[3], eHatZ[3], vVecCrosseHatZ[3];
    vunit(velCA, 3, eHatZ);
    vcross(vVec, eHatZ, vVecCrosseHatZ);
    vunit(vVecCrosseHatZ, 3, eHatY);
    vcross(eHatY, eHatZ, eHatX);
    vdot(posCA, eHatX, 3, this->mtp.x);
    vdot(posCA, eHatY, 3, this->mtp.y);
    vdot(posCA, eHatZ, 3, this->mtp.z);

    if (tMap == this->t && propSim->integBodies[i].propStm == true){
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

        real alpha = vInf*vInf;
        real *partial_alpha = new real[6];
        real *partial_vInf = new real[6];
        real *partial_a = new real[6];
        real *partial_n = new real[6];
        real temp1, temp2, rCA, rCA3, vCA;
        vnorm(posCA, 3, rCA);
        vnorm(velCA, 3, vCA);
        rCA3 = rCA*rCA*rCA;
        for (size_t k = 0; k < 6; k++) {
            vdot(velCA, partial_v_vec[k], 3, temp1);
            vdot(posCA, partial_r_vec[k], 3, temp2);
            partial_alpha[k] = 2*(temp1 + mu*temp2/rCA3);
            partial_vInf[k] = partial_alpha[k]/(2*vInf);
            partial_a[k] = -a*partial_alpha[k]/alpha;
            partial_n[k] = -3*n*partial_a[k]/2/a;
        }

        real **partial_hVec = new real*[6];
        real *partial_h = new real[6];
        for (size_t k = 0; k < 6; k++) {
            partial_hVec[k] = new real[3];
        }
        real temp1Vec[3];
        for (size_t k = 0; k < 6; k++) {
            vcross(partial_r_vec[k], velCA, partial_hVec[k]);
            vcross(posCA, partial_v_vec[k], temp1Vec);
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
        }
        for (size_t k = 0; k < 6; k++) {
            vcross(partial_v_vec[k], hVec, temp1Vec);
            for (size_t k2 = 0; k2 < 3; k2++) {
                partial_eVec[k][k2] = temp1Vec[k2] / mu;
            }
            vcross(velCA, partial_hVec[k], temp1Vec);
            for (size_t k2 = 0; k2 < 3; k2++) {
                partial_eVec[k][k2] += temp1Vec[k2] / mu;
            }
            for (size_t k2 = 0; k2 < 3; k2++) {
                partial_eVec[k][k2] -= partial_r_vec[k][k2] / r;
            }
            vdot(posCA, partial_r_vec[k], 3, temp1Vec[0]);
            for (size_t k2 = 0; k2 < 3; k2++) {
                partial_eVec[k][k2] += temp1Vec[0] * posCA[k2] / (r * r * r);
            }
            vdot(partial_eVec[k], eVec, 3, partial_e[k]);
            partial_e[k] /= e;
        }

        real *partial_F = new real[6];
        for (size_t k = 0; k < 6; k++) {
            vdot(posCA, partial_v_vec[k], 3, temp1);
            vdot(velCA, partial_r_vec[k], 3, temp2);
            partial_F[k] = (posDotVel*partial_vInf[k]/mu + vInf*(temp1+temp2)/mu - partial_e[k]*sinhF)/e/cosh(F);
        }
        std::vector<real> xCA = propSim->interpolate(tMap);
        std::vector<real> accCA = get_state_der(tMap, xCA, propSim);
        real accRel[3], accPlanet[3];
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
            accRel[k] = propSim->integBodies[i].acc[k] - accPlanet[k];
        }
        real posDotAcc, posCADotVelCA;
        vdot(posCA, accRel, 3, posDotAcc);
        vdot(posCA, velCA, 3, posCADotVelCA);

        for (size_t k = 0; k < 6; k++) {
            this->dTLinMinusT[k] = (r*partial_F[k]/a - (e*sinhF + 1)*partial_e[k]/e - (tLin-tMap)*partial_n[k])/n;
            vdot(posCA, partial_v_vec[k], 3, temp1);
            vdot(velCA, partial_r_vec[k], 3, temp2);
            this->dt[k] = -(temp1+temp2)/(posDotAcc + vCA*vCA);
        }

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
                partial_sHat[k][k2] = -partial_e[k]*pHat[k2]/e/e + partial_pHat[k][k2]/e + partial_e[k]*qHat[k2]/e/e/sqrt(e*e-1) + sqrt(e*e-1)*partial_qHat[k][k2]/e;
            }
            vcross(vVec, partial_sHat[k], temp1Vec);
            vdot(tHat, temp1Vec, 3, temp1);
            vcross(vVec, sHat, temp2Vec);
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
        real *partial_lambda = new real[6];
        for (size_t k = 0; k < 6; k++) {
            vdot(tHat, partial_hVec[k], 3, temp1);
            vdot(hVec, partial_tHat[k], 3, temp2);
            this->kizner.dx[k] = (temp1 + temp2 - this->kizner.x*partial_vInf[k])/vInf;
            vdot(rHat, partial_hVec[k], 3, temp1);
            vdot(hVec, partial_rHat[k], 3, temp2);

            partial_lambda[k] = -2*mu*partial_vInf[k]/this->gravFocusFactor/radius/vInf/vInf/vInf;
            this->kizner.dy[k] = -(temp1 + temp2 + this->kizner.y*partial_vInf[k])/vInf;
            this->scaled.dx[k] = (this->kizner.dx[k] - this->scaled.x*partial_lambda[k])/this->gravFocusFactor;
            this->scaled.dy[k] = (this->kizner.dy[k] - this->scaled.y*partial_lambda[k])/this->gravFocusFactor;
        }

        // // Acceleration of the Sun is not included in the calculation
        // // because it is negligibly small
        // real accSun[3];
        // memset(accSun, 0, 3*sizeof(real));
        // for (size_t k = 0; k < propSim->integParams.nSpice; k++) {
        //     if (propSim->spiceBodies[k].spiceId == 10) {
        //         accSun[0] = propSim->spiceBodies[k].acc[0];
        //         accSun[1] = propSim->spiceBodies[k].acc[1];
        //         accSun[2] = propSim->spiceBodies[k].acc[2];
        //         break;
        //     }
        // }
        real **partial_vel_planet = new real*[6];
        for (size_t k = 0; k < 6; k++) {
            partial_vel_planet[k] = new real[3];
            for (size_t k2 = 0; k2 < 3; k2++) {
                // partial_vel_planet[k][k2] = this->dt[k]*(accPlanet[k2] - accSun[k2]);
                partial_vel_planet[k][k2] = this->dt[k]*accPlanet[k2];
            }
        }
        real **partial_xiHat = new real*[6];
        real temp3, temp3Vec[3];
        real temp4, temp4Vec[3];
        for (size_t k = 0; k < 6; k++) {
            partial_xiHat[k] = new real[3];
            vcross(partial_v_vec[k], sHat, temp1Vec);
            vcross(vCentralBodyHelio, partial_sHat[k], temp2Vec);
            vcross(vCentralBodyHelio, sHat, temp3Vec);
            vnorm(temp3Vec, 3, temp3);
            for (size_t k2 = 0; k2 < 3; k2++) {
                temp4Vec[0] = (temp1Vec[k2]+temp2Vec[k2])*xiHat[0];
                temp4Vec[1] = (temp1Vec[k2]+temp2Vec[k2])*xiHat[1];
                temp4Vec[2] = (temp1Vec[k2]+temp2Vec[k2])*xiHat[2];
                vdot(xiHat, temp4Vec, 3, temp4);
                partial_xiHat[k][k2] = (temp1Vec[k2] + temp2Vec[k2] - temp4)/temp3;
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
        for (size_t k = 0; k < 6; k++) {
            vdot(zetaHat, partial_hVec[k], 3, temp1);
            vdot(hVec, partial_zetaHat[k], 3, temp2);
            this->opik.dx[k] = (temp1 + temp2 - this->opik.x*partial_vInf[k])/vInf;
            vdot(xiHat, partial_hVec[k], 3, temp1);
            vdot(hVec, partial_xiHat[k], 3, temp2);
            this->opik.dy[k] = -(temp1 + temp2 + this->opik.y*partial_vInf[k])/vInf;
        }

        this->mtp.dx = {eHatX[0], eHatX[1], eHatX[2], 0.0, 0.0, 0.0};
        this->mtp.dy = {eHatY[0], eHatY[1], eHatY[2], 0.0, 0.0, 0.0};
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
        delete[] partial_xiHat;
        delete[] partial_zetaHat;
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
    std::vector<real> impactRelStateInertial(6);
    impactRelStateInertial[0] = this->xRel[0]*propSim->consts.du2m/1.0e3L;
    impactRelStateInertial[1] = this->xRel[1]*propSim->consts.du2m/1.0e3L;
    impactRelStateInertial[2] = this->xRel[2]*propSim->consts.du2m/1.0e3L;
    impactRelStateInertial[3] = this->xRel[3]*propSim->consts.duptu2mps/1.0e3L;
    impactRelStateInertial[4] = this->xRel[4]*propSim->consts.duptu2mps/1.0e3L;
    impactRelStateInertial[5] = this->xRel[5]*propSim->consts.duptu2mps/1.0e3L;
    mat_vec_mul(rotMat, impactRelStateInertial, this->xRelBodyFixed);
    this->xRelBodyFixed[0] *= 1.0e3L/propSim->consts.du2m;
    this->xRelBodyFixed[1] *= 1.0e3L/propSim->consts.du2m;
    this->xRelBodyFixed[2] *= 1.0e3L/propSim->consts.du2m;
    this->xRelBodyFixed[3] *= 1.0e3L/propSim->consts.duptu2mps;
    this->xRelBodyFixed[4] *= 1.0e3L/propSim->consts.duptu2mps;
    this->xRelBodyFixed[5] *= 1.0e3L/propSim->consts.duptu2mps;
    real x, y, z, lon, lat, dist;
    x = this->xRelBodyFixed[0];
    y = this->xRelBodyFixed[1];
    z = this->xRelBodyFixed[2];
    dist = sqrt(x*x + y*y + z*z);
    lat = atan2(z, sqrt(x*x + y*y));
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
    this->lon = lon;
    this->lat = lat;
    this->alt = (dist-centralBodyRadius)*propSim->consts.du2m/1.0e3L;
}

/**
 * @param[in] prec Precision of the output.
 */
void ImpactParameters::print_summary(int prec){
    std::cout.precision(prec);
    std::cout << "MJD " << this->t << " TDB:" << std::endl;
    std::cout << "    " << this->flybyBody << " impacted " << this->centralBody << " with a relative velocity of " << this->vel << " AU/d." << std::endl;
    std::cout << "    Impact location: " << std::endl;
    std::cout << "        Longitude: " << this->lon*180.0L/PI << " deg" << std::endl;
    std::cout << "        Latitude: " << this->lat*180.0L/PI << " deg" << std::endl;
    std::cout << "        Altitude: " << this->alt << " km" << std::endl;
}
