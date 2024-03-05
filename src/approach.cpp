#include "approach.h"

void check_ca_or_impact(propSimulation *propSim, const real &tOld,
                        const std::vector<real> xIntegOld, const real &t,
                        const std::vector<real> xInteg) {
    // Check for close approach or impact
    // If a close approach is detected, the time of closest approach is
    // determined by root finding using Brent's method.
    // If an impact is detected, the simulation is terminated.
    const bool forwardProp = propSim->integParams.t0 < propSim->integParams.tf;
    const bool backwardProp = propSim->integParams.t0 > propSim->integParams.tf;
    size_t starti = 0;
    real posOld[3], pos[3], velOld[3], vel[3];
    real relPosOld[3], relPos[3], relVelOld[3], relVel[3];
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        const real radius = propSim->integBodies[i].radius;
        for (size_t j = 0; j < 3; j++) {
            posOld[j] = xIntegOld[starti + j];
            pos[j] = xInteg[starti + j];
            velOld[j] = xIntegOld[starti + 3 + j];
            vel[j] = xInteg[starti + 3 + j];
        }
        real radiusj;
        Body *bodyj;
        size_t startj = 0;
        for (size_t j = 0; j < propSim->integParams.nTotal; j++) {
            if (i != j) {
                real relDistOld, relDist;
                if (j < propSim->integParams.nInteg) {
                    bodyj = &propSim->integBodies[j];
                    radiusj = bodyj->radius;
                    for (size_t k = 0; k < 3; k++) {
                        relPosOld[k] = posOld[k] - xIntegOld[startj + k];
                        relPos[k] = pos[k] - xInteg[startj + k];
                        relVelOld[k] = velOld[k] - xIntegOld[startj + 3 + k];
                        relVel[k] = vel[k] - xInteg[startj + 3 + k];
                    }
                } else {
                    bodyj =
                        &propSim->spiceBodies[j - propSim->integParams.nInteg];
                    radiusj = bodyj->radius;
                    const real atm_offset = 100.0e3/propSim->consts.du2m;
                    radiusj += atm_offset;
                    double xSpiceOld[9];
                    get_spk_state(bodyj->spiceId, tOld, propSim->ephem,
                                  xSpiceOld);
                    for (size_t k = 0; k < 3; k++) {
                        relPosOld[k] = posOld[k] - xSpiceOld[k];
                        relPos[k] = pos[k] - bodyj->pos[k];
                        relVelOld[k] = velOld[k] - xSpiceOld[3 + k];
                        relVel[k] = vel[k] - bodyj->vel[k];
                    }
                }
                relDistOld = sqrt(relPosOld[0] * relPosOld[0] +
                                  relPosOld[1] * relPosOld[1] +
                                  relPosOld[2] * relPosOld[2]);
                relDist = sqrt(relPos[0] * relPos[0] + relPos[1] * relPos[1] +
                               relPos[2] * relPos[2]);
                // check impacts with spiceBodies first
                if (j > propSim->integParams.nInteg &&
                    relDist <= radius + radiusj && relDistOld > radius + radiusj) {
                    real tImp;
                    real xRelImp[6];
                    get_ca_or_impact_time(propSim, i, j, tOld, t, tImp, impact_r_calc);
                    get_rel_state(propSim, i, j, tImp, xRelImp);
                    ImpactParameters impact;
                    impact.t = tImp;
                    for (size_t k = 0; k < 6; k++) {
                        impact.xRel[k] = xRelImp[k];
                    }
                    impact.flybyBodyIdx = i;
                    impact.centralBodyIdx = j;
                    impact.get_ca_parameters(propSim, tImp);
                    impact.get_impact_parameters(propSim);
                    impact.impact = true;
                    propSim->impactParams.push_back(impact);
                    std::cout << "Impact detected at MJD " << tImp << " TDB. "
                              << propSim->integBodies[i].name
                              << " collided with " << bodyj->name << "!" << std::endl;
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
                const bool relDistWithinTol = relDist <= bodyj->caTol ||
                                              relDistOld <= bodyj->caTol;
                if (relDistMinimum && relDistWithinTol) {
                    real tCA;
                    real xRel[6];
                    get_ca_or_impact_time(propSim, i, j, tOld, t, tCA, ca_rdot_calc);
                    get_rel_state(propSim, i, j, tCA, xRel);
                    CloseApproachParameters ca;
                    ca.t = tCA;
                    for (size_t k = 0; k < 6; k++) {
                        ca.xRel[k] = xRel[k];
                    }
                    ca.flybyBodyIdx = i;
                    ca.centralBodyIdx = j;
                    ca.get_ca_parameters(propSim, tCA);
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

void ca_rdot_calc(propSimulation *propSim, const size_t &i, const size_t &j,
                  const real &t, real &rDot) {
    // Calculate the radial velocity between two bodies at a given time
    // This is used in get_ca_or_impact_time to find the time of closest approach
    real xRel[6];
    get_rel_state(propSim, i, j, t, xRel);
    real relDist =
        sqrt(xRel[0] * xRel[0] + xRel[1] * xRel[1] + xRel[2] * xRel[2]);
    rDot =
        (xRel[0] * xRel[3] + xRel[1] * xRel[4] + xRel[2] * xRel[5]) / relDist;
}

void impact_r_calc(propSimulation *propSim, const size_t &i, const size_t &j,
                  const real &t, real &r) {
    // Calculate the distance between two bodies at a given time, accounting for
    // the radius of the bodies
    real xRel[6];
    get_rel_state(propSim, i, j, t, xRel);
    real relDist =
        sqrt(xRel[0] * xRel[0] + xRel[1] * xRel[1] + xRel[2] * xRel[2]);
    real radius, radiusj;
    radius = propSim->integBodies[i].radius;
    if (j < propSim->integParams.nInteg) {
        radiusj = propSim->integBodies[j].radius;
    } else {
        radiusj = propSim->spiceBodies[j - propSim->integParams.nInteg].radius;
        const real atm_offset = 100.0e3/propSim->consts.du2m;
        radiusj += atm_offset;
    }
    r = relDist - radius - radiusj;
}

void get_rel_state(propSimulation *propSim, const size_t &i, const size_t &j,
                   const real &t, real xRel[6]) {
    std::vector<real> xInterp = propSim->interpolate(t);
    size_t starti = 0;
    for (size_t k = 0; k < i; k++) {
        starti += 2*propSim->integBodies[k].n2Derivs;
    }
    if (j < propSim->integParams.nInteg) {
        size_t startj = 0;
        for (size_t k = 0; k < j; k++) {
            startj += 2*propSim->integBodies[k].n2Derivs;
        }
        for (size_t k = 0; k < 6; k++) {
            xRel[k] = xInterp[starti + k] - xInterp[startj + k];
        }
    } else {
        double xSpice[9];
        get_spk_state(
            propSim->spiceBodies[j - propSim->integParams.nInteg].spiceId, t,
            propSim->ephem, xSpice);
        for (size_t k = 0; k < 6; k++) {
            xRel[k] = xInterp[starti + k] - xSpice[k];
        }
    }
}

void get_ca_or_impact_time(propSimulation *propSim, const size_t &i,
                           const size_t &j, const real &x1, const real &x2,
                           real &tCA,
                           void (*zero_func)(propSimulation *, const size_t &,
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

void CloseApproachParameters::get_ca_parameters(propSimulation *propSim, const real &tMap) {
    // Calculate the parameters of a close approach
    this->tMap = tMap;
    real xRelMap[6];
    const size_t i = this->flybyBodyIdx;
    const size_t j = this->centralBodyIdx;
    get_rel_state(propSim, i, j, tMap, xRelMap);
    for (size_t k = 0; k < 6; k++) {
        this->xRelMap[k] = xRelMap[k];
    }
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
    get_spk_state(this->centralBodySpiceId, tMap, propSim->ephem, xCentralBody);
    double xSun[9];
    get_spk_state(10, tMap, propSim->ephem, xSun);
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
    }
}

void CloseApproachParameters::print_summary(int prec){
    std::cout.precision(prec);
    std::cout << "MJD " << this->t << " TDB:" << std::endl;
    std::cout << "    " << this->flybyBody << " approached " << this->centralBody << " at " << this->dist << " AU." << std::endl;
    std::cout << "    Relative Velocity: " << this->vel << " AU/d. V-infinity: " << this->vInf << " AU/d." << std::endl;
    std::cout << "    Gravitational focusing factor: " << this->gravFocusFactor << ". Impact: " << std::boolalpha << this->impact << std::endl;
}

void ImpactParameters::get_impact_parameters(propSimulation *propSim){
    ConstSpiceChar *baseBodyFrame;
    switch (this->centralBodySpiceId) {
        case 10:
            baseBodyFrame = "IAU_SUN";
            break;
        case 1:
        case 199:
            baseBodyFrame = "IAU_MERCURY";
            break;
        case 2:
        case 299:
            baseBodyFrame = "IAU_VENUS";
            break;
        case 399:
            baseBodyFrame = "ITRF93";
            // High precision frame is not defined before 1972 JAN 01 00:00:42.183 TDB
            if (this->t < 41317.0004882291666666666L) {
                baseBodyFrame = "IAU_EARTH";
            }
            break;
        case 499:
            baseBodyFrame = "IAU_MARS";
            break;
        case 599:
            baseBodyFrame = "IAU_JUPITER";
            break;
        case 699:
            baseBodyFrame = "IAU_SATURN";
            break;
        case 799:
            baseBodyFrame = "IAU_URANUS";
            break;
        case 899:
            baseBodyFrame = "IAU_NEPTUNE";
            break;
        case 999:
            baseBodyFrame = "IAU_PLUTO";
            break;
        default:
            std::cout << "get_impact_parameters: Given impacted body: "
                      << this->centralBody << std::endl;
            throw std::invalid_argument("Given base body not supported");
            break;
    }
    real tET;
    mjd_to_et(this->t, tET);
    SpiceDouble rotMat[6][6];
    sxform_c("J2000", baseBodyFrame, tET, rotMat);
    SpiceDouble impactRelStateInertial[6];
    impactRelStateInertial[0] = this->xRel[0]*propSim->consts.du2m/1.0e3L;
    impactRelStateInertial[1] = this->xRel[1]*propSim->consts.du2m/1.0e3L;
    impactRelStateInertial[2] = this->xRel[2]*propSim->consts.du2m/1.0e3L;
    impactRelStateInertial[3] = this->xRel[3]*propSim->consts.duptu2mps/1.0e3L;
    impactRelStateInertial[4] = this->xRel[4]*propSim->consts.duptu2mps/1.0e3L;
    impactRelStateInertial[5] = this->xRel[5]*propSim->consts.duptu2mps/1.0e3L;
    SpiceDouble impactRelStateBodyFixed[6];
    mxvg_c(rotMat, impactRelStateInertial, 6, 6, impactRelStateBodyFixed);
    impactRelStateBodyFixed[0] *= 1.0e3L/propSim->consts.du2m;
    impactRelStateBodyFixed[1] *= 1.0e3L/propSim->consts.du2m;
    impactRelStateBodyFixed[2] *= 1.0e3L/propSim->consts.du2m;
    impactRelStateBodyFixed[3] *= 1.0e3L/propSim->consts.duptu2mps;
    impactRelStateBodyFixed[4] *= 1.0e3L/propSim->consts.duptu2mps;
    impactRelStateBodyFixed[5] *= 1.0e3L/propSim->consts.duptu2mps;
    for (size_t i = 0; i < 6; i++) {
        this->xRelBodyFixed[i] = (real)impactRelStateBodyFixed[i];
    }
    SpiceDouble bodyFixedPos[3];
    for (size_t i = 0; i < 3; i++) {
        bodyFixedPos[i] = impactRelStateBodyFixed[i];
    }
    SpiceDouble lon, lat, alt;
    reclat_c(bodyFixedPos, &alt, &lon, &lat);
    if (lon < 0.0) {
        lon += 2 * PI;
    }
    real radiusj;
    if ((size_t) this->centralBodyIdx < propSim->integParams.nInteg) {
        radiusj = propSim->integBodies[this->centralBodyIdx].radius;
    } else {
        radiusj = propSim->spiceBodies[this->centralBodyIdx - propSim->integParams.nInteg].radius;
    }
    alt -= radiusj;
    this->lon = lon;
    this->lat = lat;
    this->alt = alt*propSim->consts.du2m/1.0e3L;
}

void ImpactParameters::print_summary(int prec){
    std::cout.precision(prec);
    std::cout << "MJD " << this->t << " TDB:" << std::endl;
    std::cout << "    " << this->flybyBody << " impacted " << this->centralBody << " with a relative velocity of " << this->vel << " AU/d." << std::endl;
    std::cout << "    Impact location: " << std::endl;
    std::cout << "        Longitude: " << this->lon*180.0L/PI << " deg" << std::endl;
    std::cout << "        Latitude: " << this->lat*180.0L/PI << " deg" << std::endl;
    std::cout << "        Altitude: " << this->alt << " km" << std::endl;
}
