#include "approach.h"

void check_ca_or_impact(propSimulation *propSim, const real &tOld,
                        const std::vector<real> xIntegOld, const real &t,
                        const std::vector<real> xInteg, int &keepStepping) {
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
                    relDist <= radius + radiusj) {
                    std::cout << "Impact detected at MJD " << t << " TDB. Body "
                              << propSim->integBodies[i].name
                              << " collided with body " << bodyj->name
                              << ". Terminating simulation!" << std::endl;
                    keepStepping = 0;
                    return;
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
                    real xRelCA[6];
                    get_ca_time(propSim, i, j, tOld, t, tCA);
                    get_ca_state(propSim, i, j, tCA, xRelCA);
                    CloseApproachParameters ca;
                    ca.tCA = tCA;
                    for (size_t k = 0; k < 6; k++) {
                        ca.xRelCA[k] = xRelCA[k];
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
    // This is used in root_brent to find the time of closest approach
    real xRel[6];
    get_ca_state(propSim, i, j, t, xRel);
    real relDist =
        sqrt(xRel[0] * xRel[0] + xRel[1] * xRel[1] + xRel[2] * xRel[2]);
    rDot =
        (xRel[0] * xRel[3] + xRel[1] * xRel[4] + xRel[2] * xRel[5]) / relDist;
}

void get_ca_state(propSimulation *propSim, const size_t &i, const size_t &j,
                  const real &t, real xRelCA[6]) {
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
            xRelCA[k] = xInterp[starti + k] - xInterp[startj + k];
        }
    } else {
        double xSpice[9];
        get_spk_state(
            propSim->spiceBodies[j - propSim->integParams.nInteg].spiceId, t,
            propSim->ephem, xSpice);
        for (size_t k = 0; k < 6; k++) {
            xRelCA[k] = xInterp[starti + k] - xSpice[k];
        }
    }
}

void get_ca_time(propSimulation *propSim, const size_t &i, const size_t &j,
                 const real &x1, const real &x2, real &tCA) {
    // Brent's method for root finding, from Numerical Recipes in C/C++, 3rd
    // edition, p. 454
    const real tol = 1.0e-3 / 86400.0;  // 1 millisecond
    const size_t maxIter = 100;
    const real eps = std::numeric_limits<real>::epsilon();
    real a = x1;
    real b = x2;
    real c = x2;
    real d = 0.0;
    real e = 0.0;
    real min1, min2;
    real fa, fb, fc, p, q, r, s, tol1, xm;
    ca_rdot_calc(propSim, i, j, a, fa);
    ca_rdot_calc(propSim, i, j, b, fb);
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
        throw std::runtime_error("Root must be bracketed in root_brent");
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
        ca_rdot_calc(propSim, i, j, b, fb);
    }
    throw std::runtime_error(
        "Maximum number of iterations exceeded in root_brent");
}

void CloseApproachParameters::get_ca_parameters(propSimulation *propSim, const real &tMap) {
    // Calculate the parameters of a close approach
    this->tMap = tMap;
    real xRelMap[6];
    const size_t i = this->flybyBodyIdx;
    const size_t j = this->centralBodyIdx;
    get_ca_state(propSim, i, j, tMap, xRelMap);
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
    const real F = -asinh(vInf * posDotVel / mu / e);
    const real n = sqrt(-mu / a / a / a);
    this->tPeri = tMap - (vInf * posDotVel / mu - F) / n;
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
        zetaHat[k] = -zetaHat[k];
    }
    vdot(bVec, xiHat, 3, this->opik.x);
    vdot(bVec, zetaHat, 3, this->opik.y);
    vdot(bVec, sHat, 3, this->opik.z);
    real eHatX[3], eHatY[3], eHatZ[3], vVecCrosseHatZ[3];
    vunit(vel, 3, eHatZ);
    vcross(vVec, eHatZ, vVecCrosseHatZ);
    vunit(vVecCrosseHatZ, 3, eHatY);
    vcross(eHatY, eHatZ, eHatX);
    vdot(pos, eHatX, 3, this->mtp.x);
    vdot(pos, eHatY, 3, this->mtp.y);
    vdot(pos, eHatZ, 3, this->mtp.z);
}
