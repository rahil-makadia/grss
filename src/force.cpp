#include "force.h"
// #define PRINT_FORCES 1
#ifdef PRINT_FORCES
#include <iomanip>
#include <fstream>
#endif

std::vector<real> get_state_der(const real &t, const std::vector<real> &xInteg,
                                propSimulation *propSim) {
    std::vector<real> accInteg(propSim->integParams.n2Derivs, 0.0);
    size_t starti = 0;
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        propSim->integBodies[i].pos[0] = xInteg[starti];
        propSim->integBodies[i].pos[1] = xInteg[starti + 1];
        propSim->integBodies[i].pos[2] = xInteg[starti + 2];
        propSim->integBodies[i].vel[0] = xInteg[starti + 3];
        propSim->integBodies[i].vel[1] = xInteg[starti + 4];
        propSim->integBodies[i].vel[2] = xInteg[starti + 5];
        if (propSim->integBodies[i].propStm) {
            for (size_t j = 0; j < propSim->integBodies[i].stm.size(); j++) {
                propSim->integBodies[i].stm[j] = xInteg[starti + 6 + j];
            }
        }
        starti += 2*propSim->integBodies[i].n2Derivs;
    }
    double xSpice[9];
    for (size_t i = 0; i < propSim->integParams.nSpice; i++) {
        get_spk_state(propSim->spiceBodies[i].spiceId, t, propSim->ephem,
                      xSpice);
        propSim->spiceBodies[i].pos[0] = xSpice[0];
        propSim->spiceBodies[i].pos[1] = xSpice[1];
        propSim->spiceBodies[i].pos[2] = xSpice[2];
        propSim->spiceBodies[i].vel[0] = xSpice[3];
        propSim->spiceBodies[i].vel[1] = xSpice[4];
        propSim->spiceBodies[i].vel[2] = xSpice[5];
        propSim->spiceBodies[i].acc[0] = xSpice[6];
        propSim->spiceBodies[i].acc[1] = xSpice[7];
        propSim->spiceBodies[i].acc[2] = xSpice[8];
    }
    #ifdef PRINT_FORCES
    std::ofstream forceFile;
    forceFile.precision(8);
    forceFile.setf(std::ios::scientific);
    forceFile.setf(std::ios::right, std::ios::adjustfield);
    forceFile.open("cpp.11", std::ios::app);
    forceFile << "time (MJD): " << std::setw(16) << t << " state:";
    forceFile.precision(16);
    for (size_t i = 0; i < 6; i++) {
        forceFile << std::setw(25) << xInteg[i];
    }
    forceFile << std::endl;
    forceFile.close();
    #endif
    force_newton(propSim, accInteg);
    // force_ppn_simple(propSim, accInteg);
    force_ppn_eih(propSim, accInteg);
    force_J2(propSim, accInteg);
    force_nongrav(propSim, accInteg);
    force_thruster(propSim, accInteg);
    #ifdef PRINT_FORCES
    forceFile.open("cpp.11", std::ios::app);
    forceFile << std::setw(10) << "total acc:" << std::setw(25) << accInteg[0]
              << std::setw(25) << accInteg[1] << std::setw(25) << accInteg[2]
              << std::endl;
    forceFile.close();
    #endif
    starti = 0;
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        propSim->integBodies[i].acc[0] = accInteg[starti];
        propSim->integBodies[i].acc[1] = accInteg[starti+1];
        propSim->integBodies[i].acc[2] = accInteg[starti+2];
        starti += propSim->integBodies[i].n2Derivs;
    }
    return accInteg;
}

void force_newton(const propSimulation *propSim, std::vector<real> &accInteg) {
    #ifdef PRINT_FORCES
    std::ofstream forceFile;
    forceFile.precision(16);
    forceFile.setf(std::ios::scientific);
    forceFile.setf(std::ios::right, std::ios::adjustfield);
    forceFile.open("cpp.11", std::ios::app);
    #endif
    const real G = propSim->consts.G;
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        const real x = propSim->integBodies[i].pos[0];
        const real y = propSim->integBodies[i].pos[1];
        const real z = propSim->integBodies[i].pos[2];
        for (size_t j = 0; j < propSim->integParams.nTotal; j++) {
            const Body *bodyj;
            if (j < propSim->integParams.nInteg) {
                bodyj = &propSim->integBodies[j];
            } else {
                bodyj = &propSim->spiceBodies[j - propSim->integParams.nInteg];
            }
            const real massj = bodyj->mass;
            if (i != j && massj != 0.0) {
                const real dx = x - bodyj->pos[0];
                const real dy = y - bodyj->pos[1];
                const real dz = z - bodyj->pos[2];
                const real rRel = sqrt(dx * dx + dy * dy + dz * dz);
                const real fac = -G * massj / (rRel * rRel * rRel);
                accInteg[3 * i + 0] += fac * dx;
                accInteg[3 * i + 1] += fac * dy;
                accInteg[3 * i + 2] += fac * dz;
                if (propSim->integBodies[i].propStm) {
                    stm_newton(propSim->integBodies[i], G*massj, dx, dy, dz, 3*i+3,
                               accInteg);
                }
                #ifdef PRINT_FORCES
                forceFile << std::setw(10) << bodyj->spiceId << std::setw(25)
                          << G * massj << std::setw(25) << dx << std::setw(25)
                          << dy << std::setw(25) << dz << std::setw(25)
                          << fac * dx << std::setw(25)
                          << fac * dy << std::setw(25)
                          << fac * dz << std::endl;
                #endif
            }
        }
    }
    #ifdef PRINT_FORCES
    forceFile.close();
    #endif
}

void force_ppn_simple(const propSimulation *propSim,
                      std::vector<real> &accInteg) {
#ifdef PRINT_FORCES
    std::ofstream forceFile;
    forceFile.precision(16);
    forceFile.setf(std::ios::scientific);
    forceFile.setf(std::ios::right, std::ios::adjustfield);
    forceFile.open("cpp.11", std::ios::app);
#endif
    const real G = propSim->consts.G;
    const real c = propSim->consts.clight;
    const real c2 = c * c;
    const real beta = 1.0L;
    const real gamma = 1.0L;
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        const real x = propSim->integBodies[i].pos[0];
        const real y = propSim->integBodies[i].pos[1];
        const real z = propSim->integBodies[i].pos[2];
        const real vx = propSim->integBodies[i].vel[0];
        const real vy = propSim->integBodies[i].vel[1];
        const real vz = propSim->integBodies[i].vel[2];
        for (size_t j = 0; j < propSim->integParams.nTotal; j++) {
            const Body *bodyj;
            if (j < propSim->integParams.nInteg) {
                bodyj = &propSim->integBodies[j];
            } else {
                bodyj = &propSim->spiceBodies[j - propSim->integParams.nInteg];
            }
            const real massj = bodyj->mass;
            if (i != j && massj != 0.0 && bodyj->spiceId == 10) {
                const real gm = G * massj;
                const real dx = x - bodyj->pos[0];
                const real dy = y - bodyj->pos[1];
                const real dz = z - bodyj->pos[2];
                const real dvx = vx - bodyj->vel[0];
                const real dvy = vy - bodyj->vel[1];
                const real dvz = vz - bodyj->vel[2];
                const real rRel = sqrt(dx * dx + dy * dy + dz * dz);
                const real dPosDotVel = dx * dvx + dy * dvy + dz * dvz;
                const real dVelDotVel = dvx * dvx + dvy * dvy + dvz * dvz;
                // 1st order PPN approximation, equation 4-61 from Moyer (2003),
                // https://descanso.jpl.nasa.gov/monograph/series2/Descanso2_all.pdf
                const real fac1 = gm / (c2 * rRel * rRel * rRel);
                const real fac2 =
                    (2 * (beta + gamma) * gm / rRel - gamma * dVelDotVel);
                const real fac3 = 2 * (1 + gamma) * dPosDotVel;
                accInteg[3 * i + 0] += fac1 * (fac2 * dx + fac3 * dvx);
                accInteg[3 * i + 1] += fac1 * (fac2 * dy + fac3 * dvy);
                accInteg[3 * i + 2] += fac1 * (fac2 * dz + fac3 * dvz);
                if (propSim->integBodies[i].propStm) {
                    stm_ppn_simple(propSim->integBodies[i], gm, c, beta, gamma, dx,
                                   dy, dz, dvx, dvy, dvz, 3 * i + 3, accInteg);
                }
                #ifdef PRINT_FORCES
                forceFile << std::setw(10) << bodyj->spiceId << std::setw(25)
                          << G * massj << std::setw(25)
                          << fac1 * (fac2 * dx + fac3 * dvx) << std::setw(25)
                          << fac1 * (fac2 * dy + fac3 * dvy) << std::setw(25)
                          << fac1 * (fac2 * dz + fac3 * dvz) << std::endl;
                #endif
            }
        }
    }
    #ifdef PRINT_FORCES
    forceFile.close();
    #endif
}

void force_ppn_eih(const propSimulation *propSim, std::vector<real> &accInteg) {
// calculate accelerations using the Einstein-Infeld-Hoffmann (EIH) PPN
// formalism see eqn 27 in
// https://iopscience.iop.org/article/10.3847/1538-3881/abd414/pdf (without
// the factor of 1 in the first big summation)
    #ifdef PRINT_FORCES
    std::ofstream forceFile;
    forceFile.precision(16);
    forceFile.setf(std::ios::scientific);
    forceFile.setf(std::ios::right, std::ios::adjustfield);
    forceFile.open("cpp.11", std::ios::app);
    #endif
    const real G = propSim->consts.G;
    const real c = propSim->consts.clight;
    const real oneOverC2 = 1.0 / (c * c);
    const real beta = 1.0;
    const real gamma = 1.0;
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        const real xi = propSim->integBodies[i].pos[0];
        const real yi = propSim->integBodies[i].pos[1];
        const real zi = propSim->integBodies[i].pos[2];
        const real vxi = propSim->integBodies[i].vel[0];
        const real vyi = propSim->integBodies[i].vel[1];
        const real vzi = propSim->integBodies[i].vel[2];
        real axi = 0.0;
        real ayi = 0.0;
        real azi = 0.0;
        for (size_t j = 0; j < propSim->integParams.nTotal; j++) {
            const Body *bodyj;
            if (j < propSim->integParams.nInteg) {
                bodyj = &propSim->integBodies[j];
            } else {
                bodyj = &propSim->spiceBodies[j - propSim->integParams.nInteg];
            }
            const real massj = bodyj->mass;
            if (i != j && massj != 0.0 && bodyj->isPPN) {
                const real muj = G * massj;
                const real xj = bodyj->pos[0];
                const real yj = bodyj->pos[1];
                const real zj = bodyj->pos[2];
                const real vxj = bodyj->vel[0];
                const real vyj = bodyj->vel[1];
                const real vzj = bodyj->vel[2];
                const real dxij = xi - xj;
                const real dyij = yi - yj;
                const real dzij = zi - zj;
                const real dvxij = vxi - vxj;
                const real dvyij = vyi - vyj;
                const real dvzij = vzi - vzj;
                const real rRelij =
                    sqrt(dxij * dxij + dyij * dyij + dzij * dzij);
                const real rRelij3 = rRelij * rRelij * rRelij;
                const real term1c = (vxi * vxi + vyi * vyi + vzi * vzi) * oneOverC2;
                const real term1d = (vxj * vxj + vyj * vyj + vzj * vzj) * oneOverC2;
                const real term1e = vxi * vxj + vyi * vyj + vzi * vzj;
                const real rijDotVj = dxij * vxj + dyij * vyj + dzij * vzj;
                const real term1f = rijDotVj * rijDotVj / (rRelij * rRelij);
                real term1a = 0.0;
                real term1b = 0.0;
                real axj = 0.0;
                real ayj = 0.0;
                real azj = 0.0;
                for (size_t k = 0; k < propSim->integParams.nTotal; k++) {
                    const Body *bodyk;
                    if (k < propSim->integParams.nInteg) {
                        bodyk = &propSim->integBodies[k];
                    } else {
                        bodyk =
                            &propSim
                                 ->spiceBodies[k - propSim->integParams.nInteg];
                    }
                    const real massk = bodyk->mass;
                    if (massk != 0.0 && bodyk->isMajor) {
                        const real muk = G * massk;
                        const real xk = bodyk->pos[0];
                        const real yk = bodyk->pos[1];
                        const real zk = bodyk->pos[2];
                        // if (k != i){
                        const real dxik = xi - xk;
                        const real dyik = yi - yk;
                        const real dzik = zi - zk;
                        const real rRelik =
                            sqrt(dxik * dxik + dyik * dyik + dzik * dzik);
                        term1a += muk / rRelik;
                        // }
                        if (k != j) {
                            const real dxjk = xj - xk;
                            const real dyjk = yj - yk;
                            const real dzjk = zj - zk;
                            const real rReljk =
                                sqrt(dxjk * dxjk + dyjk * dyjk + dzjk * dzjk);
                            term1b += muk / rReljk;
                            const real rReljk3 = rReljk * rReljk * rReljk;
                            axj -= muk * dxjk / rReljk3;
                            ayj -= muk * dyjk / rReljk3;
                            azj -= muk * dzjk / rReljk3;
                        }
                    }
                }
                const real rijDotAj = dxij * axj + dyij * ayj + dzij * azj;
                const real term1g = -rijDotAj;
                const real term1Fac = -muj / rRelij3 *
                    (-2.0 * (beta + gamma) * oneOverC2 * term1a -
                     (2.0 * beta - 1) * oneOverC2 * term1b + gamma * term1c +
                     (1.0 + gamma) * term1d -
                     2.0 * (1.0 + gamma) * oneOverC2 * term1e -
                     1.5 * oneOverC2 * term1f + 0.5 * oneOverC2 * term1g);
                const real term1X = term1Fac * dxij;
                const real term1Y = term1Fac * dyij;
                const real term1Z = term1Fac * dzij;
                const real term2DotProduct = dxij *
                        ((2.0 + 2.0 * gamma) * vxi -
                         (1.0 + 2.0 * gamma) * vxj) +
                    dyij *
                        ((2.0 + 2.0 * gamma) * vyi -
                         (1.0 + 2.0 * gamma) * vyj) +
                    dzij *
                        ((2.0 + 2.0 * gamma) * vzi - (1.0 + 2.0 * gamma) * vzj);
                const real term2Fac =
                    oneOverC2 * muj / rRelij3 * term2DotProduct;
                const real term2X = term2Fac * dvxij;
                const real term2Y = term2Fac * dvyij;
                const real term2Z = term2Fac * dvzij;
                const real term3Fac =
                    (3.0 + 4.0 * gamma) * 0.5 * oneOverC2 * muj / rRelij;
                const real term3X = term3Fac * axj;
                const real term3Y = term3Fac * ayj;
                const real term3Z = term3Fac * azj;
                axi += term1X + term2X + term3X;
                ayi += term1Y + term2Y + term3Y;
                azi += term1Z + term2Z + term3Z;
                if (propSim->integBodies[i].propStm && bodyj->spiceId == 10) {
                    stm_ppn_simple(propSim->integBodies[i], muj, c, beta, gamma,
                                   dxij, dyij, dzij, dvxij, dvyij, dvzij,
                                   3 * i + 3, accInteg);
                }
                #ifdef PRINT_FORCES
                forceFile << std::setw(10) << bodyj->spiceId << std::setw(25)
                          << term1X + term2X + term3X << std::setw(25)
                          << term1Y + term2Y + term3Y << std::setw(25)
                          << term1Z + term2Z + term3Z << std::endl;
                #endif
            }
        }
        #ifdef PRINT_FORCES
        forceFile << std::setw(10) << "EIH" << std::setw(25) << axi
                  << std::setw(25) << ayi << std::setw(25) << azi << std::endl;
        #endif
        accInteg[3 * i + 0] += axi;
        accInteg[3 * i + 1] += ayi;
        accInteg[3 * i + 2] += azi;
    }
    #ifdef PRINT_FORCES
    forceFile.close();
    #endif
}

void force_J2(const propSimulation *propSim, std::vector<real> &accInteg) {
    #ifdef PRINT_FORCES
    std::ofstream forceFile;
    forceFile.precision(16);
    forceFile.setf(std::ios::scientific);
    forceFile.setf(std::ios::right, std::ios::adjustfield);
    forceFile.open("cpp.11", std::ios::app);
    #endif
    const real G = propSim->consts.G;
    const real smoothing_threshold = 100.0e3L/propSim->consts.du2m;
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        const real x = propSim->integBodies[i].pos[0];
        const real y = propSim->integBodies[i].pos[1];
        const real z = propSim->integBodies[i].pos[2];
        for (size_t j = 0; j < propSim->integParams.nTotal; j++) {
            const Body *bodyj;
            if (j < propSim->integParams.nInteg) {
                bodyj = &propSim->integBodies[j];
            } else {
                bodyj = &propSim->spiceBodies[j - propSim->integParams.nInteg];
            }
            const real massj = bodyj->mass;
            if (i != j && massj != 0.0 && bodyj->isJ2) {
                const real dx = x - bodyj->pos[0];
                const real dy = y - bodyj->pos[1];
                const real dz = z - bodyj->pos[2];
                const real rRel = sqrt(dx * dx + dy * dy + dz * dz);
                const real rRel2 = rRel * rRel;
                const real rRel5 = rRel2 * rRel2 * rRel;
                const real radius = bodyj->radius;
                const real poleRA = bodyj->poleRA;
                const real sinRA = sin(poleRA);
                const real cosRA = cos(poleRA);
                const real poleDec = bodyj->poleDec;
                const real sinDec = sin(poleDec);
                const real cosDec = cos(poleDec);
                const real dxBody = -dx * sinRA + dy * cosRA;
                const real dyBody =
                    -dx * cosRA * sinDec - dy * sinRA * sinDec + dz * cosDec;
                const real dzBody =
                    dx * cosRA * cosDec + dy * sinRA * cosDec + dz * sinDec;
                const real fac1 =
                    3 * G * massj * bodyj->J2 * radius * radius / (2 * rRel5);
                const real fac2 = 5 * dzBody * dzBody / rRel2 - 1;
                real axBody = fac1 * fac2 * dxBody;
                real ayBody = fac1 * fac2 * dyBody;
                real azBody = fac1 * (fac2 - 2) * dzBody;
                if (rRel <= radius+smoothing_threshold) {
                    const real depth = radius+smoothing_threshold-rRel;
                    real smoothing = cos(PI*depth/(2*smoothing_threshold));
                    if (depth > smoothing_threshold){
                        smoothing = 0.0;
                    }
                    axBody *= smoothing;
                    ayBody *= smoothing;
                    azBody *= smoothing;
                }
                accInteg[3 * i + 0] += -axBody * sinRA -
                    ayBody * cosRA * sinDec + azBody * cosRA * cosDec;
                accInteg[3 * i + 1] += axBody * cosRA -
                    ayBody * sinRA * sinDec + azBody * sinRA * cosDec;
                accInteg[3 * i + 2] += ayBody * cosDec + azBody * sinDec;
                if (propSim->integBodies[i].propStm && rRel < 0.1) {
                    stm_J2(propSim->integBodies[i], G*massj, bodyj->J2, dxBody,
                           dyBody, dzBody, radius, sinRA, cosRA, sinDec, cosDec,
                           smoothing_threshold, 3*i+3, accInteg);
                }
                #ifdef PRINT_FORCES
                forceFile << std::setw(10) << bodyj->spiceId << std::setw(25)
                          << -axBody * sinRA - ayBody * cosRA * sinDec +
                        azBody * cosRA * cosDec
                          << std::setw(25)
                          << axBody * cosRA - ayBody * sinRA * sinDec +
                        azBody * sinRA * cosDec
                          << std::setw(25) << ayBody * cosDec + azBody * sinDec
                          << std::endl;
                #endif
            }
        }
    }
    #ifdef PRINT_FORCES
    forceFile.close();
    #endif
}

void force_nongrav(const propSimulation *propSim, std::vector<real> &accInteg) {
#ifdef PRINT_FORCES
    std::ofstream forceFile;
    forceFile.precision(16);
    forceFile.setf(std::ios::scientific);
    forceFile.setf(std::ios::right, std::ios::adjustfield);
    forceFile.open("cpp.11", std::ios::app);
#endif
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        for (size_t j = 0; j < propSim->integParams.nTotal; j++) {
            const Body *bodyj;
            if (j < propSim->integParams.nInteg) {
                bodyj = &propSim->integBodies[j];
            } else {
                bodyj = &propSim->spiceBodies[j - propSim->integParams.nInteg];
            }
            if (bodyj->spiceId == 10 && propSim->integBodies[i].isNongrav) {
                const real x = propSim->integBodies[i].pos[0];
                const real y = propSim->integBodies[i].pos[1];
                const real z = propSim->integBodies[i].pos[2];
                const real vx = propSim->integBodies[i].vel[0];
                const real vy = propSim->integBodies[i].vel[1];
                const real vz = propSim->integBodies[i].vel[2];
                const real a1 = propSim->integBodies[i].ngParams.a1;
                const real a2 = propSim->integBodies[i].ngParams.a2;
                const real a3 = propSim->integBodies[i].ngParams.a3;
                const real alpha = propSim->integBodies[i].ngParams.alpha;
                const real k = propSim->integBodies[i].ngParams.k;
                const real m = propSim->integBodies[i].ngParams.m;
                const real n = propSim->integBodies[i].ngParams.n;
                const real r0 = propSim->integBodies[i].ngParams.r0_au;
                const real dx = x - bodyj->pos[0];
                const real dy = y - bodyj->pos[1];
                const real dz = z - bodyj->pos[2];
                const real dvx = vx - bodyj->vel[0];
                const real dvy = vy - bodyj->vel[1];
                const real dvz = vz - bodyj->vel[2];
                const real rRel = sqrt(dx * dx + dy * dy + dz * dz);
                const real g =
                    alpha * pow(rRel / r0, -m) * pow(1 + pow(rRel / r0, n), -k);
                real *dpos = new real[3];
                dpos[0] = dx;
                dpos[1] = dy;
                dpos[2] = dz;
                real *dvel = new real[3];
                dvel[0] = dvx;
                dvel[1] = dvy;
                dvel[2] = dvz;
                real *hRelVec = new real[3];
                memset(hRelVec, 0.0, 3 * sizeof(real));
                real *eRHat = new real[3];
                memset(eRHat, 0.0, 3 * sizeof(real));
                real *eTHat = new real[3];
                memset(eTHat, 0.0, 3 * sizeof(real));
                real *eNHat = new real[3];
                memset(eNHat, 0.0, 3 * sizeof(real));
                vunit(dpos, (size_t)3, eRHat);
                vcross(dpos, dvel, hRelVec);
                vunit(hRelVec, (size_t)3, eNHat);
                vcross(eNHat, eRHat, eTHat);
                accInteg[3 * i + 0] +=
                    g * (a1 * eRHat[0] + a2 * eTHat[0] + a3 * eNHat[0]);
                accInteg[3 * i + 1] +=
                    g * (a1 * eRHat[1] + a2 * eTHat[1] + a3 * eNHat[1]);
                accInteg[3 * i + 2] +=
                    g * (a1 * eRHat[2] + a2 * eTHat[2] + a3 * eNHat[2]);
                if (propSim->integBodies[i].propStm) {
                    stm_nongrav(propSim->integBodies[i], g,
                                propSim->integBodies[i].ngParams, dx, dy, dz,
                                dvx, dvy, dvz, dpos, hRelVec, 3 * i + 3,
                                accInteg);
                }
                #ifdef PRINT_FORCES
                forceFile << std::setw(10) << bodyj->spiceId << std::setw(25)
                          << g * (a1 * eRHat[0] + a2 * eTHat[0] + a3 * eNHat[0])
                          << std::setw(25)
                          << g * (a1 * eRHat[1] + a2 * eTHat[1] + a3 * eNHat[1])
                          << std::setw(25)
                          << g * (a1 * eRHat[2] + a2 * eTHat[2] + a3 * eNHat[2])
                          << std::endl;
                #endif
                delete[] dpos;
                delete[] dvel;
                delete[] hRelVec;
                delete[] eRHat;
                delete[] eTHat;
                delete[] eNHat;
            }
        }
    }
    #ifdef PRINT_FORCES
    forceFile.close();
    #endif
}

void force_thruster(const propSimulation *propSim,
                    std::vector<real> &accInteg) {
#ifdef PRINT_FORCES
    std::ofstream forceFile;
    forceFile.precision(16);
    forceFile.setf(std::ios::scientific);
    forceFile.setf(std::ios::right, std::ios::adjustfield);
    forceFile.open("cpp.11", std::ios::app);
#endif
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        if (propSim->integBodies[i].isThrusting) {
            real *vel = new real[3];
            vel[0] = propSim->integBodies[i].vel[0];
            vel[1] = propSim->integBodies[i].vel[1];
            vel[2] = propSim->integBodies[i].vel[2];
            real *vHat = new real[3];
            memset(vHat, 0, 3 * sizeof(real));
            const real acc_thruster =
                1.0e7L / propSim->consts.du2m;  // m/day^2 -> au/day^2
            vunit(vel, (size_t)3, vHat);
            accInteg[3 * i + 0] += acc_thruster * vHat[0];
            accInteg[3 * i + 1] += acc_thruster * vHat[1];
            accInteg[3 * i + 2] += acc_thruster * vHat[2];
            #ifdef PRINT_FORCES
            forceFile << "THRUSTER " << acc_thruster * vHat[0] << std::setw(25)
                      << acc_thruster * vHat[1] << std::setw(25)
                      << acc_thruster * vHat[2] << std::endl;
            #endif
            delete[] vel;
            delete[] vHat;
        }
    }
    #ifdef PRINT_FORCES
    forceFile.close();
    #endif
}
