#include "force.h"
// #define PRINT_FORCES 1
#ifdef PRINT_FORCES
#include <iomanip>
#include <fstream>
#endif

std::vector<real> get_state_der(const real &t, const std::vector<real> &xInteg,
                                propSimulation *propSim) {
    real *accInteg = new real[3 * propSim->integParams.nInteg];
    memset(accInteg, 0.0, 3 * propSim->integParams.nInteg * sizeof(real));
    for (size_t i = 0; i < propSim->integParams.nInteg; i++){
        propSim->integBodies[i].pos[0] = xInteg[6*i];
        propSim->integBodies[i].pos[1] = xInteg[6*i+1];
        propSim->integBodies[i].pos[2] = xInteg[6*i+2];
        propSim->integBodies[i].vel[0] = xInteg[6*i+3];
        propSim->integBodies[i].vel[1] = xInteg[6*i+4];
        propSim->integBodies[i].vel[2] = xInteg[6*i+5];
    }
    for (size_t i = 0; i < propSim->integParams.nSpice; i++) {
        double xSpice[9];
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
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        propSim->integBodies[i].acc[0] = accInteg[3*i];
        propSim->integBodies[i].acc[1] = accInteg[3*i+1];
        propSim->integBodies[i].acc[2] = accInteg[3*i+2];
    }
    // initialize acceleration vector in std::vector form
    std::vector<real> accIntegVec(3 * propSim->integParams.nInteg, 0.0);
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        accIntegVec[3*i] = accInteg[3*i];
        accIntegVec[3*i+1] = accInteg[3*i+1];
        accIntegVec[3*i+2] = accInteg[3*i+2];
    }
    delete[] accInteg;
    return accIntegVec;
}

void force_newton(const propSimulation *propSim, real *accInteg) {
    #ifdef PRINT_FORCES
    std::ofstream forceFile;
    forceFile.precision(16);
    forceFile.setf(std::ios::scientific);
    forceFile.setf(std::ios::right, std::ios::adjustfield);
    forceFile.open("cpp.11", std::ios::app);
    #endif
    real G = propSim->consts.G;
    real x, y, z;
    real dx, dy, dz;
    real rRel, rRel3;
    real ax, ay, az;
    real massj;
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        x = propSim->integBodies[i].pos[0];
        y = propSim->integBodies[i].pos[1];
        z = propSim->integBodies[i].pos[2];
        ax = 0.0;
        ay = 0.0;
        az = 0.0;
        for (size_t j = 0; j < propSim->integParams.nTotal; j++) {
            massj = propSim->forceParams.masses[j];
            if (i != j && massj != 0.0) {
                const Body *bodyj;
                if (j < propSim->integParams.nInteg) {
                    bodyj = &propSim->integBodies[j];
                } else {
                    bodyj = &propSim->spiceBodies[j-propSim->integParams.nInteg];
                }
                dx = x - bodyj->pos[0];
                dy = y - bodyj->pos[1];
                dz = z - bodyj->pos[2];
                rRel = sqrt(dx * dx + dy * dy + dz * dz);
                rRel3 = rRel * rRel * rRel;
                ax -= G * massj * dx / rRel3;
                ay -= G * massj * dy / rRel3;
                az -= G * massj * dz / rRel3;
                #ifdef PRINT_FORCES
                forceFile << std::setw(10) << propSim->forceParams.spiceIdList[j]
                          << std::setw(25) << G * massj << std::setw(25) << dx
                          << std::setw(25) << dy << std::setw(25) << dz
                          << std::setw(25) << -G * massj * dx / rRel3
                          << std::setw(25) << -G * massj * dy / rRel3
                          << std::setw(25) << -G * massj * dz / rRel3
                          << std::endl;
                #endif
            }
        }
        accInteg[3 * i + 0] += ax;
        accInteg[3 * i + 1] += ay;
        accInteg[3 * i + 2] += az;
    }
    #ifdef PRINT_FORCES
    forceFile.close();
    #endif
}

void force_ppn_simple(const propSimulation *propSim, real *accInteg) {
    #ifdef PRINT_FORCES
    std::ofstream forceFile;
    forceFile.precision(16);
    forceFile.setf(std::ios::scientific);
    forceFile.setf(std::ios::right, std::ios::adjustfield);
    forceFile.open("cpp.11", std::ios::app);
    #endif
    real G = propSim->consts.G;
    real c = propSim->consts.clight;
    real c2 = c * c;
    real x, y, z;
    real dx, dy, dz;
    real rRel, rRel3;
    real vx, vy, vz;
    real dvx, dvy, dvz;
    real ax, ay, az;
    real massj;
    real dPosDotVel, dVelDotVel;
    real gm, gmOverC2;
    real fac1, fac2, fac3;
    real beta = 1.0L;
    real gamma = 1.0L;
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        x = propSim->integBodies[i].pos[0];
        y = propSim->integBodies[i].pos[1];
        z = propSim->integBodies[i].pos[2];
        vx = propSim->integBodies[i].vel[0];
        vy = propSim->integBodies[i].vel[1];
        vz = propSim->integBodies[i].vel[2];
        ax = 0.0;
        ay = 0.0;
        az = 0.0;
        for (size_t j = 0; j < propSim->integParams.nTotal; j++) {
            massj = propSim->forceParams.masses[j];
            if (i != j && massj != 0.0 && propSim->forceParams.spiceIdList[j] == 10) {
                gm = G * massj;
                gmOverC2 = gm / c2;
                const Body *bodyj;
                if (j < propSim->integParams.nInteg) {
                    bodyj = &propSim->integBodies[j];
                } else {
                    bodyj = &propSim->spiceBodies[j-propSim->integParams.nInteg];
                }
                dx = x - bodyj->pos[0];
                dy = y - bodyj->pos[1];
                dz = z - bodyj->pos[2];
                dvx = vx - bodyj->vel[0];
                dvy = vy - bodyj->vel[1];
                dvz = vz - bodyj->vel[2];
                rRel = sqrt(dx * dx + dy * dy + dz * dz);
                rRel3 = rRel * rRel * rRel;
                dPosDotVel = dx * dvx + dy * dvy + dz * dvz;
                dVelDotVel = dvx * dvx + dvy * dvy + dvz * dvz;
                // 1st order PPN approximation, equation 4-61 from Moyer (2003),
                // https://descanso.jpl.nasa.gov/monograph/series2/Descanso2_all.pdf
                fac1 = gmOverC2 / rRel3;
                fac2 = (2 * (beta + gamma) * gm / rRel - gamma * dVelDotVel);
                fac3 = 2 * (1 + gamma) * dPosDotVel;
                ax += fac1 * (fac2 * dx + fac3 * dvx);
                ay += fac1 * (fac2 * dy + fac3 * dvy);
                az += fac1 * (fac2 * dz + fac3 * dvz);
                #ifdef PRINT_FORCES
                forceFile << std::setw(10) << propSim->forceParams.spiceIdList[j]
                          << std::setw(25) << G * massj << std::setw(25)
                          << fac1 * (fac2 * dx + fac3 * dvx) << std::setw(25)
                          << fac1 * (fac2 * dy + fac3 * dvy) << std::setw(25)
                          << fac1 * (fac2 * dz + fac3 * dvz) << std::endl;
                #endif
            }
        }
        accInteg[3 * i + 0] += ax;
        accInteg[3 * i + 1] += ay;
        accInteg[3 * i + 2] += az;
    }
    #ifdef PRINT_FORCES
    forceFile.close();
    #endif
}

void force_ppn_eih(const propSimulation *propSim, real *accInteg) {
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
    real G = propSim->consts.G;
    real c2 = propSim->consts.clight * propSim->consts.clight;
    real oneOverC2 = 1.0 / c2;
    real beta = 1.0;
    real gamma = 1.0;
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
        real axj, ayj, azj;
        for (size_t j = 0; j < propSim->integParams.nTotal; j++) {
            const real massj = propSim->forceParams.masses[j];
            if (i != j && massj != 0.0 && propSim->forceParams.isPPNList[j]) {
                const Body *bodyj;
                if (j < propSim->integParams.nInteg) {
                    bodyj = &propSim->integBodies[j];
                } else {
                    bodyj = &propSim->spiceBodies[j-propSim->integParams.nInteg];
                }
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
                const real viDotVi = vxi * vxi + vyi * vyi + vzi * vzi;
                const real term1c = viDotVi * oneOverC2;
                const real vjDotVj = vxj * vxj + vyj * vyj + vzj * vzj;
                const real term1d = vjDotVj * oneOverC2;
                const real viDotVj = vxi * vxj + vyi * vyj + vzi * vzj;
                const real term1e = viDotVj;
                const real rijDotVj = dxij * vxj + dyij * vyj + dzij * vzj;
                const real term1f = rijDotVj * rijDotVj / (rRelij * rRelij);
                real term1a = 0.0;
                real term1b = 0.0;
                axj = 0.0;
                ayj = 0.0;
                azj = 0.0;
                for (size_t k = 0; k < propSim->integParams.nTotal; k++) {
                    const real massk = propSim->forceParams.masses[k];
                    if (massk != 0.0 && propSim->forceParams.isMajorList[k]) {
                        const real muk = G * massk;
                        const Body *bodyk;
                        if (k < propSim->integParams.nInteg) {
                            bodyk = &propSim->integBodies[k];
                        } else {
                            bodyk = &propSim->spiceBodies[k-propSim->integParams.nInteg];
                        }
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

                const real term1Fac =
                    -muj / rRelij3 * (
                        -2.0 * (beta + gamma) * oneOverC2 * term1a -
                        (2.0 * beta - 1) * oneOverC2 * term1b + gamma * term1c +
                        (1.0 + gamma) * term1d -
                        2.0 * (1.0 + gamma) * oneOverC2 * term1e -
                        1.5 * oneOverC2 * term1f + 0.5 * oneOverC2 * term1g
                    );
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

                const real term3Fac = (3.0 + 4.0 * gamma) * 0.5 * oneOverC2 * muj / rRelij;
                const real term3X = term3Fac * axj;
                const real term3Y = term3Fac * ayj;
                const real term3Z = term3Fac * azj;

                axi += term1X + term2X + term3X;
                ayi += term1Y + term2Y + term3Y;
                azi += term1Z + term2Z + term3Z;
                #ifdef PRINT_FORCES
                forceFile << std::setw(10) << propSim->forceParams.spiceIdList[j]
                          << std::setw(25) << term1X + term2X + term3X
                          << std::setw(25) << term1Y + term2Y + term3Y
                          << std::setw(25) << term1Z + term2Z + term3Z
                          << std::endl;
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

void force_J2(const propSimulation *propSim, real *accInteg) {
    #ifdef PRINT_FORCES
    std::ofstream forceFile;
    forceFile.precision(16);
    forceFile.setf(std::ios::scientific);
    forceFile.setf(std::ios::right, std::ios::adjustfield);
    forceFile.open("cpp.11", std::ios::app);
    #endif
    real G = propSim->consts.G;
    real x, y, z;
    real dx, dy, dz;
    real dxBody, dyBody, dzBody;
    real rRel, rRel2, rRel5;
    real radius;
    real ax, ay, az;
    real massj, poleRA, poleDec;
    real sinRA, cosRA, sinDec, cosDec;
    std::vector<std::vector<real>> R1(3, std::vector<real>(3));
    std::vector<std::vector<real>> R2(3, std::vector<real>(3));
    std::vector<std::vector<real>> R(3, std::vector<real>(3));
    std::vector<std::vector<real>> Rinv(3, std::vector<real>(3));
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        x = propSim->integBodies[i].pos[0];
        y = propSim->integBodies[i].pos[1];
        z = propSim->integBodies[i].pos[2];
        ax = 0.0;
        ay = 0.0;
        az = 0.0;
        for (size_t j = 0; j < propSim->integParams.nTotal; j++) {
            massj = propSim->forceParams.masses[j];
            if (i != j && massj != 0.0 && propSim->forceParams.isJ2List[j]) {
                const Body *bodyj;
                if (j < propSim->integParams.nInteg) {
                    bodyj = &propSim->integBodies[j];
                } else {
                    bodyj = &propSim->spiceBodies[j-propSim->integParams.nInteg];
                }
                dx = x - bodyj->pos[0];
                dy = y - bodyj->pos[1];
                dz = z - bodyj->pos[2];
                rRel = sqrt(dx * dx + dy * dy + dz * dz);
                rRel2 = rRel * rRel;
                rRel5 = rRel2 * rRel2 * rRel;
                radius = propSim->forceParams.radii[j];
                poleRA = propSim->forceParams.poleRAList[j];
                sinRA = sin(poleRA);
                cosRA = cos(poleRA);
                poleDec = propSim->forceParams.poleDecList[j];
                sinDec = sin(poleDec);
                cosDec = cos(poleDec);
                dxBody = -dx * sinRA + dy * cosRA;
                dyBody =
                    -dx * cosRA * sinDec - dy * sinRA * sinDec + dz * cosDec;
                dzBody =
                    dx * cosRA * cosDec + dy * sinRA * cosDec + dz * sinDec;
                real fac1 = 3 * G * massj * propSim->forceParams.J2List[j] * radius *
                    radius / (2 * rRel5);
                real fac2 = 5 * dzBody * dzBody / rRel2 - 1;
                real axBody = fac1 * fac2 * dxBody;
                real ayBody = fac1 * fac2 * dyBody;
                real azBody = fac1 * (fac2 - 2) * dzBody;
                real axEquat = -axBody * sinRA - ayBody * cosRA * sinDec +
                    azBody * cosRA * cosDec;
                real ayEquat = axBody * cosRA - ayBody * sinRA * sinDec +
                    azBody * sinRA * cosDec;
                real azEquat = ayBody * cosDec + azBody * sinDec;
                ax += axEquat;
                ay += ayEquat;
                az += azEquat;
                #ifdef PRINT_FORCES
                forceFile << std::setw(10) << propSim->forceParams.spiceIdList[j]
                          << std::setw(25) << axEquat << std::setw(25)
                          << ayEquat << std::setw(25) << azEquat << std::endl;
                #endif
            }
        }
        accInteg[3 * i + 0] += ax;
        accInteg[3 * i + 1] += ay;
        accInteg[3 * i + 2] += az;
    }
    #ifdef PRINT_FORCES
    forceFile.close();
    #endif
}

void force_nongrav(const propSimulation *propSim, real *accInteg) {
    #ifdef PRINT_FORCES
    std::ofstream forceFile;
    forceFile.precision(16);
    forceFile.setf(std::ios::scientific);
    forceFile.setf(std::ios::right, std::ios::adjustfield);
    forceFile.open("cpp.11", std::ios::app);
    #endif
    real a1, a2, a3;
    real alpha, k, m, n;
    real r0;
    real g;
    real x, y, z;
    real vx, vy, vz;
    real dx, dy, dz;
    real dvx, dvy, dvz;
    real rRel;
    std::vector<real> hRelVec(3), hRelHatVec(3);
    std::vector<real> eRHat(3), eTHat(3), eNHat(3);
    real ax, ay, az;
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        x = propSim->integBodies[i].pos[0];
        y = propSim->integBodies[i].pos[1];
        z = propSim->integBodies[i].pos[2];
        vx = propSim->integBodies[i].vel[0];
        vy = propSim->integBodies[i].vel[1];
        vz = propSim->integBodies[i].vel[2];
        ax = 0.0;
        ay = 0.0;
        az = 0.0;
        for (size_t j = 0; j < propSim->integParams.nTotal; j++) {
            if (propSim->forceParams.spiceIdList[j] == 10 &&
                propSim->forceParams.isNongravList[i]) {  // j is Sun idx in spiceIdList
                a1 = propSim->forceParams.ngParamsList[i].a1;
                a2 = propSim->forceParams.ngParamsList[i].a2;
                a3 = propSim->forceParams.ngParamsList[i].a3;
                alpha = propSim->forceParams.ngParamsList[i].alpha;
                k = propSim->forceParams.ngParamsList[i].k;
                m = propSim->forceParams.ngParamsList[i].m;
                n = propSim->forceParams.ngParamsList[i].n;
                r0 = propSim->forceParams.ngParamsList[i].r0_au * 1.495978707e11 /
                    propSim->consts.du2m;
                const Body *bodyj;
                if (j < propSim->integParams.nInteg) {
                    bodyj = &propSim->integBodies[j];
                } else {
                    bodyj = &propSim->spiceBodies[j-propSim->integParams.nInteg];
                }
                dx = x - bodyj->pos[0];
                dy = y - bodyj->pos[1];
                dz = z - bodyj->pos[2];
                dvx = vx - bodyj->vel[0];
                dvy = vy - bodyj->vel[1];
                dvz = vz - bodyj->vel[2];
                rRel = sqrt(dx * dx + dy * dy + dz * dz);
                g = alpha * pow(rRel / r0, -m) * pow(1 + pow(rRel / r0, n), -k);
                vunit({dx, dy, dz}, eRHat);
                vcross({dx, dy, dz}, {dvx, dvy, dvz}, hRelVec);
                vunit(hRelVec, eNHat);
                vcross(eNHat, eRHat, eTHat);
                ax += g * (a1 * eRHat[0] + a2 * eTHat[0] + a3 * eNHat[0]);
                ay += g * (a1 * eRHat[1] + a2 * eTHat[1] + a3 * eNHat[1]);
                az += g * (a1 * eRHat[2] + a2 * eTHat[2] + a3 * eNHat[2]);
                #ifdef PRINT_FORCES
                forceFile << std::setw(10) << propSim->forceParams.spiceIdList[j]
                          << std::setw(25)
                          << g * (a1 * eRHat[0] + a2 * eTHat[0] + a3 * eNHat[0])
                          << std::setw(25)
                          << g * (a1 * eRHat[1] + a2 * eTHat[1] + a3 * eNHat[1])
                          << std::setw(25)
                          << g * (a1 * eRHat[2] + a2 * eTHat[2] + a3 * eNHat[2])
                          << std::endl;
                #endif
            }
        }
        accInteg[3 * i + 0] += ax;
        accInteg[3 * i + 1] += ay;
        accInteg[3 * i + 2] += az;
    }
    #ifdef PRINT_FORCES
    forceFile.close();
    #endif
}

void force_thruster(const propSimulation *propSim, real *accInteg) {
#ifdef PRINT_FORCES
    std::ofstream forceFile;
    forceFile.precision(16);
    forceFile.setf(std::ios::scientific);
    forceFile.setf(std::ios::right, std::ios::adjustfield);
    forceFile.open("cpp.11", std::ios::app);
#endif
    real vx, vy, vz;
    real ax, ay, az;
    std::vector<real> vHat(3);
    real acc_thruster = 1.0e7L /propSim->consts.du2m;  // m/day^2 -> au/day^2
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        vx = propSim->integBodies[i].vel[0];
        vy = propSim->integBodies[i].vel[1];
        vz = propSim->integBodies[i].vel[2];
        ax = 0.0;
        ay = 0.0;
        az = 0.0;
        if (propSim->forceParams.isThrustingList[i]) {
            vunit({vx, vy, vz}, vHat);
            ax += acc_thruster * vHat[0];
            ay += acc_thruster * vHat[1];
            az += acc_thruster * vHat[2];
            #ifdef PRINT_FORCES
            forceFile << "THRUSTER " << acc_thruster * vHat[0] << std::setw(25)
                      << acc_thruster * vHat[1] << std::setw(25)
                      << acc_thruster * vHat[2] << std::endl;
            #endif
        }
        accInteg[3 * i + 0] += ax;
        accInteg[3 * i + 1] += ay;
        accInteg[3 * i + 2] += az;
    }
    #ifdef PRINT_FORCES
    forceFile.close();
    #endif
}
