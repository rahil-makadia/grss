#include "force.h"
// #define PRINT_FORCES 1
#ifdef PRINT_FORCES
#include <iomanip>
#include <fstream>
#endif

std::vector<real> get_state_der(const real &t, const std::vector<real> &xInteg,
                                const ForceParameters &forceParams,
                                const IntegrationParameters &integParams,
                                const Constants &consts) {
    std::vector<real> posAll(3 * integParams.nTotal, 0.0);
    std::vector<real> velAll(3 * integParams.nTotal, 0.0);
    std::vector<real> accInteg(3 * integParams.nInteg, 0.0);
    for (size_t i = 0; i < integParams.nInteg; i++) {
        posAll[3 * i] = xInteg[6 * i];
        posAll[3 * i + 1] = xInteg[6 * i + 1];
        posAll[3 * i + 2] = xInteg[6 * i + 2];
        velAll[3 * i] = xInteg[6 * i + 3];
        velAll[3 * i + 1] = xInteg[6 * i + 4];
        velAll[3 * i + 2] = xInteg[6 * i + 5];
    }
    for (size_t i = integParams.nInteg; i < integParams.nTotal; i++) {
        double xSpice_i[6];
        double lt;
        get_spice_state_lt(forceParams.spiceIdList[i], t, consts, xSpice_i, lt);
        posAll[3 * i] = xSpice_i[0];
        posAll[3 * i + 1] = xSpice_i[1];
        posAll[3 * i + 2] = xSpice_i[2];
        velAll[3 * i] = xSpice_i[3];
        velAll[3 * i + 1] = xSpice_i[4];
        velAll[3 * i + 2] = xSpice_i[5];
    }
#ifdef PRINT_FORCES
    std::ofstream forceFile;
    forceFile.precision(8);
    forceFile.setf(std::ios::scientific);
    forceFile.setf(std::ios::right, std::ios::adjustfield);
    forceFile.open("cpp.11", std::ios::app);
    forceFile << "time (MJD): " << std::setw(16) << t << " state:";
    forceFile.precision(16);
    for (size_t i = 0; i<6; i++){
        forceFile << std::setw(25) << xInteg[i];
    }
    forceFile << std::endl;
    forceFile.close();
#endif
    force_newton(posAll, accInteg, forceParams, integParams, consts);
    // force_ppn_simple(posAll, velAll, accInteg, forceParams, integParams,
    // consts);
    force_ppn_eih(posAll, velAll, accInteg, forceParams, integParams, consts);
    force_J2(posAll, accInteg, forceParams, integParams, consts);
    force_nongrav(posAll, velAll, accInteg, forceParams, integParams, consts);
    force_thruster(velAll, accInteg, forceParams, integParams, consts);
#ifdef PRINT_FORCES
    forceFile.open("cpp.11", std::ios::app);
    forceFile << std::setw(10) << "total acc:" << std::setw(25) << accInteg[0] << std::setw(25) << accInteg[1] << std::setw(25)
              << accInteg[2] << std::endl;
    forceFile.close();
#endif
    return accInteg;
}

void force_newton(const std::vector<real> &posAll, std::vector<real> &accInteg,
                  const ForceParameters &forceParams,
                  const IntegrationParameters &integParams,
                  const Constants &consts) {
#ifdef PRINT_FORCES
    std::ofstream forceFile;
    forceFile.precision(16);
    forceFile.setf(std::ios::scientific);
    forceFile.setf(std::ios::right, std::ios::adjustfield);
    forceFile.open("cpp.11", std::ios::app);
#endif
    real G = consts.G;
    real x, y, z;
    real dx, dy, dz;
    real rRel, rRel3;
    real ax, ay, az;
    real massj;
    for (size_t i = 0; i < integParams.nInteg; i++) {
        x = posAll[3 * i];
        y = posAll[3 * i + 1];
        z = posAll[3 * i + 2];
        ax = 0.0;
        ay = 0.0;
        az = 0.0;
        for (size_t j = 0; j < integParams.nTotal; j++) {
            massj = forceParams.masses[j];
            if (i != j && massj != 0.0) {
                dx = x - posAll[3 * j];
                dy = y - posAll[3 * j + 1];
                dz = z - posAll[3 * j + 2];
                rRel = sqrt(dx * dx + dy * dy + dz * dz);
                rRel3 = rRel * rRel * rRel;
                ax -= G * massj * dx / rRel3;
                ay -= G * massj * dy / rRel3;
                az -= G * massj * dz / rRel3;
#ifdef PRINT_FORCES
                forceFile << std::setw(10) << forceParams.spiceIdList[j] << std::setw(25) << G * massj
                          << std::setw(25) << dx << std::setw(25) << dy << std::setw(25) << dz << std::setw(25)
                          << -G * massj * dx / rRel3 << std::setw(25)
                          << -G * massj * dy / rRel3 << std::setw(25)
                          << -G * massj * dz / rRel3 << std::endl;
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

void force_ppn_simple(const std::vector<real> &posAll,
                      const std::vector<real> &velAll,
                      std::vector<real> &accInteg,
                      const ForceParameters &forceParams,
                      const IntegrationParameters &integParams,
                      const Constants &consts) {
#ifdef PRINT_FORCES
    std::ofstream forceFile;
    forceFile.precision(16);
    forceFile.setf(std::ios::scientific);
    forceFile.setf(std::ios::right, std::ios::adjustfield);
    forceFile.open("cpp.11", std::ios::app);
#endif
    real G = consts.G;
    real c = consts.clight;
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
    for (size_t i = 0; i < integParams.nInteg; i++) {
        x = posAll[3 * i];
        y = posAll[3 * i + 1];
        z = posAll[3 * i + 2];
        vx = velAll[3 * i];
        vy = velAll[3 * i + 1];
        vz = velAll[3 * i + 2];
        ax = 0.0;
        ay = 0.0;
        az = 0.0;
        for (size_t j = 0; j < integParams.nTotal; j++) {
            massj = forceParams.masses[j];
            gm = G * massj;
            gmOverC2 = gm / c2;
            if (i != j && massj != 0.0 && forceParams.spiceIdList[j] == 10) {
                dx = x - posAll[3 * j];
                dy = y - posAll[3 * j + 1];
                dz = z - posAll[3 * j + 2];
                rRel = sqrt(dx * dx + dy * dy + dz * dz);
                rRel3 = rRel * rRel * rRel;
                dvx = vx - velAll[3 * j];
                dvy = vy - velAll[3 * j + 1];
                dvz = vz - velAll[3 * j + 2];
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
                forceFile << std::setw(10) << forceParams.spiceIdList[j] << std::setw(25) << G * massj << std::setw(25)
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

void force_ppn_eih(const std::vector<real> &posAll,
                   const std::vector<real> &velAll, std::vector<real> &accInteg,
                   const ForceParameters &forceParams,
                   const IntegrationParameters &integParams,
                   const Constants &consts) {
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
    real G = consts.G;
    real c2 = consts.clight * consts.clight;
    real oneOverC2 = 1.0 / c2;
    real beta = 1.0;
    real gamma = 1.0;
    for (size_t i = 0; i < integParams.nInteg; i++) {
        const real xi = posAll[3 * i];
        const real yi = posAll[3 * i + 1];
        const real zi = posAll[3 * i + 2];
        const real vxi = velAll[3 * i];
        const real vyi = velAll[3 * i + 1];
        const real vzi = velAll[3 * i + 2];
        real axi = 0.0;
        real ayi = 0.0;
        real azi = 0.0;
        real axj, ayj, azj;
        for (size_t j = 0; j < integParams.nTotal; j++) {
            const real massj = forceParams.masses[j];
            if (i != j && massj != 0.0 && forceParams.isPPNList[j]) {
                const real muj = G * massj;
                const real xj = posAll[3 * j];
                const real yj = posAll[3 * j + 1];
                const real zj = posAll[3 * j + 2];
                const real vxj = velAll[3 * j];
                const real vyj = velAll[3 * j + 1];
                const real vzj = velAll[3 * j + 2];
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
                for (size_t k = 0; k < integParams.nTotal; k++) {
                    const real massk = forceParams.masses[k];
                    if (massk != 0.0 && forceParams.isMajorList[k]) {
                        const real muk = G * massk;
                        const real xk = posAll[3 * k];
                        const real yk = posAll[3 * k + 1];
                        const real zk = posAll[3 * k + 2];
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
                    -2.0 * (beta + gamma) * oneOverC2 * term1a -
                    (2.0 * beta - 1) * oneOverC2 * term1b + gamma * term1c +
                    (1.0 + gamma) * term1d -
                    2.0 * (1.0 + gamma) * oneOverC2 * term1e -
                    1.5 * oneOverC2 * term1f + 0.5 * oneOverC2 * term1g;
                const real term1X = -muj / rRelij3 * dxij * term1Fac;
                const real term1Y = -muj / rRelij3 * dyij * term1Fac;
                const real term1Z = -muj / rRelij3 * dzij * term1Fac;

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

                const real term3X =
                    (3.0 + 4.0 * gamma) * 0.5 * oneOverC2 * muj / rRelij * axj;
                const real term3Y =
                    (3.0 + 4.0 * gamma) * 0.5 * oneOverC2 * muj / rRelij * ayj;
                const real term3Z =
                    (3.0 + 4.0 * gamma) * 0.5 * oneOverC2 * muj / rRelij * azj;

                axi += term1X + term2X + term3X;
                ayi += term1Y + term2Y + term3Y;
                azi += term1Z + term2Z + term3Z;
#ifdef PRINT_FORCES
                forceFile << std::setw(10) << forceParams.spiceIdList[j] << std::setw(25)
                          << term1X + term2X + term3X << std::setw(25)
                          << term1Y + term2Y + term3Y << std::setw(25)
                          << term1Z + term2Z + term3Z << std::endl;
#endif
            }
        }
#ifdef PRINT_FORCES
        forceFile << std::setw(10) << "EIH" << std::setw(25) << axi << std::setw(25)
        << ayi << std::setw(25) << azi << std::endl;
#endif
        accInteg[3 * i + 0] += axi;
        accInteg[3 * i + 1] += ayi;
        accInteg[3 * i + 2] += azi;
    }
#ifdef PRINT_FORCES
    forceFile.close();
#endif
}

void force_J2(const std::vector<real> &posAll, std::vector<real> &accInteg,
              const ForceParameters &forceParams,
              const IntegrationParameters &integParams,
              const Constants &consts) {
#ifdef PRINT_FORCES
    std::ofstream forceFile;
    forceFile.precision(16);
    forceFile.setf(std::ios::scientific);
    forceFile.setf(std::ios::right, std::ios::adjustfield);
    forceFile.open("cpp.11", std::ios::app);
#endif
    real G = consts.G;
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
    for (size_t i = 0; i < integParams.nInteg; i++) {
        x = posAll[3 * i];
        y = posAll[3 * i + 1];
        z = posAll[3 * i + 2];
        ax = 0.0;
        ay = 0.0;
        az = 0.0;
        for (size_t j = 0; j < integParams.nTotal; j++) {
            massj = forceParams.masses[j];
            if (i != j && massj != 0.0 && forceParams.isJ2List[j]) {
                dx = x - posAll[3 * j];
                dy = y - posAll[3 * j + 1];
                dz = z - posAll[3 * j + 2];
                rRel = sqrt(dx * dx + dy * dy + dz * dz);
                rRel2 = rRel * rRel;
                rRel5 = rRel2 * rRel2 * rRel;
                radius = forceParams.radii[j];
                poleRA = forceParams.poleRAList[j];
                sinRA = sin(poleRA);
                cosRA = cos(poleRA);
                poleDec = forceParams.poleDecList[j];
                sinDec = sin(poleDec);
                cosDec = cos(poleDec);
                dxBody = -dx*sinRA + dy*cosRA;
                dyBody = -dx*cosRA*sinDec - dy*sinRA*sinDec + dz*cosDec;
                dzBody =  dx*cosRA*cosDec + dy*sinRA*cosDec + dz*sinDec;
                real fac1 = 3 * G * massj * forceParams.J2List[j] * radius *
                    radius / (2 * rRel5);
                real fac2 = 5 * dzBody * dzBody / rRel2 - 1;
                real axBody = fac1 * fac2 * dxBody;
                real ayBody = fac1 * fac2 * dyBody;
                real azBody = fac1 * (fac2 - 2) * dzBody;
                real axEquat = -axBody*sinRA - ayBody*cosRA*sinDec + azBody*cosRA*cosDec;
                real ayEquat =  axBody*cosRA - ayBody*sinRA*sinDec + azBody*sinRA*cosDec;
                real azEquat =  ayBody*cosDec + azBody*sinDec;
                ax += axEquat;
                ay += ayEquat;
                az += azEquat;
#ifdef PRINT_FORCES
                forceFile << std::setw(10) << forceParams.spiceIdList[j] << std::setw(25) << axEquat
                          << std::setw(25) << ayEquat << std::setw(25) << azEquat << std::endl;
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

void force_nongrav(const std::vector<real> &posAll,
                   const std::vector<real> &velAll, std::vector<real> &accInteg,
                   const ForceParameters &forceParams,
                   const IntegrationParameters &integParams,
                   const Constants &consts) {
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
    for (size_t i = 0; i < integParams.nInteg; i++) {
        x = posAll[3 * i];
        y = posAll[3 * i + 1];
        z = posAll[3 * i + 2];
        vx = velAll[3 * i];
        vy = velAll[3 * i + 1];
        vz = velAll[3 * i + 2];
        ax = 0.0;
        ay = 0.0;
        az = 0.0;
        for (size_t j = 0; j < integParams.nTotal; j++) {
            if (forceParams.spiceIdList[j] == 10 &&
                forceParams.isNongravList[i]) {  // j is Sun idx in spiceIdList
                a1 = forceParams.ngParamsList[i].a1;
                a2 = forceParams.ngParamsList[i].a2;
                a3 = forceParams.ngParamsList[i].a3;
                alpha = forceParams.ngParamsList[i].alpha;
                k = forceParams.ngParamsList[i].k;
                m = forceParams.ngParamsList[i].m;
                n = forceParams.ngParamsList[i].n;
                r0 = forceParams.ngParamsList[i].r0_au * 1.495978707e11 /
                    consts.du2m;
                dx = x - posAll[3 * j];
                dy = y - posAll[3 * j + 1];
                dz = z - posAll[3 * j + 2];
                dvx = vx - velAll[3 * j];
                dvy = vy - velAll[3 * j + 1];
                dvz = vz - velAll[3 * j + 2];
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
                forceFile << std::setw(10) << forceParams.spiceIdList[j] << std::setw(25)
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

void force_thruster(const std::vector<real> &velAll,
                    std::vector<real> &accInteg,
                    const ForceParameters &forceParams,
                    const IntegrationParameters &integParams,
                    const Constants &consts) {
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
    real acc_thruster = 1.0e7L / consts.du2m;  // m/day^2 -> au/day^2
    for (size_t i = 0; i < integParams.nInteg; i++) {
        vx = velAll[3 * i];
        vy = velAll[3 * i + 1];
        vz = velAll[3 * i + 2];
        ax = 0.0;
        ay = 0.0;
        az = 0.0;
        if (forceParams.isThrustingList[i]) {
            vunit({vx, vy, vz}, vHat);
            ax += acc_thruster * vHat[0];
            ay += acc_thruster * vHat[1];
            az += acc_thruster * vHat[2];
#ifdef PRINT_FORCES
            forceFile << "THRUSTER " << acc_thruster * vHat[0] << std::setw(25)
                      << acc_thruster * vHat[1] << std::setw(25) << acc_thruster * vHat[2]
                      << std::endl;
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
