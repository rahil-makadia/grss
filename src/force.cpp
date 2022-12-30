#include "force.h"

std::vector<real> get_state_der(const real &t, const std::vector<real> &xInteg, const ForceParameters &forceParams, const IntegrationParameters &integParams, const Constants &consts){
    std::vector<real> xDotInteg(6*integParams.nInteg, 0.0);
    std::vector<real> posAll;
    std::vector<real> velAll;
    for (size_t i=0; i<integParams.nInteg; i++){
        for (size_t j=0; j<3; j++){
            posAll.push_back(xInteg[6*i+j]);
            velAll.push_back(xInteg[6*i+j+3]);
        }
    }
    for (size_t i=integParams.nInteg; i<integParams.nTotal; i++){
        double xSpice_i[6];
        double lt;
        get_spice_state_lt(forceParams.spiceIdList[i], t, consts, xSpice_i, lt);
        for (size_t j=0; j<3; j++){
            posAll.push_back(xSpice_i[j]);
            velAll.push_back(xSpice_i[j+3]);
        }
    }
    for (size_t i=0; i<integParams.nInteg; i++){
        xDotInteg[6*i] = velAll[3*i];
        xDotInteg[6*i+1] = velAll[3*i+1];
        xDotInteg[6*i+2] = velAll[3*i+2];
    }

    force_newton(posAll, xDotInteg, forceParams, integParams, consts);
    force_ppn(posAll, velAll, xDotInteg, forceParams, integParams, consts);
    force_J2(posAll, xDotInteg, forceParams, integParams, consts);
    force_nongrav(posAll, velAll, xDotInteg, forceParams, integParams, consts);

    std::vector<real> accInteg(3*integParams.nInteg, 0.0);
    for (size_t i=0; i<integParams.nInteg; i++){
        accInteg[3*i] = xDotInteg[6*i+3];
        accInteg[3*i+1] = xDotInteg[6*i+4];
        accInteg[3*i+2] = xDotInteg[6*i+5];
    }
    return accInteg;
};

void force_newton(const std::vector<real> &posAll, std::vector<real> &xDotInteg, const ForceParameters &forceParams, const IntegrationParameters &integParams, const Constants &consts){
    real G = consts.G;
    real x, y, z;
    real dx, dy, dz;
    real rRel, rRel3;
    real ax, ay, az;
    real massJ;
    for (size_t i=0; i<integParams.nInteg; i++){
        x = posAll[3*i];
        y = posAll[3*i+1];
        z = posAll[3*i+2];
        ax = 0.0;
        ay = 0.0;
        az = 0.0;
        for (size_t j=integParams.nInteg; j<integParams.nTotal; j++){
            if (i != j){
                massJ = forceParams.masses[j];
                dx = x - posAll[3*j];
                dy = y - posAll[3*j+1];
                dz = z - posAll[3*j+2];
                rRel = sqrt(dx*dx + dy*dy + dz*dz);
                rRel3 = rRel*rRel*rRel;
                ax -= G*massJ*dx/rRel3;
                ay -= G*massJ*dy/rRel3;
                az -= G*massJ*dz/rRel3;
            }
        }
        xDotInteg[6*i+3] += ax;
        xDotInteg[6*i+4] += ay;
        xDotInteg[6*i+5] += az;
    }
};

void force_ppn(const std::vector<real> &posAll, const std::vector<real> &velAll, std::vector<real> &xDotInteg, const ForceParameters &forceParams, const IntegrationParameters &integParams, const Constants &consts){
    real G = consts.G;
    real c = consts.clight;
    real c2 = c*c;
    real x, y, z;
    real dx, dy, dz;
    real rRel, rRel3;
    real vx, vy, vz;
    real dvx, dvy, dvz;
    real ax, ay, az;
    real massJ;
    real dPosDotVel, dVelDotVel;
    real beta = 1.0L;
    real gamma = 1.0L;
    for (size_t i=0; i<integParams.nInteg; i++){
        x = posAll[3*i];
        y = posAll[3*i+1];
        z = posAll[3*i+2];
        vx = velAll[3*i];
        vy = velAll[3*i+1];
        vz = velAll[3*i+2];
        ax = 0.0;
        ay = 0.0;
        az = 0.0;
        for (size_t j=0; j<integParams.nTotal; j++){
            if (i != j && forceParams.isPPNList[j]){
                massJ = forceParams.masses[j];
                dx = x - posAll[3*j];
                dy = y - posAll[3*j+1];
                dz = z - posAll[3*j+2];
                rRel = sqrt(dx*dx + dy*dy + dz*dz);
                rRel3 = rRel*rRel*rRel;
                dvx = vx - velAll[3*j];
                dvy = vy - velAll[3*j+1];
                dvz = vz - velAll[3*j+2];
                dPosDotVel = dx*dvx + dy*dvy + dz*dvz;
                dVelDotVel = dvx*dvx + dvy*dvy + dvz*dvz;
                // 1st order PPN approximation, equation 4-61 from Moyer (2003), https://descanso.jpl.nasa.gov/monograph/series2/Descanso2_all.pdf
                ax += G*massJ/(rRel3*c2)*((2*(beta+gamma)*G*massJ/rRel - dVelDotVel)*dx + 2*(1+gamma)*dPosDotVel*dvx);
                ay += G*massJ/(rRel3*c2)*((2*(beta+gamma)*G*massJ/rRel - dVelDotVel)*dy + 2*(1+gamma)*dPosDotVel*dvy);
                az += G*massJ/(rRel3*c2)*((2*(beta+gamma)*G*massJ/rRel - dVelDotVel)*dz + 2*(1+gamma)*dPosDotVel*dvz);
            }
        }
        xDotInteg[6*i+3] += ax;
        xDotInteg[6*i+4] += ay;
        xDotInteg[6*i+5] += az;
    }
};

void force_J2(const std::vector<real> &posAll, std::vector<real> &xDotInteg, const ForceParameters &forceParams, const IntegrationParameters &integParams, const Constants &consts){
    real G = consts.G;
    real x, y, z;
    real dx, dy, dz;
    real rRel, rRel2, rRel5;
    real radius;
    real ax, ay, az;
    real massJ;
    std::vector< std::vector<real> > R1(3, std::vector<real>(3));;
    std::vector< std::vector<real> > R2(3, std::vector<real>(3));;
    std::vector< std::vector<real> > R(3, std::vector<real>(3));;
    std::vector<real> dPosBody(3);
    for (size_t i=0; i<integParams.nInteg; i++){
        x = posAll[3*i];
        y = posAll[3*i+1];
        z = posAll[3*i+2];
        ax = 0.0;
        ay = 0.0;
        az = 0.0;
        for (size_t j=0; j<integParams.nTotal; j++){
            if (i != j && forceParams.isJ2List[j]){
                massJ = forceParams.masses[j];
                dx = x - posAll[3*j];
                dy = y - posAll[3*j+1];
                dz = z - posAll[3*j+2];
                rRel = sqrt(dx*dx + dy*dy + dz*dz);
                rRel2 = rRel*rRel;
                rRel5 = rRel2*rRel2*rRel;
                radius = forceParams.radii[j];
                rot_mat_x(EARTH_OBLIQUITY, R1); // equatorial to ecliptic
                rot_mat_x(-forceParams.obliquityList[j], R2); // ecliptic to body equator
                mat_mat_mul(R1, R2, R); // equatorial to body equator
                mat_vec_mul(R, {dx, dy, dz}, dPosBody);
                dx = dPosBody[0];
                dy = dPosBody[1];
                dz = dPosBody[2];
                ax -= 3*G*massJ*forceParams.J2List[j]*radius*radius*(1-5*dz*dz/rRel2)*dx/(2*rRel5);
                ay -= 3*G*massJ*forceParams.J2List[j]*radius*radius*(1-5*dz*dz/rRel2)*dy/(2*rRel5);
                az -= 3*G*massJ*forceParams.J2List[j]*radius*radius*(3-5*dz*dz/rRel2)*dz/(2*rRel5);
            }
        }
        xDotInteg[6*i+3] += ax;
        xDotInteg[6*i+4] += ay;
        xDotInteg[6*i+5] += az;
    }
};

void force_nongrav(const std::vector<real> &posAll, const std::vector<real> &velAll, std::vector<real> &xDotInteg, const ForceParameters &forceParams, const IntegrationParameters &integParams, const Constants &consts){
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
    for (size_t i=0; i<integParams.nInteg; i++){
        x = posAll[3*i];
        y = posAll[3*i+1];
        z = posAll[3*i+2];
        vx = velAll[3*i];
        vy = velAll[3*i+1];
        vz = velAll[3*i+2];
        ax = 0.0;
        ay = 0.0;
        az = 0.0;
        for (size_t j=0; j<integParams.nTotal; j++){
            if (forceParams.spiceIdList[j]==10 && forceParams.isNongravList[i]){ // j is Sun index (value 10) in spiceIdList
                // std::cout<< "Sun found at j = " << j << std::endl;
                a1 = forceParams.ngParamsList[i].a1;
                a2 = forceParams.ngParamsList[i].a2;
                a3 = forceParams.ngParamsList[i].a3;
                alpha = forceParams.ngParamsList[i].alpha;
                k = forceParams.ngParamsList[i].k;
                m = forceParams.ngParamsList[i].m;
                n = forceParams.ngParamsList[i].n;
                r0 = forceParams.ngParamsList[i].r0_au*1.495978707e11/consts.du2m;
                dx = x - posAll[3*j];
                dy = y - posAll[3*j+1];
                dz = z - posAll[3*j+2];
                dvx = vx - velAll[3*j];
                dvy = vy - velAll[3*j+1];
                dvz = vz - velAll[3*j+2]; 
                rRel = sqrt(dx*dx + dy*dy + dz*dz);
                g = pow(alpha*(rRel/r0), -m) * pow(1+pow(rRel/r0, n), -k);
                vunit({dx, dy, dz}, eRHat);
                vcross({dx, dy, dz}, {dvx, dvy, dvz}, hRelVec);
                vunit(hRelVec, eNHat);
                vcross(eNHat, eRHat, eTHat);
                ax += g*(a1*eRHat[0] + a2*eTHat[0] + a3*eNHat[0]);
                ay += g*(a1*eRHat[1] + a2*eTHat[1] + a3*eNHat[1]);
                az += g*(a1*eRHat[2] + a2*eTHat[2] + a3*eNHat[2]);
            }
        }
        xDotInteg[6*i+3] += ax;
        xDotInteg[6*i+4] += ay;
        xDotInteg[6*i+5] += az;
    }
};