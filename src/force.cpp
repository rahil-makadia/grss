#include "force.h"

void get_state_der(const real t, const std::vector<real> xInteg, std::vector<real> &xDotInteg, const ForceParameters forceParams, const IntegrationParameters integParams, const Constants consts){
    std::vector<real> posAll;
    std::vector<real> velAll;
    for (size_t i=0; i<integParams.nInteg; i++){
        for (size_t j=0; j<3; j++){
            posAll.push_back(xInteg[6*i+j]);
            velAll.push_back(xInteg[6*i+j+3]);
        }
    }
    for (size_t i=integParams.nInteg; i<integParams.nSpice; i++){
        double xSpice_i[6];
        double lt;
        get_spice_state_lt(forceParams.spiceIdList[i], t, consts, xSpice_i, lt);
        for (size_t j=0; j<3; j++){
            posAll.push_back(xSpice_i[j]);
            velAll.push_back(xSpice_i[j+3]);
        }
    }
    // print posAll and velAll
    std::cout << "posAll: " << std::endl;
    for (size_t i=0; i<posAll.size(); i+=3){
        std::cout << posAll[i] << " " << posAll[i+1] << " " << posAll[i+2] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "velAll: " << std::endl;
    for (size_t i=0; i<velAll.size(); i+=3){
        std::cout << velAll[i] << " " << velAll[i+1] << " " << velAll[i+2] << std::endl;
    }
    std::cout << std::endl;
    std::cout<< "xDotInteg0:" << std::endl;
    for (size_t i=0; i<xDotInteg.size(); i+=6){
        std::cout << xDotInteg[i] << " " << xDotInteg[i+1] << " " << xDotInteg[i+2] << " " << xDotInteg[i+3] << " " << xDotInteg[i+4] << " " << xDotInteg[i+5] << std::endl;
    }
    for (size_t i=0; i<integParams.nInteg; i++){
        xDotInteg[6*i] = velAll[3*i];
        xDotInteg[6*i+1] = velAll[3*i+1];
        xDotInteg[6*i+2] = velAll[3*i+2];
    }
    std::cout<< "xDotInteg1:" << std::endl;
    for (size_t i=0; i<xDotInteg.size(); i+=6){
        std::cout << xDotInteg[i] << " " << xDotInteg[i+1] << " " << xDotInteg[i+2] << " " << xDotInteg[i+3] << " " << xDotInteg[i+4] << " " << xDotInteg[i+5] << std::endl;
    }

    force_newton(t, posAll, velAll, xDotInteg, forceParams, integParams, consts);
    std::cout<< "xDotInteg2:" << std::endl;
    for (size_t i=0; i<xDotInteg.size(); i+=6){
        std::cout << xDotInteg[i] << " " << xDotInteg[i+1] << " " << xDotInteg[i+2] << " " << xDotInteg[i+3] << " " << xDotInteg[i+4] << " " << xDotInteg[i+5] << std::endl;
    }
    force_ppn(t, posAll, velAll, xDotInteg, forceParams, integParams, consts);
    std::cout<< "xDotInteg3:" << std::endl;
    for (size_t i=0; i<xDotInteg.size(); i+=6){
        std::cout << xDotInteg[i] << " " << xDotInteg[i+1] << " " << xDotInteg[i+2] << " " << xDotInteg[i+3] << " " << xDotInteg[i+4] << " " << xDotInteg[i+5] << std::endl;
    }
    force_J2(t, posAll, velAll, xDotInteg, forceParams, integParams, consts);
    std::cout<< "xDotInteg4:" << std::endl;
    for (size_t i=0; i<xDotInteg.size(); i+=6){
        std::cout << xDotInteg[i] << " " << xDotInteg[i+1] << " " << xDotInteg[i+2] << " " << xDotInteg[i+3] << " " << xDotInteg[i+4] << " " << xDotInteg[i+5] << std::endl;
    }
    force_nongrav(t, posAll, velAll, xDotInteg, forceParams, integParams, consts);
    std::cout<< "xDotInteg5:" << std::endl;
    for (size_t i=0; i<xDotInteg.size(); i+=6){
        std::cout << xDotInteg[i] << " " << xDotInteg[i+1] << " " << xDotInteg[i+2] << " " << xDotInteg[i+3] << " " << xDotInteg[i+4] << " " << xDotInteg[i+5] << std::endl;
    }
};

void force_newton(const real t, const std::vector<real> posAll, const std::vector<real> velAll, std::vector<real> &xDotInteg, const ForceParameters forceParams, const IntegrationParameters integParams, const Constants consts){
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
        for (size_t j=0; j<integParams.nTotal; j++){
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

void force_ppn(const real t, const std::vector<real> posAll, const std::vector<real> velAll, std::vector<real> &xDotInteg, const ForceParameters forceParams, const IntegrationParameters integParams, const Constants consts){
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
                // 1st order PPN correction
                ax -= G*massJ/(rRel3*c2)*((4*G*massJ/rRel - dVelDotVel)*dx + 4*dPosDotVel*dvx);
                ay -= G*massJ/(rRel3*c2)*((4*G*massJ/rRel - dVelDotVel)*dy + 4*dPosDotVel*dvy);
                az -= G*massJ/(rRel3*c2)*((4*G*massJ/rRel - dVelDotVel)*dz + 4*dPosDotVel*dvz);
            }
        }
        xDotInteg[6*i+3] += ax;
        xDotInteg[6*i+4] += ay;
        xDotInteg[6*i+5] += az;
    }
};

void force_J2(const real t, const std::vector<real> posAll, const std::vector<real> velAll, std::vector<real> &xDotInteg, const ForceParameters forceParams, const IntegrationParameters integParams, const Constants consts){
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
                rot_mat_x(forceParams.obliquityList[j], R2); // ecliptic to body equator
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

void force_nongrav(const real t, const std::vector<real> posAll, const std::vector<real> velAll, std::vector<real> &xDotInteg, const ForceParameters forceParams, const IntegrationParameters integParams, const Constants consts){
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
            if (forceParams.spiceIdList[j]==10 && forceParams.isNongravList[i]){ // j is Sun index in spiceIdList
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
                // real t1, t2, t3;
                // std::cout<< "eRHat: " << eRHat[0] << " " << eRHat[1] << " " << eRHat[2] << std::endl;
                // vnorm(eRHat, t1);
                // std::cout<< "eRHat mag: " << t1 << std::endl;
                // std::cout<< "eTHat: " << eTHat[0] << " " << eTHat[1] << " " << eTHat[2] << std::endl;
                // vnorm(eTHat, t2);
                // std::cout<< "eTHat mag: " << t2 << std::endl;
                // std::cout<< "eNHat: " << eNHat[0] << " " << eNHat[1] << " " << eNHat[2] << std::endl;
                // vnorm(eNHat, t3);
                // std::cout<< "eNHat mag: " << t3 << std::endl;
                ax = g*(a1*eRHat[0] + a2*eTHat[0] + a3*eNHat[0]);
                ay = g*(a1*eRHat[1] + a2*eTHat[1] + a3*eNHat[1]);
                az = g*(a1*eRHat[2] + a2*eTHat[2] + a3*eNHat[2]);
            }
        }
        xDotInteg[6*i+3] += ax;
        xDotInteg[6*i+4] += ay;
        xDotInteg[6*i+5] += az;
    }
};
