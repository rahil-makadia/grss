#include "stm.h"

void bcd_and_dot(const std::vector<real> &stm, real *B, real *Bdot, real *C,
                 real *Cdot, real *D, real *Ddot) {
    B[0] = stm[0];
    B[1] = stm[1];
    B[2] = stm[2];
    C[0] = stm[3];
    C[1] = stm[4];
    C[2] = stm[5];
    B[3] = stm[6];
    B[4] = stm[7];
    B[5] = stm[8];
    C[3] = stm[9];
    C[4] = stm[10];
    C[5] = stm[11];
    B[6] = stm[12];
    B[7] = stm[13];
    B[8] = stm[14];
    C[6] = stm[15];
    C[7] = stm[16];
    C[8] = stm[17];

    Bdot[0] = stm[18];
    Bdot[1] = stm[19];
    Bdot[2] = stm[20];
    Cdot[0] = stm[21];
    Cdot[1] = stm[22];
    Cdot[2] = stm[23];
    Bdot[3] = stm[24];
    Bdot[4] = stm[25];
    Bdot[5] = stm[26];
    Cdot[3] = stm[27];
    Cdot[4] = stm[28];
    Cdot[5] = stm[29];
    Bdot[6] = stm[30];
    Bdot[7] = stm[31];
    Bdot[8] = stm[32];
    Cdot[6] = stm[33];
    Cdot[7] = stm[34];
    Cdot[8] = stm[35];

    const size_t numParams = (stm.size() - 36) / 6;
    size_t startd = 0;
    size_t startStm = 36;
    for (size_t param = 0; param < numParams; param++) {
        D[startd] = stm[startStm];
        D[startd + 1] = stm[startStm + 1];
        D[startd + 2] = stm[startStm + 2];
        Ddot[startd] = stm[startStm + 3];
        Ddot[startd + 1] = stm[startStm + 4];
        Ddot[startd + 2] = stm[startStm + 5];
        startd += 3;
        startStm += 6;
    }
}

void bcd_2dot(const real *B, const real *Bdot, const real *C, const real *Cdot,
              const real *D, const real *Ddot, const real *dfdpos,
              const real *dfdvel, const real *dfdpar, size_t numParams,
              size_t stmStarti, std::vector<real> &accInteg) {
    real *termB1 = new real[9];
    real *termB2 = new real[9];
    real *B2dot = new real[9];

    mat3_mat3_mul(dfdpos, B, termB1);
    mat3_mat3_mul(dfdvel, Bdot, termB2);
    mat3_mat3_add(termB1, termB2, B2dot);

    real *termC1 = new real[9];
    real *termC2 = new real[9];
    real *C2dot = new real[9];

    mat3_mat3_mul(dfdpos, C, termC1);
    mat3_mat3_mul(dfdvel, Cdot, termC2);
    mat3_mat3_add(termC1, termC2, C2dot);

    accInteg[stmStarti + 0] += B2dot[0];
    accInteg[stmStarti + 1] += B2dot[1];
    accInteg[stmStarti + 2] += B2dot[2];
    accInteg[stmStarti + 3] += C2dot[0];
    accInteg[stmStarti + 4] += C2dot[1];
    accInteg[stmStarti + 5] += C2dot[2];
    accInteg[stmStarti + 6] += B2dot[3];
    accInteg[stmStarti + 7] += B2dot[4];
    accInteg[stmStarti + 8] += B2dot[5];
    accInteg[stmStarti + 9] += C2dot[3];
    accInteg[stmStarti + 10] += C2dot[4];
    accInteg[stmStarti + 11] += C2dot[5];
    accInteg[stmStarti + 12] += B2dot[6];
    accInteg[stmStarti + 13] += B2dot[7];
    accInteg[stmStarti + 14] += B2dot[8];
    accInteg[stmStarti + 15] += C2dot[6];
    accInteg[stmStarti + 16] += C2dot[7];
    accInteg[stmStarti + 17] += C2dot[8];

    if (numParams > 0) {
        real *D2dot = new real[3 * numParams];
        for (size_t param = 0; param < numParams; param++) {
            D2dot[3 * param + 0] = dfdpos[0] * D[3 * param + 0] +
                dfdpos[1] * D[3 * param + 1] + dfdpos[2] * D[3 * param + 2] +
                dfdvel[0] * Ddot[3 * param + 0] +
                dfdvel[1] * Ddot[3 * param + 1] +
                dfdvel[2] * Ddot[3 * param + 2] + dfdpar[3 * param + 0];
            D2dot[3 * param + 1] = dfdpos[3] * D[3 * param + 0] +
                dfdpos[4] * D[3 * param + 1] + dfdpos[5] * D[3 * param + 2] +
                dfdvel[3] * Ddot[3 * param + 0] +
                dfdvel[4] * Ddot[3 * param + 1] +
                dfdvel[5] * Ddot[3 * param + 2] + dfdpar[3 * param + 1];
            D2dot[3 * param + 2] = dfdpos[6] * D[3 * param + 0] +
                dfdpos[7] * D[3 * param + 1] + dfdpos[8] * D[3 * param + 2] +
                dfdvel[6] * Ddot[3 * param + 0] +
                dfdvel[7] * Ddot[3 * param + 1] +
                dfdvel[8] * Ddot[3 * param + 2] + dfdpar[3 * param + 2];
            accInteg[stmStarti + 18 + 3 * param + 0] += D2dot[3 * param + 0];
            accInteg[stmStarti + 18 + 3 * param + 1] += D2dot[3 * param + 1];
            accInteg[stmStarti + 18 + 3 * param + 2] += D2dot[3 * param + 2];
        }
    }
}

void stm_newton(const IntegBody &bodyi, const real &gm, const real &dx,
                const real &dy, const real &dz, const size_t &stmStarti,
                std::vector<real> &accInteg) {
    real *dfdpos = new real[9];
    real *dfdvel = new real[9];
    memset(dfdvel, 0, 9 * sizeof(real));
    const size_t numParams = (bodyi.stm.size() - 36) / 6;
    real *dfdpar = new real[3 * numParams];
    memset(dfdpar, 0, 3 * numParams * sizeof(real));

    real *B = new real[9];
    real *Bdot = new real[9];
    real *C = new real[9];
    real *Cdot = new real[9];
    real *D = new real[3 * numParams];
    real *Ddot = new real[3 * numParams];
    bcd_and_dot(bodyi.stm, B, Bdot, C, Cdot, D, Ddot);

    const real r = sqrt(dx * dx + dy * dy + dz * dz);
    const real r3 = r * r * r;
    const real r5 = r3 * r * r;
    dfdpos[0] = gm * (3 * dx * dx / r5 - 1 / r3);
    dfdpos[1] = gm * (3 * dx * dy / r5);
    dfdpos[2] = gm * (3 * dx * dz / r5);
    dfdpos[3] = gm * (3 * dx * dy / r5);
    dfdpos[4] = gm * (3 * dy * dy / r5 - 1 / r3);
    dfdpos[5] = gm * (3 * dy * dz / r5);
    dfdpos[6] = gm * (3 * dx * dz / r5);
    dfdpos[7] = gm * (3 * dy * dz / r5);
    dfdpos[8] = gm * (3 * dz * dz / r5 - 1 / r3);

    bcd_2dot(B, Bdot, C, Cdot, D, Ddot, dfdpos, dfdvel, dfdpar, numParams,
             stmStarti, accInteg);

    delete[] dfdpos;
    delete[] dfdvel;
    delete[] dfdpar;
    delete[] B;
    delete[] Bdot;
    delete[] C;
    delete[] Cdot;
    delete[] D;
    delete[] Ddot;
}

void stm_ppn_simple(const IntegBody &bodyi, const real &gm, const real &c,
                    const real &beta, const real &gamma, const real &dx,
                    const real &dy, const real &dz, const real &dvx,
                    const real &dvy, const real &dvz, const size_t &stmStarti,
                    std::vector<real> &accInteg) {
    real *dfdpos = new real[9];
    real *dfdvel = new real[9];
    const size_t numParams = (bodyi.stm.size() - 36) / 6;
    real *dfdpar = new real[3 * numParams];
    memset(dfdpar, 0, 3 * numParams * sizeof(real));

    real *B = new real[9];
    real *Bdot = new real[9];
    real *C = new real[9];
    real *Cdot = new real[9];
    real *D = new real[3 * numParams];
    real *Ddot = new real[3 * numParams];
    bcd_and_dot(bodyi.stm, B, Bdot, C, Cdot, D, Ddot);

    const real c2 = c * c;
    const real r = sqrt(dx * dx + dy * dy + dz * dz);
    const real r3 = r * r * r;

    const real rDotV = dx * dvx + dy * dvy + dz * dvz;
    const real v2 = dvx * dvx + dvy * dvy + dvz * dvz;

    const real fac1 = gm / (c2 * r3);
    const real dfac1dr = -3 * gm / (c2 * r3 * r);
    const real fac2 =
        (2 * (beta + gamma) * gm / r - gamma * v2);
    const real dfac2dr = -2 * (beta + gamma) * gm / (r * r);
    const real fac3 = 2 * (1 + gamma) * rDotV;

    dfdpos[0] = dfac1dr*dx/r * (fac2*dx + fac3*dvx) + fac1*(fac2 + dfac2dr*dx*dx/r + 2*(1+gamma)*dvx*dvx);
    dfdpos[1] = dfac1dr*dy/r * (fac2*dx + fac3*dvx) + fac1*(dfac2dr*dx*dy/r + 2*(1+gamma)*dvx*dvy);
    dfdpos[2] = dfac1dr*dz/r * (fac2*dx + fac3*dvx) + fac1*(dfac2dr*dx*dz/r + 2*(1+gamma)*dvx*dvz);
    dfdpos[3] = dfac1dr*dx/r * (fac2*dy + fac3*dvy) + fac1*(dfac2dr*dy*dx/r + 2*(1+gamma)*dvy*dvx);
    dfdpos[4] = dfac1dr*dy/r * (fac2*dy + fac3*dvy) + fac1*(fac2 + dfac2dr*dy*dy/r + 2*(1+gamma)*dvy*dvy);
    dfdpos[5] = dfac1dr*dz/r * (fac2*dy + fac3*dvy) + fac1*(dfac2dr*dy*dz/r + 2*(1+gamma)*dvy*dvz);
    dfdpos[6] = dfac1dr*dx/r * (fac2*dz + fac3*dvz) + fac1*(dfac2dr*dz*dx/r + 2*(1+gamma)*dvz*dvx);
    dfdpos[7] = dfac1dr*dy/r * (fac2*dz + fac3*dvz) + fac1*(dfac2dr*dz*dy/r + 2*(1+gamma)*dvz*dvy);
    dfdpos[8] = dfac1dr*dz/r * (fac2*dz + fac3*dvz) + fac1*(fac2 + dfac2dr*dz*dz/r + 2*(1+gamma)*dvz*dvz);

    dfdvel[0] = fac1*(-2*gamma*dvx*dx + fac3 + 2*(1+gamma)*dx*dvx);
    dfdvel[1] = fac1*(-2*gamma*dvy*dx + 2*(1+gamma)*dy*dvx);
    dfdvel[2] = fac1*(-2*gamma*dvz*dx + 2*(1+gamma)*dz*dvx);
    dfdvel[3] = fac1*(-2*gamma*dvx*dy + 2*(1+gamma)*dx*dvy);
    dfdvel[4] = fac1*(-2*gamma*dvy*dy + fac3 + 2*(1+gamma)*dy*dvy);
    dfdvel[5] = fac1*(-2*gamma*dvz*dy + 2*(1+gamma)*dz*dvy);
    dfdvel[6] = fac1*(-2*gamma*dvx*dz + 2*(1+gamma)*dx*dvz);
    dfdvel[7] = fac1*(-2*gamma*dvy*dz + 2*(1+gamma)*dy*dvz);
    dfdvel[8] = fac1*(-2*gamma*dvz*dz + fac3 + 2*(1+gamma)*dz*dvz);

    bcd_2dot(B, Bdot, C, Cdot, D, Ddot, dfdpos, dfdvel, dfdpar, numParams,
             stmStarti, accInteg);

    delete[] dfdpos;
    delete[] dfdvel;
    delete[] dfdpar;
    delete[] B;
    delete[] Bdot;
    delete[] C;
    delete[] Cdot;
    delete[] D;
    delete[] Ddot;
}

void stm_J2(const IntegBody &bodyi, const real &gm, const real &J2,
            const real &dxBody, const real &dyBody, const real &dzBody,
            const real &radius, const real &sinRA, const real &cosRA,
            const real &sinDec, const real &cosDec,
            const real &smoothing_threshold, const size_t &stmStarti,
            std::vector<real> &accInteg) {
    real *dfdpos = new real[9];
    real *dfdvel = new real[9];
    memset(dfdvel, 0, 9 * sizeof(real));
    const size_t numParams = (bodyi.stm.size() - 36) / 6;
    real *dfdpar = new real[3 * numParams];
    memset(dfdpar, 0, 3 * numParams * sizeof(real));

    real *B = new real[9];
    real *Bdot = new real[9];
    real *C = new real[9];
    real *Cdot = new real[9];
    real *D = new real[3 * numParams];
    real *Ddot = new real[3 * numParams];
    bcd_and_dot(bodyi.stm, B, Bdot, C, Cdot, D, Ddot);

    const real r2 = dxBody * dxBody + dyBody * dyBody + dzBody * dzBody;
    const real r = sqrt(r2);
    const real r4 = r2 * r2;
    const real r5 = r2 * r2 * r;
    const real r7 = r5 * r2;

    const real fac1 = 3 * gm * J2 * radius * radius / (2 * r5);
    const real dfac1Fac = -(15*gm*J2*radius*radius)/(2*r7);
    const real dfac1dxBody = dfac1Fac*dxBody;
    const real dfac1dyBody = dfac1Fac*dyBody;
    const real dfac1dzBody = dfac1Fac*dzBody;

    const real fac2 = 5 * dzBody * dzBody / r2 - 1;
    const real dfac2dxBody = -10*dzBody*dzBody*dxBody/r4;
    const real dfac2dyBody = -10*dzBody*dzBody*dyBody/r4;
    const real dfac2dzBody = 10*dzBody/r2 - 10*dzBody*dzBody*dzBody/r4;

    real *dfBodydposBody = new real[9];
    dfBodydposBody[0] = dfac1dxBody*fac2*dxBody + fac1*(dfac2dxBody*dxBody + fac2);
    dfBodydposBody[1] = dfac1dyBody*fac2*dxBody + fac1*dfac2dyBody*dxBody;
    dfBodydposBody[2] = dfac1dzBody*fac2*dxBody + fac1*dfac2dzBody*dxBody;
    dfBodydposBody[3] = dfac1dxBody*fac2*dyBody + fac1*dfac2dxBody*dyBody;
    dfBodydposBody[4] = dfac1dyBody*fac2*dyBody + fac1*(dfac2dyBody*dyBody + fac2);
    dfBodydposBody[5] = dfac1dzBody*fac2*dyBody + fac1*dfac2dzBody*dyBody;
    dfBodydposBody[6] = dfac1dxBody*(fac2 - 2)*dzBody + fac1*dfac2dxBody*dzBody;
    dfBodydposBody[7] = dfac1dyBody*(fac2 - 2)*dzBody + fac1*dfac2dyBody*dzBody;
    dfBodydposBody[8] = dfac1dzBody*(fac2 - 2)*dzBody + fac1*(dfac2dzBody*dzBody + fac2 - 2);

    if (r <= radius+smoothing_threshold) {
        const real depth = radius+smoothing_threshold-r;
        real smoothing = cos(PI*depth/(2*smoothing_threshold));
        if (depth > smoothing_threshold){
            smoothing = 0.0;
        }
        if (smoothing != 0.0){
            const real dsmoothingdxBody =
                sin(PI * depth / (2 * smoothing_threshold)) * PI * dxBody /
                (2 * smoothing_threshold * r);
            const real dsmoothingdyBody =
                sin(PI * depth / (2 * smoothing_threshold)) * PI * dyBody /
                (2 * smoothing_threshold * r);
            const real dsmoothingdzBody =
                sin(PI * depth / (2 * smoothing_threshold)) * PI * dzBody /
                (2 * smoothing_threshold * r);
            dfBodydposBody[0] *= smoothing;
            dfBodydposBody[0] += fac1*fac2*dxBody*dsmoothingdxBody;
            dfBodydposBody[1] *= smoothing;
            dfBodydposBody[1] += fac1*fac2*dxBody*dsmoothingdyBody;
            dfBodydposBody[2] *= smoothing;
            dfBodydposBody[2] += fac1*fac2*dxBody*dsmoothingdzBody;
            dfBodydposBody[3] *= smoothing;
            dfBodydposBody[3] += fac1*fac2*dyBody*dsmoothingdxBody;
            dfBodydposBody[4] *= smoothing;
            dfBodydposBody[4] += fac1*fac2*dyBody*dsmoothingdyBody;
            dfBodydposBody[5] *= smoothing;
            dfBodydposBody[5] += fac1*fac2*dyBody*dsmoothingdzBody;
            dfBodydposBody[6] *= smoothing;
            dfBodydposBody[6] += fac1*(fac2-2)*dzBody*dsmoothingdxBody;
            dfBodydposBody[7] *= smoothing;
            dfBodydposBody[7] += fac1*(fac2-2)*dzBody*dsmoothingdyBody;
            dfBodydposBody[8] *= smoothing;
            dfBodydposBody[8] += fac1*(fac2-2)*dzBody*dsmoothingdzBody;
        }
    }

    real *dposBodydpos = new real[9];
    dposBodydpos[0] = -sinRA;
    dposBodydpos[1] = cosRA;
    dposBodydpos[2] = 0.0;
    dposBodydpos[3] = -cosRA*sinDec;
    dposBodydpos[4] = -sinRA*sinDec;
    dposBodydpos[5] = cosDec;
    dposBodydpos[6] = cosRA*cosDec;
    dposBodydpos[7] = sinRA*cosDec;
    dposBodydpos[8] = sinDec;

    real* dfdfBody = new real[9];
    dfdfBody[0] = -sinRA;
    dfdfBody[1] = -cosRA*sinDec;
    dfdfBody[2] = cosDec*cosRA;
    dfdfBody[3] = cosRA;
    dfdfBody[4] = -sinRA*sinDec;
    dfdfBody[5] = cosDec*sinRA;
    dfdfBody[6] = 0.0;
    dfdfBody[7] = cosDec;
    dfdfBody[8] = sinDec;

    real *dfBodydpos = new real[9];
    mat3_mat3_mul(dfBodydposBody, dposBodydpos, dfBodydpos);
    mat3_mat3_mul(dfdfBody, dfBodydpos, dfdpos);

    bcd_2dot(B, Bdot, C, Cdot, D, Ddot, dfdpos, dfdvel, dfdpar, numParams,
             stmStarti, accInteg);

    delete[] dfdpos;
    delete[] dfdvel;
    delete[] dfdpar;
    delete[] B;
    delete[] Bdot;
    delete[] C;
    delete[] Cdot;
    delete[] D;
    delete[] Ddot;
    delete[] dfBodydposBody;
    delete[] dposBodydpos;
    delete[] dfdfBody;
    delete[] dfBodydpos;
}

void stm_nongrav(const IntegBody &bodyi, const real &g,
                 const NongravParamaters &ngParams, const real &dx, const real &dy,
                 const real &dz, const real &dvx, const real &dvy, const real &dvz,
                 real *rVec, real *nVec, const size_t &stmStarti,
                 std::vector<real> &accInteg) {
    // maybe pull this common stuff out into a function
    // and return an STM struct in the force functions
    // that gets passed into the STM functions
    real *dfdpos = new real[9];
    real *dfdvel = new real[9];
    const size_t numParams = (bodyi.stm.size() - 36) / 6;
    real *dfdpar = new real[3 * numParams];

    real *B = new real[9];
    real *Bdot = new real[9];
    real *C = new real[9];
    real *Cdot = new real[9];
    real *D = new real[3 * numParams];
    real *Ddot = new real[3 * numParams];
    bcd_and_dot(bodyi.stm, B, Bdot, C, Cdot, D, Ddot);

    const real a1 = ngParams.a1;
    const real a2 = ngParams.a2;
    const real a3 = ngParams.a3;
    const real alpha = ngParams.alpha;
    const real k = ngParams.k;
    const real m = ngParams.m;
    const real n = ngParams.n;
    const real r0 = ngParams.r0_au;

    const real rNorm = sqrt(dx * dx + dy * dy + dz * dz);
    const real r3 = rNorm * rNorm * rNorm;
    const real v = sqrt(dvx * dvx + dvy * dvy + dvz * dvz);
    const real vDotV = v * v;
    const real rDotV = dx * dvx + dy * dvy + dz * dvz;

    real rHat[3] = {dx / rNorm, dy / rNorm, dz / rNorm};
    real tNorm;
    real *tVec = new real[3];
    vcross(nVec, rVec, tVec);
    const real vDotT = dvx * tVec[0] + dvy * tVec[1] + dvz * tVec[2];
    vnorm(tVec, 3, tNorm);
    real tHat[3] = {tVec[0] / tNorm, tVec[1] / tNorm, tVec[2] / tNorm};
    real nNorm;
    vnorm(nVec, 3, nNorm);
    real nHat[3] = {nVec[0] / nNorm, nVec[1] / nNorm, nVec[2] / nNorm};

    const real rxr3 = dx / r3;
    const real ryr3 = dy / r3;
    const real rzr3 = dz / r3;
    const real t3 = tNorm * tNorm * tNorm;
    const real txt3 = tVec[0] / t3;
    const real tyt3 = tVec[1] / t3;
    const real tzt3 = tVec[2] / t3;
    const real n3 = nNorm * nNorm * nNorm;
    const real nxn3 = nVec[0] / n3;
    const real nyn3 = nVec[1] / n3;
    const real nzn3 = nVec[2] / n3;

    const real dgdr =
        - alpha * m / r0 * pow(rNorm / r0, -m - 1) * pow(1 + pow(rNorm / r0, n), -k)
        - alpha * k * n / r0 * pow(rNorm / r0, -m) *
            pow(1 + pow(rNorm / r0, n), -k - 1) * pow(rNorm / r0, n - 1);

    const real dgdx = dgdr * rHat[0];
    const real dgdy = dgdr * rHat[1];
    const real dgdz = dgdr * rHat[2];

    dfdpos[0] = a1*(dgdx*rHat[0] - g*rxr3*dx + g/rNorm)
                + a2*(dgdx*tHat[0] - g*txt3*(2*dx*vDotT - tVec[0]*rDotV) + g*(dx*dvx - rDotV)/tNorm)
                + a3*(dgdx*nHat[0] - g*nxn3*(dx*vDotV - dvx*rDotV));
    dfdpos[1] = a1*(dgdy*rHat[0] - g*rxr3*dy)
                + a2*(dgdy*tHat[0] - g*txt3*(2*dy*vDotT - tVec[1]*rDotV) + g*(2*dy*dvx - dx*dvy)/tNorm)
                + a3*(dgdy*nHat[0] - g*nxn3*(dy*vDotV - dvy*rDotV) + g*dvz/nNorm);
    dfdpos[2] = a1*(dgdz*rHat[0] - g*rxr3*dz)
                + a2*(dgdz*tHat[0] - g*txt3*(2*dz*vDotT - tVec[2]*rDotV) + g*(2*dz*dvx - dx*dvz)/tNorm)
                + a3*(dgdz*nHat[0] - g*nxn3*(dz*vDotV - dvz*rDotV) - g*dvy/nNorm);
    dfdpos[3] = a1*(dgdx*rHat[1] - g*ryr3*dx)
                + a2*(dgdx*tHat[1] - g*tyt3*(2*dx*vDotT - tVec[0]*rDotV) + g*(2*dx*dvy - dy*dvx)/tNorm)
                + a3*(dgdx*nHat[1] - g*nyn3*(dx*vDotV - dvx*rDotV) - g*dvz/nNorm);
    dfdpos[4] = a1*(dgdy*rHat[1] - g*ryr3*dy + g/rNorm)
                + a2*(dgdy*tHat[1] - g*tyt3*(2*dy*vDotT - tVec[1]*rDotV) + g*(dy*dvy - rDotV)/tNorm)
                + a3*(dgdy*nHat[1] - g*nyn3*(dy*vDotV - dvy*rDotV));
    dfdpos[5] = a1*(dgdz*rHat[1] - g*ryr3*dz)
                + a2*(dgdz*tHat[1] - g*tyt3*(2*dz*vDotT - tVec[2]*rDotV) + g*(2*dz*dvy - dy*dvz)/tNorm)
                + a3*(dgdz*nHat[1] - g*nyn3*(dz*vDotV - dvz*rDotV) + g*dvx/nNorm);
    dfdpos[6] = a1*(dgdx*rHat[2] - g*rzr3*dx)
                + a2*(dgdx*tHat[2] - g*tzt3*(2*dx*vDotT - tVec[0]*rDotV) + g*(2*dx*dvz - dz*dvx)/tNorm)
                + a3*(dgdx*nHat[2] - g*nzn3*(dx*vDotV - dvx*rDotV) + g*dvy/nNorm);
    dfdpos[7] = a1*(dgdy*rHat[2] - g*rzr3*dy)
                + a2*(dgdy*tHat[2] - g*tzt3*(2*dy*vDotT - tVec[1]*rDotV) + g*(2*dy*dvz - dz*dvy)/tNorm)
                + a3*(dgdy*nHat[2] - g*nzn3*(dy*vDotV - dvy*rDotV) - g*dvx/nNorm);
    dfdpos[8] = a1*(dgdz*rHat[2] - g*rzr3*dz + g/rNorm)
                + a2*(dgdz*tHat[2] - g*tzt3*(2*dz*vDotT - tVec[2]*rDotV) + g*(dz*dvz - rDotV)/tNorm)
                + a3*(dgdz*nHat[2] - g*nzn3*(dz*vDotV - dvz*rDotV));
    
    const real r2 = rNorm * rNorm;
    dfdvel[0] = g*(a2*((dy*dy + dz*dz)/tNorm - txt3*r2*tVec[0]) - a3*(nxn3*(r2*dvx - dx*rDotV)));
    dfdvel[1] = g*(a2*(-dx*dy/tNorm - tyt3*r2*tVec[0]) - a3*(nxn3*(r2*dvy - dy*rDotV) + dz/nNorm));
    dfdvel[2] = g*(a2*(-dx*dz/tNorm - tzt3*r2*tVec[0]) - a3*(nxn3*(r2*dvz - dz*rDotV) - dy/nNorm));
    dfdvel[3] = g*(a2*(-dy*dx/tNorm - txt3*r2*tVec[1]) - a3*(nyn3*(r2*dvx - dx*rDotV) - dz/nNorm));
    dfdvel[4] = g*(a2*((dx*dx + dz*dz)/tNorm - tyt3*r2*tVec[1]) - a3*(nyn3*(r2*dvy - dy*rDotV)));
    dfdvel[5] = g*(a2*(-dy*dz/tNorm - tzt3*r2*tVec[1]) - a3*(nyn3*(r2*dvz - dz*rDotV) + dx/nNorm));
    dfdvel[6] = g*(a2*(-dz*dx/tNorm - txt3*r2*tVec[2]) - a3*(nzn3*(r2*dvx - dx*rDotV) + dy/nNorm));
    dfdvel[7] = g*(a2*(-dz*dy/tNorm - tyt3*r2*tVec[2]) - a3*(nzn3*(r2*dvy - dy*rDotV) - dx/nNorm));
    dfdvel[8] = g*(a2*((dx*dx + dy*dy)/tNorm - tzt3*r2*tVec[2]) - a3*(nzn3*(r2*dvz - dz*rDotV)));

    // Case 1, only one of a1, a2, or a3 is non-zero
    if (a1 != 0 && a2 == 0 && a3 == 0) {
        dfdpar[0] = g*rHat[0];
        dfdpar[1] = g*rHat[1];
        dfdpar[2] = g*rHat[2];
    } else if (a1 == 0 && a2 != 0 && a3 == 0) {
        dfdpar[0] = g*tHat[0];
        dfdpar[1] = g*tHat[1];
        dfdpar[2] = g*tHat[2];
    } else if (a1 == 0 && a2 == 0 && a3 != 0) {
        dfdpar[0] = g*nHat[0];
        dfdpar[1] = g*nHat[1];
        dfdpar[2] = g*nHat[2];
    }
    // Case 2, two of a1, a2, or a3 are non-zero
    if (a1 != 0 && a2 != 0 && a3 == 0) {
        dfdpar[0] = g*rHat[0];
        dfdpar[1] = g*rHat[1];
        dfdpar[2] = g*rHat[2];
        dfdpar[3] = g*tHat[0];
        dfdpar[4] = g*tHat[1];
        dfdpar[5] = g*tHat[2];
    } else if (a1 != 0 && a2 == 0 && a3 != 0) {
        dfdpar[0] = g*rHat[0];
        dfdpar[1] = g*rHat[1];
        dfdpar[2] = g*rHat[2];
        dfdpar[3] = g*nHat[0];
        dfdpar[4] = g*nHat[1];
        dfdpar[5] = g*nHat[2];
    } else if (a1 == 0 && a2 != 0 && a3 != 0) {
        dfdpar[0] = g*tHat[0];
        dfdpar[1] = g*tHat[1];
        dfdpar[2] = g*tHat[2];
        dfdpar[3] = g*nHat[0];
        dfdpar[4] = g*nHat[1];
        dfdpar[5] = g*nHat[2];
    }
    // Case 3, all of a1, a2, and a3 are non-zero
    if (a1 != 0 && a2 != 0 && a3 != 0) {
        dfdpar[0] = g*rHat[0];
        dfdpar[1] = g*rHat[1];
        dfdpar[2] = g*rHat[2];
        dfdpar[3] = g*tHat[0];
        dfdpar[4] = g*tHat[1];
        dfdpar[5] = g*tHat[2];
        dfdpar[6] = g*nHat[0];
        dfdpar[7] = g*nHat[1];
        dfdpar[8] = g*nHat[2];
    }

    bcd_2dot(B, Bdot, C, Cdot, D, Ddot, dfdpos, dfdvel, dfdpar, numParams,
             stmStarti, accInteg);

    delete[] dfdpos;
    delete[] dfdvel;
    delete[] dfdpar;
    delete[] B;
    delete[] Bdot;
    delete[] C;
    delete[] Cdot;
    delete[] D;
    delete[] Ddot;
    delete[] tVec;
}
