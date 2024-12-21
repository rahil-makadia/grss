/**
 * @file    stm.cpp
 * @brief   Source file for the STM functions.
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

#include "stm.h"

// declared in simulation.h
/**
 * @param[in] stm The flattened STM vector.
 * @return The reconstructed STM matrix.
 */
std::vector<std::vector<real>> reconstruct_stm(const std::vector<real> &stm) {
    const size_t numParams = (stm.size() - 36) / 6;
    std::vector<std::vector<real>> stmMatrix = std::vector<std::vector<real>>(
        6 + numParams, std::vector<real>(6 + numParams, 0.0));
    real *B = new real[9];
    real *Bdot = new real[9];
    real *C = new real[9];
    real *Cdot = new real[9];
    real *D = new real[3 * numParams];
    real *Ddot = new real[3 * numParams];
    bcd_and_dot(stm, B, Bdot, C, Cdot, D, Ddot);
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            stmMatrix[i][j] = B[3 * i + j];
            stmMatrix[i][j + 3] = C[3 * i + j];
            stmMatrix[i + 3][j] = Bdot[3 * i + j];
            stmMatrix[i + 3][j + 3] = Cdot[3 * i + j];
        }
    }
    for (size_t j = 6; j < 6 + numParams; j++) {
        for (size_t i = 0; i < 3; i++) {
            stmMatrix[i][j] = D[3 * (j - 6) + i];
            stmMatrix[i + 3][j] = Ddot[3 * (j - 6) + i];
        }
        stmMatrix[j][j] = 1.0;
    }
    return stmMatrix;
}

/**
 * @param[in] stm The full STM.
 * @param[out] B The B submatrix (dpos_dpos).
 * @param[out] Bdot The Bdot submatrix (dpos_dvel).
 * @param[out] C The C submatrix (dvel_dpos).
 * @param[out] Cdot The Cdot submatrix (dvel_dvel).
 * @param[out] D The D submatrix (dpos_dpar).
 * @param[out] Ddot The Ddot submatrix (dvel_dpar).
 */
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

/**
 * @param[in] stmParams Structure containing the STM submatrices and their derivatives.
 * @param[in] numParams Number of parameters.
 * @param[in] stmStarti Index of the first element of the STM in the second derivative vector.
 * @param[out] accInteg Second derivative vector.
 */
void bcd_2dot(STMParameters &stmParams, size_t numParams,
              size_t stmStarti, std::vector<real> &accInteg) {
    real *termB1 = new real[9];
    real *termB2 = new real[9];
    real *B2dot = new real[9];

    mat3_mat3_mul(stmParams.dfdpos, stmParams.B, termB1);
    mat3_mat3_mul(stmParams.dfdvel, stmParams.Bdot, termB2);
    mat3_mat3_add(termB1, termB2, B2dot);

    real *termC1 = new real[9];
    real *termC2 = new real[9];
    real *C2dot = new real[9];

    mat3_mat3_mul(stmParams.dfdpos, stmParams.C, termC1);
    mat3_mat3_mul(stmParams.dfdvel, stmParams.Cdot, termC2);
    mat3_mat3_add(termC1, termC2, C2dot);

    accInteg[stmStarti + 0] = B2dot[0];
    accInteg[stmStarti + 1] = B2dot[1];
    accInteg[stmStarti + 2] = B2dot[2];
    accInteg[stmStarti + 3] = C2dot[0];
    accInteg[stmStarti + 4] = C2dot[1];
    accInteg[stmStarti + 5] = C2dot[2];
    accInteg[stmStarti + 6] = B2dot[3];
    accInteg[stmStarti + 7] = B2dot[4];
    accInteg[stmStarti + 8] = B2dot[5];
    accInteg[stmStarti + 9] = C2dot[3];
    accInteg[stmStarti + 10] = C2dot[4];
    accInteg[stmStarti + 11] = C2dot[5];
    accInteg[stmStarti + 12] = B2dot[6];
    accInteg[stmStarti + 13] = B2dot[7];
    accInteg[stmStarti + 14] = B2dot[8];
    accInteg[stmStarti + 15] = C2dot[6];
    accInteg[stmStarti + 16] = C2dot[7];
    accInteg[stmStarti + 17] = C2dot[8];

    if (numParams > 0) {
        real *D2dot = new real[3 * numParams];
        for (size_t param = 0; param < numParams; param++) {
            D2dot[3 * param + 0] = stmParams.dfdpos[0] * stmParams.D[3 * param + 0] +
                stmParams.dfdpos[1] * stmParams.D[3 * param + 1] + 
                stmParams.dfdpos[2] * stmParams.D[3 * param + 2] +
                stmParams.dfdvel[0] * stmParams.Ddot[3 * param + 0] +
                stmParams.dfdvel[1] * stmParams.Ddot[3 * param + 1] +
                stmParams.dfdvel[2] * stmParams.Ddot[3 * param + 2] + 
                stmParams.dfdpar[3 * param + 0];
            D2dot[3 * param + 1] = stmParams.dfdpos[3] * stmParams.D[3 * param + 0] +
                stmParams.dfdpos[4] * stmParams.D[3 * param + 1] + 
                stmParams.dfdpos[5] * stmParams.D[3 * param + 2] +
                stmParams.dfdvel[3] * stmParams.Ddot[3 * param + 0] +
                stmParams.dfdvel[4] * stmParams.Ddot[3 * param + 1] +
                stmParams.dfdvel[5] * stmParams.Ddot[3 * param + 2] + 
                stmParams.dfdpar[3 * param + 1];
            D2dot[3 * param + 2] = stmParams.dfdpos[6] * stmParams.D[3 * param + 0] +
                stmParams.dfdpos[7] * stmParams.D[3 * param + 1] + 
                stmParams.dfdpos[8] * stmParams.D[3 * param + 2] +
                stmParams.dfdvel[6] * stmParams.Ddot[3 * param + 0] +
                stmParams.dfdvel[7] * stmParams.Ddot[3 * param + 1] +
                stmParams.dfdvel[8] * stmParams.Ddot[3 * param + 2] + 
                stmParams.dfdpar[3 * param + 2];
            accInteg[stmStarti + 18 + 3 * param + 0] = D2dot[3 * param + 0];
            accInteg[stmStarti + 18 + 3 * param + 1] = D2dot[3 * param + 1];
            accInteg[stmStarti + 18 + 3 * param + 2] = D2dot[3 * param + 2];
        }
    }
}

/**
 * @param[inout] stmParams Structure containing the STM submatrices and their derivatives.
 * @param[in] gm Gravitational parameter.
 * @param[in] dx Position difference in the x direction [AU].
 * @param[in] dy Position difference in the y direction [AU].
 * @param[in] dz Position difference in the z direction [AU].
 */
void stm_newton(STMParameters &stmParams, const real &gm, const real &dx,
                const real &dy, const real &dz) {
    const real r = sqrt(dx * dx + dy * dy + dz * dz);
    const real r3 = r * r * r;
    const real r5 = r3 * r * r;
    stmParams.dfdpos[0] += gm * (3 * dx * dx / r5 - 1 / r3);
    stmParams.dfdpos[1] += gm * (3 * dx * dy / r5);
    stmParams.dfdpos[2] += gm * (3 * dx * dz / r5);
    stmParams.dfdpos[3] += gm * (3 * dx * dy / r5);
    stmParams.dfdpos[4] += gm * (3 * dy * dy / r5 - 1 / r3);
    stmParams.dfdpos[5] += gm * (3 * dy * dz / r5);
    stmParams.dfdpos[6] += gm * (3 * dx * dz / r5);
    stmParams.dfdpos[7] += gm * (3 * dy * dz / r5);
    stmParams.dfdpos[8] += gm * (3 * dz * dz / r5 - 1 / r3);
}

/**
 * @param[inout] stmParams Structure containing the STM submatrices and their derivatives.
 * @param[in] gm Gravitational parameter.
 * @param[in] c Speed of light [AU/day].
 * @param[in] beta PPN parameter.
 * @param[in] gamma PPN parameter.
 * @param[in] dx Position difference in the x direction [AU].
 * @param[in] dy Position difference in the y direction [AU].
 * @param[in] dz Position difference in the z direction [AU].
 * @param[in] dvx Velocity difference in the x direction [AU/day].
 * @param[in] dvy Velocity difference in the y direction [AU/day].
 * @param[in] dvz Velocity difference in the z direction [AU/day].
 */
void stm_ppn_simple(STMParameters &stmParams, const real &gm, const real &c,
                    const real &beta, const real &gamma, const real &dx,
                    const real &dy, const real &dz, const real &dvx,
                    const real &dvy, const real &dvz) {
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

    stmParams.dfdpos[0] += dfac1dr*dx/r * (fac2*dx + fac3*dvx) + fac1*(fac2 + dfac2dr*dx*dx/r + 2*(1+gamma)*dvx*dvx);
    stmParams.dfdpos[1] += dfac1dr*dy/r * (fac2*dx + fac3*dvx) + fac1*(dfac2dr*dx*dy/r + 2*(1+gamma)*dvx*dvy);
    stmParams.dfdpos[2] += dfac1dr*dz/r * (fac2*dx + fac3*dvx) + fac1*(dfac2dr*dx*dz/r + 2*(1+gamma)*dvx*dvz);
    stmParams.dfdpos[3] += dfac1dr*dx/r * (fac2*dy + fac3*dvy) + fac1*(dfac2dr*dy*dx/r + 2*(1+gamma)*dvy*dvx);
    stmParams.dfdpos[4] += dfac1dr*dy/r * (fac2*dy + fac3*dvy) + fac1*(fac2 + dfac2dr*dy*dy/r + 2*(1+gamma)*dvy*dvy);
    stmParams.dfdpos[5] += dfac1dr*dz/r * (fac2*dy + fac3*dvy) + fac1*(dfac2dr*dy*dz/r + 2*(1+gamma)*dvy*dvz);
    stmParams.dfdpos[6] += dfac1dr*dx/r * (fac2*dz + fac3*dvz) + fac1*(dfac2dr*dz*dx/r + 2*(1+gamma)*dvz*dvx);
    stmParams.dfdpos[7] += dfac1dr*dy/r * (fac2*dz + fac3*dvz) + fac1*(dfac2dr*dz*dy/r + 2*(1+gamma)*dvz*dvy);
    stmParams.dfdpos[8] += dfac1dr*dz/r * (fac2*dz + fac3*dvz) + fac1*(fac2 + dfac2dr*dz*dz/r + 2*(1+gamma)*dvz*dvz);

    stmParams.dfdvel[0] += fac1*(-2*gamma*dvx*dx + fac3 + 2*(1+gamma)*dx*dvx);
    stmParams.dfdvel[1] += fac1*(-2*gamma*dvy*dx + 2*(1+gamma)*dy*dvx);
    stmParams.dfdvel[2] += fac1*(-2*gamma*dvz*dx + 2*(1+gamma)*dz*dvx);
    stmParams.dfdvel[3] += fac1*(-2*gamma*dvx*dy + 2*(1+gamma)*dx*dvy);
    stmParams.dfdvel[4] += fac1*(-2*gamma*dvy*dy + fac3 + 2*(1+gamma)*dy*dvy);
    stmParams.dfdvel[5] += fac1*(-2*gamma*dvz*dy + 2*(1+gamma)*dz*dvy);
    stmParams.dfdvel[6] += fac1*(-2*gamma*dvx*dz + 2*(1+gamma)*dx*dvz);
    stmParams.dfdvel[7] += fac1*(-2*gamma*dvy*dz + 2*(1+gamma)*dy*dvz);
    stmParams.dfdvel[8] += fac1*(-2*gamma*dvz*dz + fac3 + 2*(1+gamma)*dz*dvz);
}

/**
 * @param[inout] stmParams Structure containing the STM submatrices and their derivatives.
 * @param[in] gm Gravitational parameter.
 * @param[in] J2 Oblateness coefficient.
 * @param[in] dxBody Body fixed position difference in the x direction [AU].
 * @param[in] dyBody Body fixed position difference in the y direction [AU].
 * @param[in] dzBody Body fixed position difference in the z direction [AU].
 * @param[in] radius Radius of the body [AU].
 * @param[in] sinRA Sine of the right ascension of the body pole.
 * @param[in] cosRA Cosine of the right ascension of the body pole.
 * @param[in] sinDec Sine of the declination of the body pole.
 * @param[in] cosDec Cosine of the declination of the body pole.
 * @param[in] smoothing_threshold Threshold for the J2 smoothing function inside the body.
 */
void stm_J2(STMParameters &stmParams, const real &gm, const real &J2,
            const real &dxBody, const real &dyBody, const real &dzBody,
            const real &radius, const real &sinRA, const real &cosRA,
            const real &sinDec, const real &cosDec,
            const real &smoothing_threshold) {
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
    real *dfdposJ2 = new real[9];
    mat3_mat3_mul(dfBodydposBody, dposBodydpos, dfBodydpos);
    mat3_mat3_mul(dfdfBody, dfBodydpos, dfdposJ2);

    stmParams.dfdpos[0] += dfdposJ2[0];
    stmParams.dfdpos[1] += dfdposJ2[1];
    stmParams.dfdpos[2] += dfdposJ2[2];
    stmParams.dfdpos[3] += dfdposJ2[3];
    stmParams.dfdpos[4] += dfdposJ2[4];
    stmParams.dfdpos[5] += dfdposJ2[5];
    stmParams.dfdpos[6] += dfdposJ2[6];
    stmParams.dfdpos[7] += dfdposJ2[7];
    stmParams.dfdpos[8] += dfdposJ2[8];
}

/**
 * @param[inout] stmParams Structure containing the STM submatrices and their derivatives.
 * @param[in] g Non-gravitational acceleration scaling function.
 * @param[in] ngParams Structure containing the non-gravitational acceleration parameters.
 * @param[in] dx Heliocentric position difference in the x direction [AU].
 * @param[in] dy Heliocentric position difference in the y direction [AU].
 * @param[in] dz Heliocentric position difference in the z direction [AU].
 * @param[in] dvx Heliocentric velocity difference in the x direction [AU/day].
 * @param[in] dvy Heliocentric velocity difference in the y direction [AU/day].
 * @param[in] dvz Heliocentric velocity difference in the z direction [AU/day].
 * @param[in] rVec Relative position vector.
 * @param[in] nVec Relative angular momentum vector.
 */
void stm_nongrav(STMParameters &stmParams, const real &g,
                 const NongravParameters &ngParams, const real &dx, const real &dy,
                 const real &dz, const real &dvx, const real &dvy, const real &dvz,
                 real *rVec, real *nVec) {
    const real a1 = ngParams.a1;
    const real a2 = ngParams.a2;
    const real a3 = ngParams.a3;
    const bool a1Est = ngParams.a1Est;
    const bool a2Est = ngParams.a2Est;
    const bool a3Est = ngParams.a3Est;
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

    stmParams.dfdpos[0] += a1*(dgdx*rHat[0] - g*rxr3*dx + g/rNorm)
                + a2*(dgdx*tHat[0] - g*txt3*(2*dx*vDotT - tVec[0]*rDotV) + g*(dx*dvx - rDotV)/tNorm)
                + a3*(dgdx*nHat[0] - g*nxn3*(dx*vDotV - dvx*rDotV));
    stmParams.dfdpos[1] += a1*(dgdy*rHat[0] - g*rxr3*dy)
                + a2*(dgdy*tHat[0] - g*txt3*(2*dy*vDotT - tVec[1]*rDotV) + g*(2*dy*dvx - dx*dvy)/tNorm)
                + a3*(dgdy*nHat[0] - g*nxn3*(dy*vDotV - dvy*rDotV) + g*dvz/nNorm);
    stmParams.dfdpos[2] += a1*(dgdz*rHat[0] - g*rxr3*dz)
                + a2*(dgdz*tHat[0] - g*txt3*(2*dz*vDotT - tVec[2]*rDotV) + g*(2*dz*dvx - dx*dvz)/tNorm)
                + a3*(dgdz*nHat[0] - g*nxn3*(dz*vDotV - dvz*rDotV) - g*dvy/nNorm);
    stmParams.dfdpos[3] += a1*(dgdx*rHat[1] - g*ryr3*dx)
                + a2*(dgdx*tHat[1] - g*tyt3*(2*dx*vDotT - tVec[0]*rDotV) + g*(2*dx*dvy - dy*dvx)/tNorm)
                + a3*(dgdx*nHat[1] - g*nyn3*(dx*vDotV - dvx*rDotV) - g*dvz/nNorm);
    stmParams.dfdpos[4] += a1*(dgdy*rHat[1] - g*ryr3*dy + g/rNorm)
                + a2*(dgdy*tHat[1] - g*tyt3*(2*dy*vDotT - tVec[1]*rDotV) + g*(dy*dvy - rDotV)/tNorm)
                + a3*(dgdy*nHat[1] - g*nyn3*(dy*vDotV - dvy*rDotV));
    stmParams.dfdpos[5] += a1*(dgdz*rHat[1] - g*ryr3*dz)
                + a2*(dgdz*tHat[1] - g*tyt3*(2*dz*vDotT - tVec[2]*rDotV) + g*(2*dz*dvy - dy*dvz)/tNorm)
                + a3*(dgdz*nHat[1] - g*nyn3*(dz*vDotV - dvz*rDotV) + g*dvx/nNorm);
    stmParams.dfdpos[6] += a1*(dgdx*rHat[2] - g*rzr3*dx)
                + a2*(dgdx*tHat[2] - g*tzt3*(2*dx*vDotT - tVec[0]*rDotV) + g*(2*dx*dvz - dz*dvx)/tNorm)
                + a3*(dgdx*nHat[2] - g*nzn3*(dx*vDotV - dvx*rDotV) + g*dvy/nNorm);
    stmParams.dfdpos[7] += a1*(dgdy*rHat[2] - g*rzr3*dy)
                + a2*(dgdy*tHat[2] - g*tzt3*(2*dy*vDotT - tVec[1]*rDotV) + g*(2*dy*dvz - dz*dvy)/tNorm)
                + a3*(dgdy*nHat[2] - g*nzn3*(dy*vDotV - dvy*rDotV) - g*dvx/nNorm);
    stmParams.dfdpos[8] += a1*(dgdz*rHat[2] - g*rzr3*dz + g/rNorm)
                + a2*(dgdz*tHat[2] - g*tzt3*(2*dz*vDotT - tVec[2]*rDotV) + g*(dz*dvz - rDotV)/tNorm)
                + a3*(dgdz*nHat[2] - g*nzn3*(dz*vDotV - dvz*rDotV));
    
    const real r2 = rNorm * rNorm;
    stmParams.dfdvel[0] += g*(a2*((dy*dy + dz*dz)/tNorm - txt3*r2*tVec[0]) - a3*(nxn3*(r2*dvx - dx*rDotV)));
    stmParams.dfdvel[1] += g*(a2*(-dx*dy/tNorm - tyt3*r2*tVec[0]) - a3*(nxn3*(r2*dvy - dy*rDotV) + dz/nNorm));
    stmParams.dfdvel[2] += g*(a2*(-dx*dz/tNorm - tzt3*r2*tVec[0]) - a3*(nxn3*(r2*dvz - dz*rDotV) - dy/nNorm));
    stmParams.dfdvel[3] += g*(a2*(-dy*dx/tNorm - txt3*r2*tVec[1]) - a3*(nyn3*(r2*dvx - dx*rDotV) - dz/nNorm));
    stmParams.dfdvel[4] += g*(a2*((dx*dx + dz*dz)/tNorm - tyt3*r2*tVec[1]) - a3*(nyn3*(r2*dvy - dy*rDotV)));
    stmParams.dfdvel[5] += g*(a2*(-dy*dz/tNorm - tzt3*r2*tVec[1]) - a3*(nyn3*(r2*dvz - dz*rDotV) + dx/nNorm));
    stmParams.dfdvel[6] += g*(a2*(-dz*dx/tNorm - txt3*r2*tVec[2]) - a3*(nzn3*(r2*dvx - dx*rDotV) + dy/nNorm));
    stmParams.dfdvel[7] += g*(a2*(-dz*dy/tNorm - tyt3*r2*tVec[2]) - a3*(nzn3*(r2*dvy - dy*rDotV) - dx/nNorm));
    stmParams.dfdvel[8] += g*(a2*((dx*dx + dy*dy)/tNorm - tzt3*r2*tVec[2]) - a3*(nzn3*(r2*dvz - dz*rDotV)));

    // Case 1, only one of a1, a2, or a3 is estimated
    if (a1Est && !a2Est && !a3Est) {
        stmParams.dfdpar[0] += g*rHat[0];
        stmParams.dfdpar[1] += g*rHat[1];
        stmParams.dfdpar[2] += g*rHat[2];
    } else if (!a1Est && a2Est && !a3Est) {
        stmParams.dfdpar[0] += g*tHat[0];
        stmParams.dfdpar[1] += g*tHat[1];
        stmParams.dfdpar[2] += g*tHat[2];
    } else if (!a1Est && !a2Est && a3Est) {
        stmParams.dfdpar[0] += g*nHat[0];
        stmParams.dfdpar[1] += g*nHat[1];
        stmParams.dfdpar[2] += g*nHat[2];
    }
    // Case 2, two of a1, a2, or a3 are estimated
    if (a1Est && a2Est && !a3Est) {
        stmParams.dfdpar[0] += g*rHat[0];
        stmParams.dfdpar[1] += g*rHat[1];
        stmParams.dfdpar[2] += g*rHat[2];
        stmParams.dfdpar[3] += g*tHat[0];
        stmParams.dfdpar[4] += g*tHat[1];
        stmParams.dfdpar[5] += g*tHat[2];
    } else if (a1Est && !a2Est && a3Est) {
        stmParams.dfdpar[0] += g*rHat[0];
        stmParams.dfdpar[1] += g*rHat[1];
        stmParams.dfdpar[2] += g*rHat[2];
        stmParams.dfdpar[3] += g*nHat[0];
        stmParams.dfdpar[4] += g*nHat[1];
        stmParams.dfdpar[5] += g*nHat[2];
    } else if (!a1Est && a2Est && a3Est) {
        stmParams.dfdpar[0] += g*tHat[0];
        stmParams.dfdpar[1] += g*tHat[1];
        stmParams.dfdpar[2] += g*tHat[2];
        stmParams.dfdpar[3] += g*nHat[0];
        stmParams.dfdpar[4] += g*nHat[1];
        stmParams.dfdpar[5] += g*nHat[2];
    }
    // Case 3, all of a1, a2, and a3 are estimated
    if (a1Est && a2Est && a3Est) {
        stmParams.dfdpar[0] += g*rHat[0];
        stmParams.dfdpar[1] += g*rHat[1];
        stmParams.dfdpar[2] += g*rHat[2];
        stmParams.dfdpar[3] += g*tHat[0];
        stmParams.dfdpar[4] += g*tHat[1];
        stmParams.dfdpar[5] += g*tHat[2];
        stmParams.dfdpar[6] += g*nHat[0];
        stmParams.dfdpar[7] += g*nHat[1];
        stmParams.dfdpar[8] += g*nHat[2];
    }
}

/**
 * @param[inout] stmParams Structure containing the STM submatrices and their derivatives.
 * @param[in] propSim PropSimulation object.
 * @param[in] eventIdx Index of the continuous event.
 * @param[in] tPastEvent Time since the continuous event started [days].
 * @param[in] postFac Exponential decay factor for the continuous event.
 */
void stm_continuous_event(STMParameters &stmParams,
                          const PropSimulation *propSim, const size_t &eventIdx,
                          const real &tPastEvent, const real &postFac){
    const Event *thisEvent = &propSim->eventMngr.continuousEvents[eventIdx];
    const IntegBody *body = &propSim->integBodies[thisEvent->bodyIndex];
    const size_t numNongravs = body->ngParams.a1Est + body->ngParams.a2Est + body->ngParams.a3Est;
    size_t dfdparIdx = 3*numNongravs;
    for (size_t i = 0; i < propSim->eventMngr.impulsiveEvents.size(); i++) {
        const bool sameBody = propSim->eventMngr.impulsiveEvents[i].bodyIndex == thisEvent->bodyIndex;
        const bool hasStarted = propSim->eventMngr.impulsiveEvents[i].hasStarted;
        const bool deltaVEst = propSim->eventMngr.impulsiveEvents[i].deltaVEst;
        const bool multiplierEst = propSim->eventMngr.impulsiveEvents[i].multiplierEst;
        if (sameBody && hasStarted && multiplierEst) {
            dfdparIdx += 3*1; // only multiplier
        } else if (sameBody && hasStarted && deltaVEst) {
            dfdparIdx += 3*3; // full deltaV vector
        }
    }
    for (size_t i = 0; i < eventIdx; i++) {
        const bool sameBody = propSim->eventMngr.continuousEvents[i].bodyIndex == thisEvent->bodyIndex;
        const bool hasStarted = propSim->eventMngr.continuousEvents[i].hasStarted;
        const bool expAccel0Est = propSim->eventMngr.continuousEvents[i].expAccel0Est;
        const bool tauEst = propSim->eventMngr.continuousEvents[i].tauEst;
        if (sameBody && hasStarted && expAccel0Est) {
            dfdparIdx += 3*3; // expAccel0 vector
        }
        if (sameBody && hasStarted && tauEst) {
            dfdparIdx += 3*1; // time constant
        }
    }
    if (thisEvent->isContinuous && thisEvent->expAccel0Est) {
        stmParams.dfdpar[dfdparIdx + 0] = postFac;
        stmParams.dfdpar[dfdparIdx + 4] = postFac;
        stmParams.dfdpar[dfdparIdx + 8] = postFac;
        dfdparIdx += 9;
    }
    if (thisEvent->isContinuous && thisEvent->tauEst) {
        const real postFac2 = tPastEvent/thisEvent->tau/thisEvent->tau;
        for (size_t i = 0; i < 3; i++) {
            stmParams.dfdpar[dfdparIdx + i] = -thisEvent->expAccel0[i]*postFac*postFac2;
        }
    }
}
