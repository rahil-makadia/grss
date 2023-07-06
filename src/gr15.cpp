#include "gr15.h"

real get_initial_timestep(const real &t, const std::vector<real> &xInteg0,
                          const ForceParameters &forceParams,
                          IntegrationParameters &integParams,
                          const Constants &consts) {
    real dt;
    if (integParams.dt0 != 0.0) {
        dt = fabs(integParams.dt0);
        if (integParams.tf < integParams.t0) {
            dt *= -1.0;
            integParams.dtMax = -fabs(integParams.dtMax);
            integParams.dtMin = -fabs(integParams.dtMin);
        }
        return dt;
    }
    int order = 15;
    real absMaxPos0, absMaxAcc0, absMaxAcc1Minus0;
    real dtTemp0, dtTemp1;
    std::vector<real> posInteg0(3 * integParams.nInteg, 0.0);
    std::vector<real> accInteg1Minus0(3 * integParams.nInteg, 0.0);
    std::vector<real> xIntegNext(6 * integParams.nInteg, 0.0);

    std::vector<real> accInteg0 =
        get_state_der(t, xInteg0, forceParams, integParams, consts);
    for (size_t i = 0; i < integParams.nInteg; i++) {
        for (size_t j = 0; j < 3; j++) {
            posInteg0[3 * i + j] = xInteg0[6 * i + j];
        }
    }
    vabs_max(posInteg0, absMaxPos0);
    vabs_max(accInteg0, absMaxAcc0);
    if (absMaxPos0 < 1.0e-5 || absMaxAcc0 < 1e-5) {
        dtTemp0 = 1.0e-6;
    } else {
        dtTemp0 = 0.01 * (absMaxPos0 / absMaxAcc0);
    }
    // propagate xInteg0 to xIntegNext using an Euler step and a timestep of
    // dtTemp0
    for (size_t i = 0; i < integParams.nInteg; i++) {
        for (size_t j = 0; j < 3; j++) {
            xIntegNext[6 * i + j] =
                xInteg0[6 * i + j] + dtTemp0 * xInteg0[6 * i + j + 3];
            xIntegNext[6 * i + j + 3] =
                xInteg0[6 * i + j + 3] + dtTemp0 * accInteg0[3 * i + j];
        }
    }
    std::vector<real> accInteg1 = get_state_der(
        t + dtTemp0, xIntegNext, forceParams, integParams, consts);
    vsub(accInteg1, accInteg0, accInteg1Minus0);
    vabs_max(accInteg1Minus0, absMaxAcc1Minus0);
    if (fmax(absMaxAcc0, absMaxAcc1Minus0) <= 1e-15) {
        dtTemp1 = fmax(1.0e-6, dtTemp0 * 1e-3);
    } else {
        dtTemp1 =
            pow(0.01 / fmax(absMaxAcc0, absMaxAcc1Minus0), 1.0 / (order + 1));
    }
    dt = fmin(100 * dtTemp0, dtTemp1);
    if (fabs(integParams.tf - integParams.t0) < dt) {
        dt = fabs(integParams.tf - integParams.t0);
    }
    if (integParams.tf < integParams.t0) {
        dt *= -1.0;
        integParams.dtMax = -fabs(integParams.dtMax);
        integParams.dtMin = -fabs(integParams.dtMin);
    }
    return dt;
}

void approx_xInteg(const std::vector<real> &xInteg0,
                   const std::vector<real> &accInteg0,
                   std::vector<real> &xIntegNext, const real &dt, const real &h,
                   const std::vector<std::vector<real>> &b,
                   const size_t &nInteg) {
    for (size_t i = 0; i < nInteg; i++) {
        xIntegNext[6*i]   = xInteg0[6*i]   + dt * h * (xInteg0[6*i+3] + dt * h * (accInteg0[3*i]   + h * (b[0][3*i]   / 0.3e1 + h * (b[1][3*i]   / 0.6e1 + h * (b[2][3*i]   / 0.10e2 + h * (b[3][3*i]   / 0.15e2 + h * (b[4][3*i]   / 0.21e2 + h * (b[5][3*i]   / 0.28e2 + h * b[6][3*i]   / 0.36e2))))))) / 0.2e1);
        xIntegNext[6*i+1] = xInteg0[6*i+1] + dt * h * (xInteg0[6*i+4] + dt * h * (accInteg0[3*i+1] + h * (b[0][3*i+1] / 0.3e1 + h * (b[1][3*i+1] / 0.6e1 + h * (b[2][3*i+1] / 0.10e2 + h * (b[3][3*i+1] / 0.15e2 + h * (b[4][3*i+1] / 0.21e2 + h * (b[5][3*i+1] / 0.28e2 + h * b[6][3*i+1] / 0.36e2))))))) / 0.2e1);
        xIntegNext[6*i+2] = xInteg0[6*i+2] + dt * h * (xInteg0[6*i+5] + dt * h * (accInteg0[3*i+2] + h * (b[0][3*i+2] / 0.3e1 + h * (b[1][3*i+2] / 0.6e1 + h * (b[2][3*i+2] / 0.10e2 + h * (b[3][3*i+2] / 0.15e2 + h * (b[4][3*i+2] / 0.21e2 + h * (b[5][3*i+2] / 0.28e2 + h * b[6][3*i+2] / 0.36e2))))))) / 0.2e1);

        xIntegNext[6*i+3] = xInteg0[6*i+3] + dt * h * (accInteg0[3*i]   + h * (b[0][3*i]   / 0.2e1 + h * (b[1][3*i]   / 0.3e1 + h * (b[2][3*i]   / 0.4e1 + h * (b[3][3*i]   / 0.5e1 + h * (b[4][3*i]   / 0.6e1 + h * (b[5][3*i]   / 0.7e1 + h * b[6][3*i]   / 0.8e1)))))));
        xIntegNext[6*i+4] = xInteg0[6*i+4] + dt * h * (accInteg0[3*i+1] + h * (b[0][3*i+1] / 0.2e1 + h * (b[1][3*i+1] / 0.3e1 + h * (b[2][3*i+1] / 0.4e1 + h * (b[3][3*i+1] / 0.5e1 + h * (b[4][3*i+1] / 0.6e1 + h * (b[5][3*i+1] / 0.7e1 + h * b[6][3*i+1] / 0.8e1)))))));
        xIntegNext[6*i+5] = xInteg0[6*i+5] + dt * h * (accInteg0[3*i+2] + h * (b[0][3*i+2] / 0.2e1 + h * (b[1][3*i+2] / 0.3e1 + h * (b[2][3*i+2] / 0.4e1 + h * (b[3][3*i+2] / 0.5e1 + h * (b[4][3*i+2] / 0.6e1 + h * (b[5][3*i+2] / 0.7e1 + h * b[6][3*i+2] / 0.8e1)))))));
    }
}

void compute_g_and_b(const std::vector<std::vector<real>> &AccIntegArr,
                     const size_t &hIdx, std::vector<std::vector<real>> &g,
                     std::vector<std::vector<real>> &b, const size_t &dim) {
    const std::vector<real> Acc1 = AccIntegArr[0];
    const std::vector<real> Acc2 = AccIntegArr[1];
    const std::vector<real> Acc3 = AccIntegArr[2];
    const std::vector<real> Acc4 = AccIntegArr[3];
    const std::vector<real> Acc5 = AccIntegArr[4];
    const std::vector<real> Acc6 = AccIntegArr[5];
    const std::vector<real> Acc7 = AccIntegArr[6];
    const std::vector<real> Acc8 = AccIntegArr[7];

    for (size_t i=0; i<dim; i++){
        if (hIdx == 1) {
            g[0][i] = (Acc2[i] - Acc1[i]) * rMat[1][0];
        } else if (hIdx == 2) {
            g[0][i] = (Acc2[i] - Acc1[i]) * rMat[1][0];
            g[1][i] = ((Acc3[i] - Acc1[i]) * rMat[2][0] - g[0][i]) * rMat[2][1];
        } else if (hIdx == 3) {
            g[0][i] = (Acc2[i] - Acc1[i]) * rMat[1][0];
            g[1][i] = ((Acc3[i] - Acc1[i]) * rMat[2][0] - g[0][i]) * rMat[2][1];
            g[2][i] = (((Acc4[i] - Acc1[i]) * rMat[3][0] - g[0][i]) * rMat[3][1] - g[1][i]) * rMat[3][2];
        } else if (hIdx == 4) {
            g[0][i] = (Acc2[i] - Acc1[i]) * rMat[1][0];
            g[1][i] = ((Acc3[i] - Acc1[i]) * rMat[2][0] - g[0][i]) * rMat[2][1];
            g[2][i] = (((Acc4[i] - Acc1[i]) * rMat[3][0] - g[0][i]) * rMat[3][1] - g[1][i]) * rMat[3][2];
            g[3][i] = ((((Acc5[i] - Acc1[i]) * rMat[4][0] - g[0][i]) * rMat[4][1] - g[1][i]) * rMat[4][2] - g[2][i]) * rMat[4][3];
        } else if (hIdx == 5) {
            g[0][i] = (Acc2[i] - Acc1[i]) * rMat[1][0];
            g[1][i] = ((Acc3[i] - Acc1[i]) * rMat[2][0] - g[0][i]) * rMat[2][1];
            g[2][i] = (((Acc4[i] - Acc1[i]) * rMat[3][0] - g[0][i]) * rMat[3][1] - g[1][i]) * rMat[3][2];
            g[3][i] = ((((Acc5[i] - Acc1[i]) * rMat[4][0] - g[0][i]) * rMat[4][1] - g[1][i]) * rMat[4][2] - g[2][i]) * rMat[4][3];
            g[4][i] = (((((Acc6[i] - Acc1[i]) * rMat[5][0] - g[0][i]) * rMat[5][1] - g[1][i]) * rMat[5][2] - g[2][i]) * rMat[5][3] - g[3][i]) * rMat[5][4];
        } else if (hIdx == 6) {
            g[0][i] = (Acc2[i] - Acc1[i]) * rMat[1][0];
            g[1][i] = ((Acc3[i] - Acc1[i]) * rMat[2][0] - g[0][i]) * rMat[2][1];
            g[2][i] = (((Acc4[i] - Acc1[i]) * rMat[3][0] - g[0][i]) * rMat[3][1] - g[1][i]) * rMat[3][2];
            g[3][i] = ((((Acc5[i] - Acc1[i]) * rMat[4][0] - g[0][i]) * rMat[4][1] - g[1][i]) * rMat[4][2] - g[2][i]) * rMat[4][3];
            g[4][i] = (((((Acc6[i] - Acc1[i]) * rMat[5][0] - g[0][i]) * rMat[5][1] - g[1][i]) * rMat[5][2] - g[2][i]) * rMat[5][3] - g[3][i]) * rMat[5][4];
            g[5][i] = ((((((Acc7[i] - Acc1[i]) * rMat[6][0] - g[0][i]) * rMat[6][1] - g[1][i]) * rMat[6][2] - g[2][i]) * rMat[6][3] - g[3][i]) * rMat[6][4] - g[4][i]) * rMat[6][5];
        } else if (hIdx == 7) {
            g[0][i] = (Acc2[i] - Acc1[i]) * rMat[1][0];
            g[1][i] = ((Acc3[i] - Acc1[i]) * rMat[2][0] - g[0][i]) * rMat[2][1];
            g[2][i] = (((Acc4[i] - Acc1[i]) * rMat[3][0] - g[0][i]) * rMat[3][1] - g[1][i]) * rMat[3][2];
            g[3][i] = ((((Acc5[i] - Acc1[i]) * rMat[4][0] - g[0][i]) * rMat[4][1] - g[1][i]) * rMat[4][2] - g[2][i]) * rMat[4][3];
            g[4][i] = (((((Acc6[i] - Acc1[i]) * rMat[5][0] - g[0][i]) * rMat[5][1] - g[1][i]) * rMat[5][2] - g[2][i]) * rMat[5][3] - g[3][i]) * rMat[5][4];
            g[5][i] = ((((((Acc7[i] - Acc1[i]) * rMat[6][0] - g[0][i]) * rMat[6][1] - g[1][i]) * rMat[6][2] - g[2][i]) * rMat[6][3] - g[3][i]) * rMat[6][4] - g[4][i]) * rMat[6][5];
            g[6][i] = (((((((Acc8[i] - Acc1[i]) * rMat[7][0] - g[0][i]) * rMat[7][1] - g[1][i]) * rMat[7][2] - g[2][i]) * rMat[7][3] - g[3][i]) * rMat[7][4] - g[4][i]) * rMat[7][5] - g[5][i]) * rMat[7][6];
        }
    }
    for (size_t i=0; i<dim; i++){
        if (hIdx == 1) {
            b[0][i] = cMat[0][0]*g[0][i] + cMat[1][0]*g[1][i] + cMat[2][0]*g[2][i] + cMat[3][0]*g[3][i] + cMat[4][0]*g[4][i] + cMat[5][0]*g[5][i] + cMat[6][0]*g[6][i];
        } else if (hIdx == 2) {
            b[0][i] = cMat[0][0]*g[0][i] + cMat[1][0]*g[1][i] + cMat[2][0]*g[2][i] + cMat[3][0]*g[3][i] + cMat[4][0]*g[4][i] + cMat[5][0]*g[5][i] + cMat[6][0]*g[6][i];
            b[1][i] =                    + cMat[1][1]*g[1][i] + cMat[2][1]*g[2][i] + cMat[3][1]*g[3][i] + cMat[4][1]*g[4][i] + cMat[5][1]*g[5][i] + cMat[6][1]*g[6][i];
        } else if (hIdx == 3) {
            b[0][i] = cMat[0][0]*g[0][i] + cMat[1][0]*g[1][i] + cMat[2][0]*g[2][i] + cMat[3][0]*g[3][i] + cMat[4][0]*g[4][i] + cMat[5][0]*g[5][i] + cMat[6][0]*g[6][i];
            b[1][i] =                    + cMat[1][1]*g[1][i] + cMat[2][1]*g[2][i] + cMat[3][1]*g[3][i] + cMat[4][1]*g[4][i] + cMat[5][1]*g[5][i] + cMat[6][1]*g[6][i];
            b[2][i] =                                         + cMat[2][2]*g[2][i] + cMat[3][2]*g[3][i] + cMat[4][2]*g[4][i] + cMat[5][2]*g[5][i] + cMat[6][2]*g[6][i];
        } else if (hIdx == 4) {
            b[0][i] = cMat[0][0]*g[0][i] + cMat[1][0]*g[1][i] + cMat[2][0]*g[2][i] + cMat[3][0]*g[3][i] + cMat[4][0]*g[4][i] + cMat[5][0]*g[5][i] + cMat[6][0]*g[6][i];
            b[1][i] =                    + cMat[1][1]*g[1][i] + cMat[2][1]*g[2][i] + cMat[3][1]*g[3][i] + cMat[4][1]*g[4][i] + cMat[5][1]*g[5][i] + cMat[6][1]*g[6][i];
            b[2][i] =                                         + cMat[2][2]*g[2][i] + cMat[3][2]*g[3][i] + cMat[4][2]*g[4][i] + cMat[5][2]*g[5][i] + cMat[6][2]*g[6][i];
            b[3][i] =                                                              + cMat[3][3]*g[3][i] + cMat[4][3]*g[4][i] + cMat[5][3]*g[5][i] + cMat[6][3]*g[6][i];
        } else if (hIdx == 5) {
            b[0][i] = cMat[0][0]*g[0][i] + cMat[1][0]*g[1][i] + cMat[2][0]*g[2][i] + cMat[3][0]*g[3][i] + cMat[4][0]*g[4][i] + cMat[5][0]*g[5][i] + cMat[6][0]*g[6][i];
            b[1][i] =                    + cMat[1][1]*g[1][i] + cMat[2][1]*g[2][i] + cMat[3][1]*g[3][i] + cMat[4][1]*g[4][i] + cMat[5][1]*g[5][i] + cMat[6][1]*g[6][i];
            b[2][i] =                                         + cMat[2][2]*g[2][i] + cMat[3][2]*g[3][i] + cMat[4][2]*g[4][i] + cMat[5][2]*g[5][i] + cMat[6][2]*g[6][i];
            b[3][i] =                                                              + cMat[3][3]*g[3][i] + cMat[4][3]*g[4][i] + cMat[5][3]*g[5][i] + cMat[6][3]*g[6][i];
            b[4][i] =                                                                                   + cMat[4][4]*g[4][i] + cMat[5][4]*g[5][i] + cMat[6][4]*g[6][i];
        } else if (hIdx == 6) {
            b[0][i] = cMat[0][0]*g[0][i] + cMat[1][0]*g[1][i] + cMat[2][0]*g[2][i] + cMat[3][0]*g[3][i] + cMat[4][0]*g[4][i] + cMat[5][0]*g[5][i] + cMat[6][0]*g[6][i];
            b[1][i] =                    + cMat[1][1]*g[1][i] + cMat[2][1]*g[2][i] + cMat[3][1]*g[3][i] + cMat[4][1]*g[4][i] + cMat[5][1]*g[5][i] + cMat[6][1]*g[6][i];
            b[2][i] =                                         + cMat[2][2]*g[2][i] + cMat[3][2]*g[3][i] + cMat[4][2]*g[4][i] + cMat[5][2]*g[5][i] + cMat[6][2]*g[6][i];
            b[3][i] =                                                              + cMat[3][3]*g[3][i] + cMat[4][3]*g[4][i] + cMat[5][3]*g[5][i] + cMat[6][3]*g[6][i];
            b[4][i] =                                                                                   + cMat[4][4]*g[4][i] + cMat[5][4]*g[5][i] + cMat[6][4]*g[6][i];
            b[5][i] =                                                                                                        + cMat[5][5]*g[5][i] + cMat[6][5]*g[6][i];
        } else if (hIdx == 7) {
            b[0][i] = cMat[0][0]*g[0][i] + cMat[1][0]*g[1][i] + cMat[2][0]*g[2][i] + cMat[3][0]*g[3][i] + cMat[4][0]*g[4][i] + cMat[5][0]*g[5][i] + cMat[6][0]*g[6][i];
            b[1][i] =                    + cMat[1][1]*g[1][i] + cMat[2][1]*g[2][i] + cMat[3][1]*g[3][i] + cMat[4][1]*g[4][i] + cMat[5][1]*g[5][i] + cMat[6][1]*g[6][i];
            b[2][i] =                                         + cMat[2][2]*g[2][i] + cMat[3][2]*g[3][i] + cMat[4][2]*g[4][i] + cMat[5][2]*g[5][i] + cMat[6][2]*g[6][i];
            b[3][i] =                                                              + cMat[3][3]*g[3][i] + cMat[4][3]*g[4][i] + cMat[5][3]*g[5][i] + cMat[6][3]*g[6][i];
            b[4][i] =                                                                                   + cMat[4][4]*g[4][i] + cMat[5][4]*g[5][i] + cMat[6][4]*g[6][i];
            b[5][i] =                                                                                                        + cMat[5][5]*g[5][i] + cMat[6][5]*g[6][i];
            b[6][i] =                                                                                                                             + cMat[6][6]*g[6][i];
        }
    }
}

void refine_b(std::vector<std::vector<real>> &b,
              std::vector<std::vector<real>> &e, const real &dtRatio,
              const size_t &dim, const size_t &timestepCounter) {
    std::vector<std::vector<real>> bDiff(7, std::vector<real>(dim, 0.0L));
    if (timestepCounter > 1){
        for (size_t i = 0; i < dim; i++){
            bDiff[0][i] = b[0][i] - e[0][i];
            bDiff[1][i] = b[1][i] - e[1][i];
            bDiff[2][i] = b[2][i] - e[2][i];
            bDiff[3][i] = b[3][i] - e[3][i];
            bDiff[4][i] = b[4][i] - e[4][i];
            bDiff[5][i] = b[5][i] - e[5][i];
            bDiff[6][i] = b[6][i] - e[6][i];
        }
    }

    real q = dtRatio;
    real q2 = q * q;
    real q3 = q2 * q;
    real q4 = q2 * q2;
    real q5 = q2 * q3;
    real q6 = q3 * q3;
    real q7 = q2 * q5;

    for (size_t i = 0; i < dim; i++) {
        e[0][i] = q  * (b[6][i] * 7.0  + b[5][i] * 6.0  + b[4][i] * 5.0  + b[3][i] * 4.0 + b[2][i] * 3.0 + b[1][i] * 2.0 + b[0][i]);
        e[1][i] = q2 * (b[6][i] * 21.0 + b[5][i] * 15.0 + b[4][i] * 10.0 + b[3][i] * 6.0 + b[2][i] * 3.0 + b[1][i]);
        e[2][i] = q3 * (b[6][i] * 35.0 + b[5][i] * 20.0 + b[4][i] * 10.0 + b[3][i] * 4.0 + b[2][i]);
        e[3][i] = q4 * (b[6][i] * 35.0 + b[5][i] * 15.0 + b[4][i] * 5.0  + b[3][i]);
        e[4][i] = q5 * (b[6][i] * 21.0 + b[5][i] * 6.0  + b[4][i]);
        e[5][i] = q6 * (b[6][i] * 7.0  + b[5][i]);
        e[6][i] = q7 * (b[6][i]);
    }

    for (size_t i = 0; i < dim; i++) {
        b[0][i] = e[0][i] + bDiff[0][i];
        b[1][i] = e[1][i] + bDiff[1][i];
        b[2][i] = e[2][i] + bDiff[2][i];
        b[3][i] = e[3][i] + bDiff[3][i];
        b[4][i] = e[4][i] + bDiff[4][i];
        b[5][i] = e[5][i] + bDiff[5][i];
        b[6][i] = e[6][i] + bDiff[6][i];
    }
}

void check_and_apply_events(propSimulation *propSim, const real &t,
                            real &tNextEvent, size_t &nextEventIdx,
                            std::vector<real> &xInteg) {
    while (nextEventIdx < propSim->events.size() && t == tNextEvent) {
        // apply events for the state just reached by the integrator
        real propDir;
        if (propSim->integParams.t0 < propSim->integParams.tf) {
            propDir = 1.0L;
        } else {
            propDir = -1.0L;
        }
        propSim->events[nextEventIdx].apply(t, xInteg, propDir);
        // update next event index and time
        nextEventIdx += 1;
        if (nextEventIdx < propSim->events.size()) {
            tNextEvent = propSim->events[nextEventIdx].t;
        } else {
            tNextEvent = propSim->integParams.tf;
        }
    }
}

void gr15(propSimulation *propSim) {
    real t = propSim->t;
    std::vector<real> xInteg0 = propSim->xInteg;
    // ForceParameters &forceParams = propSim->forceParams;
    // IntegrationParameters &integParams = propSim->integParams;
    size_t nh = 8;
    size_t dim = 3 * propSim->integParams.nInteg;
    // Constants &consts = propSim->consts;
    real dt = get_initial_timestep(t, xInteg0, propSim->forceParams,
                                   propSim->integParams, propSim->consts);
    propSim->integParams.timestepCounter = 0;
    std::vector<real> accInteg0;
    std::vector<real> xInteg(2 * dim, 0.0);
    std::vector<std::vector<real> > b_old(7, std::vector<real>(dim, 0.0));
    std::vector<std::vector<real> > b(7, std::vector<real>(dim, 0.0));
    std::vector<std::vector<real> > g(7, std::vector<real>(dim, 0.0));
    std::vector<std::vector<real> > e(7, std::vector<real>(dim, 0.0));
    std::vector<std::vector<real> > accIntegArr(nh,
                                                std::vector<real>(dim, 0.0));
    std::vector<real> b6Tilde(dim, 0.0);
    real b6TildeMax, accIntegArr7Max;
    real b6TildeEstim, b6Max, accIntegNextMax;
    real relError, dtReq;

    real tNextEvent = propSim->integParams.tf;
    static size_t nextEventIdx = 0;
    if (t == propSim->integParams.t0) {
        nextEventIdx = 0;
    }
    if (propSim->events.size() != 0) {
        tNextEvent = propSim->events[0].t;
    }
    check_and_apply_events(propSim, t, tNextEvent, nextEventIdx, xInteg0);
    if ((propSim->integParams.tf > propSim->integParams.t0 &&
         t + dt > tNextEvent) ||
        (propSim->integParams.tf < propSim->integParams.t0 &&
         t + dt < tNextEvent)) {
        dt = tNextEvent - t;
    }

    size_t PCmaxIter = 12;
    int maxLoops = 25;
    int loopCounter = 0;
    int keepStepping = 1;
    int oneStepDone = 0;
    if (propSim->integParams.t0 == propSim->integParams.tf) {
        keepStepping = 0;
    }
    while (keepStepping) {
        t = propSim->t;
        xInteg0 = propSim->xInteg;
        oneStepDone = 0;
        while (!oneStepDone) {
            accInteg0 = get_state_der(t, xInteg0, propSim->forceParams,
                                      propSim->integParams, propSim->consts);
            for (size_t PCidx = 1; PCidx < PCmaxIter; PCidx++) {
                for (size_t hIdx = 0; hIdx < nh; hIdx++) {
                    approx_xInteg(xInteg0, accInteg0, xInteg, dt, hVec[hIdx], b,
                                  propSim->integParams.nInteg);
                    accIntegArr[hIdx] = get_state_der(
                        t + hVec[hIdx] * dt, xInteg, propSim->forceParams,
                        propSim->integParams, propSim->consts);
                    compute_g_and_b(accIntegArr, hIdx, g, b, dim);
                }
                for (size_t i = 0; i < dim; i++) {
                    b6Tilde[i] = b[6][i] - b_old[6][i];
                }
                vabs_max(b6Tilde, b6TildeMax);
                vabs_max(accIntegArr[7], accIntegArr7Max);
                if (b6TildeMax / accIntegArr7Max < propSim->integParams.tolPC) {
                    break;
                }
                b_old = b;
            }
            approx_xInteg(xInteg0, accInteg0, xInteg, dt, 1.0, b,
                          propSim->integParams.nInteg);
            std::vector<real> accIntegNext =
                get_state_der(t + dt, xInteg, propSim->forceParams,
                              propSim->integParams, propSim->consts);

            vabs_max(b[6], b6Max);
            vabs_max(accIntegNext, accIntegNextMax);
            b6TildeEstim = b6Max / accIntegNextMax;
            if (propSim->integParams.adaptiveTimestep) {
                relError = pow(b6TildeEstim / propSim->integParams.tolInteg,
                               1.0L / 7.0L);
            } else {
                relError =
                    pow(b6TildeEstim / propSim->integParams.tolInteg, 0.0L);
            }
            dtReq = dt / relError;

            if (relError <= 1 || loopCounter > maxLoops) {
                if (loopCounter > maxLoops) {
                    std::cout << "Warning: maxLoops ("<< maxLoops <<") reached without meeting "
                                 "integrator tolerance at time "
                              << t << ". Timestep: " << dt << std::endl;
                }
                oneStepDone = 1;
                propSim->integParams.timestepCounter += 1;
                loopCounter = 0;
                if (propSim->tEval.size() != propSim->xIntegEval.size()) {
                    interpolate(t, dt, xInteg0, accInteg0, b, propSim);
                }
                t += dt;
                if ((propSim->integParams.tf > propSim->integParams.t0 &&
                     t >= propSim->integParams.tf) ||
                    (propSim->integParams.tf < propSim->integParams.t0 &&
                     t <= propSim->integParams.tf)) {
                    keepStepping = 0;
                }
                b_old = b;
                refine_b(b, e, dtReq / dt, dim,
                         propSim->integParams.timestepCounter);
                check_and_apply_events(propSim, t, tNextEvent, nextEventIdx,
                                       xInteg);
                // add close approaach/impact handling here
                // xInteg0 = xInteg;
                propSim->t = t;
                propSim->xInteg = xInteg;
                propSim->tStep.push_back(t);
                propSim->xIntegStep.push_back(xInteg);
            } else {
                loopCounter += 1;
            }
            if (dtReq / dt > 1.0 / propSim->integParams.dtChangeFactor) {
                dt /= propSim->integParams.dtChangeFactor;
            } else if (fabs(dtReq) < 1.0e-12L) {
                dt *= propSim->integParams.dtChangeFactor;
            } else {
                dt = dtReq;
            }
            if ((propSim->integParams.tf > propSim->integParams.t0 &&
                 dt > propSim->integParams.dtMax) ||
                (propSim->integParams.tf < propSim->integParams.t0 &&
                 dt < propSim->integParams.dtMax)) {
                dt = propSim->integParams.dtMax;
            }
            if ((propSim->integParams.tf > propSim->integParams.t0 &&
                 dt < propSim->integParams.dtMin) ||
                (propSim->integParams.tf < propSim->integParams.t0 &&
                 dt > propSim->integParams.dtMin)) {
                dt = propSim->integParams.dtMin;
            }
            if ((propSim->integParams.tf > propSim->integParams.t0 &&
                 t + dt > tNextEvent) ||
                (propSim->integParams.tf < propSim->integParams.t0 &&
                 t + dt < tNextEvent)) {
                dt = tNextEvent - t;
            }
        }
    }
}
