#include "gr15.h"

real get_initial_timestep(const real &t, const std::vector<real> &xInteg0,
                          propSimulation *propSim) {
    real dt;
    if (propSim->integParams.dt0 != 0.0) {
        dt = fabs(propSim->integParams.dt0);
        if (propSim->integParams.tf < propSim->integParams.t0) {
            dt *= -1.0;
        }
        return dt;
    }
    int order = 15;
    real absMaxPos0, absMaxAcc0, absMaxAcc1Minus0;
    real dtTemp0, dtTemp1;
    std::vector<real> posInteg0(3 * propSim->integParams.nInteg, 0.0);
    std::vector<real> accInteg1Minus0(3 * propSim->integParams.nInteg, 0.0);
    std::vector<real> xIntegNext(6 * propSim->integParams.nInteg, 0.0);
    std::vector<real> accInteg0 =
        get_state_der(t, xInteg0, propSim);
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
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
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        for (size_t j = 0; j < 3; j++) {
            xIntegNext[6 * i + j] =
                xInteg0[6 * i + j] + dtTemp0 * xInteg0[6 * i + j + 3];
            xIntegNext[6 * i + j + 3] =
                xInteg0[6 * i + j + 3] + dtTemp0 * accInteg0[3 * i + j];
        }
    }
    std::vector<real> accInteg1 = get_state_der(
        t + dtTemp0, xIntegNext, propSim);
    vsub(accInteg1, accInteg0, accInteg1Minus0);
    vabs_max(accInteg1Minus0, absMaxAcc1Minus0);
    if (fmax(absMaxAcc0, absMaxAcc1Minus0) <= 1e-15) {
        dtTemp1 = fmax(1.0e-6, dtTemp0 * 1e-3);
    } else {
        dtTemp1 =
            pow(0.01 / fmax(absMaxAcc0, absMaxAcc1Minus0), 1.0 / (order + 1));
    }
    dt = fmin(100 * dtTemp0, dtTemp1);
    if (fabs(propSim->integParams.tf - propSim->integParams.t0) < dt) {
        dt = fabs(propSim->integParams.tf - propSim->integParams.t0);
    }
    if (propSim->integParams.tf < propSim->integParams.t0) {
        dt *= -1.0;
    }
    return dt;
}

void approx_xInteg(const std::vector<real> &xInteg0,
                   const std::vector<real> &accInteg0,
                   std::vector<real> &xIntegNext, const real &dt, const real &h,
                   const std::vector<std::vector<real>> &b,
                   const size_t &nInteg) {
    for (size_t i = 0; i < nInteg; i++) {
        xIntegNext[6*i]   = xInteg0[6*i]   + ((((((((b[6][3*i]*7.*h/9.   + b[5][3*i])*3.*h/4.   + b[4][3*i])*5.*h/7.   + b[3][3*i])*2.*h/3.   + b[2][3*i])*3.*h/5.   + b[1][3*i])*h/2.      + b[0][3*i])*h/3.   + accInteg0[3*i]  )*dt*h/2. + xInteg0[6*i+3])*dt*h;
        xIntegNext[6*i+1] = xInteg0[6*i+1] + ((((((((b[6][3*i+1]*7.*h/9. + b[5][3*i+1])*3.*h/4. + b[4][3*i+1])*5.*h/7. + b[3][3*i+1])*2.*h/3. + b[2][3*i+1])*3.*h/5. + b[1][3*i+1])*h/2.    + b[0][3*i+1])*h/3. + accInteg0[3*i+1])*dt*h/2. + xInteg0[6*i+4])*dt*h;
        xIntegNext[6*i+2] = xInteg0[6*i+2] + ((((((((b[6][3*i+2]*7.*h/9. + b[5][3*i+2])*3.*h/4. + b[4][3*i+2])*5.*h/7. + b[3][3*i+2])*2.*h/3. + b[2][3*i+2])*3.*h/5. + b[1][3*i+2])*h/2.    + b[0][3*i+2])*h/3. + accInteg0[3*i+2])*dt*h/2. + xInteg0[6*i+5])*dt*h;
        xIntegNext[6*i+3] = xInteg0[6*i+3] + (((((((b[6][3*i]*7.*h/8.    + b[5][3*i])*6.*h/7.   + b[4][3*i])*5.*h/6.   + b[3][3*i])*4.*h/5.   + b[2][3*i])*3.*h/4.   + b[1][3*i])*2.*h/3.   + b[0][3*i])*h/2.   + accInteg0[3*i]  )*dt*h;
        xIntegNext[6*i+4] = xInteg0[6*i+4] + (((((((b[6][3*i+1]*7.*h/8.  + b[5][3*i+1])*6.*h/7. + b[4][3*i+1])*5.*h/6. + b[3][3*i+1])*4.*h/5. + b[2][3*i+1])*3.*h/4. + b[1][3*i+1])*2.*h/3. + b[0][3*i+1])*h/2. + accInteg0[3*i+1])*dt*h;
        xIntegNext[6*i+5] = xInteg0[6*i+5] + (((((((b[6][3*i+2]*7.*h/8.  + b[5][3*i+2])*6.*h/7. + b[4][3*i+2])*5.*h/6. + b[3][3*i+2])*4.*h/5. + b[2][3*i+2])*3.*h/4. + b[1][3*i+2])*2.*h/3. + b[0][3*i+2])*h/2. + accInteg0[3*i+2])*dt*h;
    }
}

void compute_g_and_b(const std::vector<std::vector<real>> &AccIntegArr,
                     const size_t &hIdx, real *g,
                     std::vector<std::vector<real>> &b, const size_t &dim) {
    for (size_t i=0; i<dim; i++){
        if (hIdx == 1) {
            g[0*dim+i] = (AccIntegArr[1][i] - AccIntegArr[0][i]) * rMat[1][0];
        } else if (hIdx == 2) {
            g[1*dim+i] = ((AccIntegArr[2][i] - AccIntegArr[0][i]) * rMat[2][0] - g[0*dim+i]) * rMat[2][1];
        } else if (hIdx == 3) {
            g[2*dim+i] = (((AccIntegArr[3][i] - AccIntegArr[0][i]) * rMat[3][0] - g[0*dim+i]) * rMat[3][1] - g[1*dim+i]) * rMat[3][2];
        } else if (hIdx == 4) {
            g[3*dim+i] = ((((AccIntegArr[4][i] - AccIntegArr[0][i]) * rMat[4][0] - g[0*dim+i]) * rMat[4][1] - g[1*dim+i]) * rMat[4][2] - g[2*dim+i]) * rMat[4][3];
        } else if (hIdx == 5) {
            g[4*dim+i] = (((((AccIntegArr[5][i] - AccIntegArr[0][i]) * rMat[5][0] - g[0*dim+i]) * rMat[5][1] - g[1*dim+i]) * rMat[5][2] - g[2*dim+i]) * rMat[5][3] - g[3*dim+i]) * rMat[5][4];
        } else if (hIdx == 6) {
            g[5*dim+i] = ((((((AccIntegArr[6][i] - AccIntegArr[0][i]) * rMat[6][0] - g[0*dim+i]) * rMat[6][1] - g[1*dim+i]) * rMat[6][2] - g[2*dim+i]) * rMat[6][3] - g[3*dim+i]) * rMat[6][4] - g[4*dim+i]) * rMat[6][5];
        } else if (hIdx == 7) {
            g[6*dim+i] = (((((((AccIntegArr[7][i] - AccIntegArr[0][i]) * rMat[7][0] - g[0*dim+i]) * rMat[7][1] - g[1*dim+i]) * rMat[7][2] - g[2*dim+i]) * rMat[7][3] - g[3*dim+i]) * rMat[7][4] - g[4*dim+i]) * rMat[7][5] - g[5*dim+i]) * rMat[7][6];
        }
    }
    for (size_t i=0; i<dim; i++){
        if (hIdx == 1) {
            b[0][i] = cMat[0][0]*g[0*dim+i] + cMat[1][0]*g[1*dim+i] + cMat[2][0]*g[2*dim+i] + cMat[3][0]*g[3*dim+i] + cMat[4][0]*g[4*dim+i] + cMat[5][0]*g[5*dim+i] + cMat[6][0]*g[6*dim+i];
        } else if (hIdx == 2) {
            b[0][i] = cMat[0][0]*g[0*dim+i] + cMat[1][0]*g[1*dim+i] + cMat[2][0]*g[2*dim+i] + cMat[3][0]*g[3*dim+i] + cMat[4][0]*g[4*dim+i] + cMat[5][0]*g[5*dim+i] + cMat[6][0]*g[6*dim+i];
            b[1][i] =                       + cMat[1][1]*g[1*dim+i] + cMat[2][1]*g[2*dim+i] + cMat[3][1]*g[3*dim+i] + cMat[4][1]*g[4*dim+i] + cMat[5][1]*g[5*dim+i] + cMat[6][1]*g[6*dim+i];
        } else if (hIdx == 3) {
            b[0][i] = cMat[0][0]*g[0*dim+i] + cMat[1][0]*g[1*dim+i] + cMat[2][0]*g[2*dim+i] + cMat[3][0]*g[3*dim+i] + cMat[4][0]*g[4*dim+i] + cMat[5][0]*g[5*dim+i] + cMat[6][0]*g[6*dim+i];
            b[1][i] =                       + cMat[1][1]*g[1*dim+i] + cMat[2][1]*g[2*dim+i] + cMat[3][1]*g[3*dim+i] + cMat[4][1]*g[4*dim+i] + cMat[5][1]*g[5*dim+i] + cMat[6][1]*g[6*dim+i];
            b[2][i] =                                               + cMat[2][2]*g[2*dim+i] + cMat[3][2]*g[3*dim+i] + cMat[4][2]*g[4*dim+i] + cMat[5][2]*g[5*dim+i] + cMat[6][2]*g[6*dim+i];
        } else if (hIdx == 4) {
            b[0][i] = cMat[0][0]*g[0*dim+i] + cMat[1][0]*g[1*dim+i] + cMat[2][0]*g[2*dim+i] + cMat[3][0]*g[3*dim+i] + cMat[4][0]*g[4*dim+i] + cMat[5][0]*g[5*dim+i] + cMat[6][0]*g[6*dim+i];
            b[1][i] =                       + cMat[1][1]*g[1*dim+i] + cMat[2][1]*g[2*dim+i] + cMat[3][1]*g[3*dim+i] + cMat[4][1]*g[4*dim+i] + cMat[5][1]*g[5*dim+i] + cMat[6][1]*g[6*dim+i];
            b[2][i] =                                               + cMat[2][2]*g[2*dim+i] + cMat[3][2]*g[3*dim+i] + cMat[4][2]*g[4*dim+i] + cMat[5][2]*g[5*dim+i] + cMat[6][2]*g[6*dim+i];
            b[3][i] =                                                                       + cMat[3][3]*g[3*dim+i] + cMat[4][3]*g[4*dim+i] + cMat[5][3]*g[5*dim+i] + cMat[6][3]*g[6*dim+i];
        } else if (hIdx == 5) {
            b[0][i] = cMat[0][0]*g[0*dim+i] + cMat[1][0]*g[1*dim+i] + cMat[2][0]*g[2*dim+i] + cMat[3][0]*g[3*dim+i] + cMat[4][0]*g[4*dim+i] + cMat[5][0]*g[5*dim+i] + cMat[6][0]*g[6*dim+i];
            b[1][i] =                       + cMat[1][1]*g[1*dim+i] + cMat[2][1]*g[2*dim+i] + cMat[3][1]*g[3*dim+i] + cMat[4][1]*g[4*dim+i] + cMat[5][1]*g[5*dim+i] + cMat[6][1]*g[6*dim+i];
            b[2][i] =                                               + cMat[2][2]*g[2*dim+i] + cMat[3][2]*g[3*dim+i] + cMat[4][2]*g[4*dim+i] + cMat[5][2]*g[5*dim+i] + cMat[6][2]*g[6*dim+i];
            b[3][i] =                                                                       + cMat[3][3]*g[3*dim+i] + cMat[4][3]*g[4*dim+i] + cMat[5][3]*g[5*dim+i] + cMat[6][3]*g[6*dim+i];
            b[4][i] =                                                                                               + cMat[4][4]*g[4*dim+i] + cMat[5][4]*g[5*dim+i] + cMat[6][4]*g[6*dim+i];
        } else if (hIdx == 6) {
            b[0][i] = cMat[0][0]*g[0*dim+i] + cMat[1][0]*g[1*dim+i] + cMat[2][0]*g[2*dim+i] + cMat[3][0]*g[3*dim+i] + cMat[4][0]*g[4*dim+i] + cMat[5][0]*g[5*dim+i] + cMat[6][0]*g[6*dim+i];
            b[1][i] =                       + cMat[1][1]*g[1*dim+i] + cMat[2][1]*g[2*dim+i] + cMat[3][1]*g[3*dim+i] + cMat[4][1]*g[4*dim+i] + cMat[5][1]*g[5*dim+i] + cMat[6][1]*g[6*dim+i];
            b[2][i] =                                               + cMat[2][2]*g[2*dim+i] + cMat[3][2]*g[3*dim+i] + cMat[4][2]*g[4*dim+i] + cMat[5][2]*g[5*dim+i] + cMat[6][2]*g[6*dim+i];
            b[3][i] =                                                                       + cMat[3][3]*g[3*dim+i] + cMat[4][3]*g[4*dim+i] + cMat[5][3]*g[5*dim+i] + cMat[6][3]*g[6*dim+i];
            b[4][i] =                                                                                               + cMat[4][4]*g[4*dim+i] + cMat[5][4]*g[5*dim+i] + cMat[6][4]*g[6*dim+i];
            b[5][i] =                                                                                                                       + cMat[5][5]*g[5*dim+i] + cMat[6][5]*g[6*dim+i];
        } else if (hIdx == 7) {
            b[0][i] = cMat[0][0]*g[0*dim+i] + cMat[1][0]*g[1*dim+i] + cMat[2][0]*g[2*dim+i] + cMat[3][0]*g[3*dim+i] + cMat[4][0]*g[4*dim+i] + cMat[5][0]*g[5*dim+i] + cMat[6][0]*g[6*dim+i];
            b[1][i] =                       + cMat[1][1]*g[1*dim+i] + cMat[2][1]*g[2*dim+i] + cMat[3][1]*g[3*dim+i] + cMat[4][1]*g[4*dim+i] + cMat[5][1]*g[5*dim+i] + cMat[6][1]*g[6*dim+i];
            b[2][i] =                                               + cMat[2][2]*g[2*dim+i] + cMat[3][2]*g[3*dim+i] + cMat[4][2]*g[4*dim+i] + cMat[5][2]*g[5*dim+i] + cMat[6][2]*g[6*dim+i];
            b[3][i] =                                                                       + cMat[3][3]*g[3*dim+i] + cMat[4][3]*g[4*dim+i] + cMat[5][3]*g[5*dim+i] + cMat[6][3]*g[6*dim+i];
            b[4][i] =                                                                                               + cMat[4][4]*g[4*dim+i] + cMat[5][4]*g[5*dim+i] + cMat[6][4]*g[6*dim+i];
            b[5][i] =                                                                                                                       + cMat[5][5]*g[5*dim+i] + cMat[6][5]*g[6*dim+i];
            b[6][i] =                                                                                                                                               + cMat[6][6]*g[6*dim+i];
        }
    }
}

void refine_b(std::vector<std::vector<real>> &b,
              real *e, const real &dtRatio,
              const size_t &dim) {
    std::vector<std::vector<real>> bDiff(7, std::vector<real>(dim, 0.0L));
    for (size_t i = 0; i < dim; i++){
        bDiff[0][i] = b[0][i] - e[0*dim+i];
        bDiff[1][i] = b[1][i] - e[1*dim+i];
        bDiff[2][i] = b[2][i] - e[2*dim+i];
        bDiff[3][i] = b[3][i] - e[3*dim+i];
        bDiff[4][i] = b[4][i] - e[4*dim+i];
        bDiff[5][i] = b[5][i] - e[5*dim+i];
        bDiff[6][i] = b[6][i] - e[6*dim+i];
    }

    real q = dtRatio;
    real q2 = q * q;
    real q3 = q2 * q;
    real q4 = q2 * q2;
    real q5 = q2 * q3;
    real q6 = q3 * q3;
    real q7 = q2 * q5;

    for (size_t i = 0; i < dim; i++) {
        e[0*dim+i] = q  * (b[6][i] * 7.0  + b[5][i] * 6.0  + b[4][i] * 5.0  + b[3][i] * 4.0 + b[2][i] * 3.0 + b[1][i] * 2.0 + b[0][i]);
        e[1*dim+i] = q2 * (b[6][i] * 21.0 + b[5][i] * 15.0 + b[4][i] * 10.0 + b[3][i] * 6.0 + b[2][i] * 3.0 + b[1][i]);
        e[2*dim+i] = q3 * (b[6][i] * 35.0 + b[5][i] * 20.0 + b[4][i] * 10.0 + b[3][i] * 4.0 + b[2][i]);
        e[3*dim+i] = q4 * (b[6][i] * 35.0 + b[5][i] * 15.0 + b[4][i] * 5.0  + b[3][i]);
        e[4*dim+i] = q5 * (b[6][i] * 21.0 + b[5][i] * 6.0  + b[4][i]);
        e[5*dim+i] = q6 * (b[6][i] * 7.0  + b[5][i]);
        e[6*dim+i] = q7 * (b[6][i]);
    }

    for (size_t i = 0; i < dim; i++) {
        b[0][i] = e[0*dim+i] + bDiff[0][i];
        b[1][i] = e[1*dim+i] + bDiff[1][i];
        b[2][i] = e[2*dim+i] + bDiff[2][i];
        b[3][i] = e[3*dim+i] + bDiff[3][i];
        b[4][i] = e[4*dim+i] + bDiff[4][i];
        b[5][i] = e[5*dim+i] + bDiff[5][i];
        b[6][i] = e[6*dim+i] + bDiff[6][i];
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

static real root7(real num){
    real fac = 1.0;
    if (num<1e-7){
        fac = 0.1;
        num *= 1e7;
    }
    if (num>1e2){
        fac = 10.0;
        num *= 1e-7;
    }
    real root = 1.0;
    for (size_t i = 0; i < 20; i++) {
        root += (num/(root*root*root*root*root*root)-root)/7.0;
    }
    return fac*root;
}

void gr15(propSimulation *propSim) {
    real t = propSim->t;
    std::vector<real> xInteg0 = propSim->xInteg;
    size_t nh = 8;
    size_t dim = propSim->integParams.n2Derivs;
    real dt = get_initial_timestep(t, xInteg0, propSim);
    propSim->integParams.timestepCounter = 0;
    std::vector<real> accInteg0 = get_state_der(t, xInteg0, propSim);
    std::vector<real> accIntegNext = std::vector<real>(accInteg0.size(), 0.0);
    std::vector<real> xInteg(2 * dim, 0.0);
    std::vector<std::vector<real> > bOld(7, std::vector<real>(dim, 0.0));
    std::vector<std::vector<real> > b(7, std::vector<real>(dim, 0.0));
    real *g = new real[7 * dim];
    memset(g, 0.0, 7 * dim * sizeof(real));
    real *e = new real[7 * dim];
    memset(e, 0.0, 7 * dim * sizeof(real));
    std::vector<std::vector<real> > accIntegArr(nh,
                                                std::vector<real>(dim, 0.0));
    real *b6Tilde = new real[dim];
    memset(b6Tilde, 0.0, dim * sizeof(real));
    real b6TildeMax, accIntegArr7Max;
    real b6Max, accIntegNextMax;
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
            real PCerr = 1.0/propSim->integParams.tolPC;
            real PCerrPrev = PCerr + 1.0;
            accIntegArr[0] = accInteg0;
            size_t PCIter = 0;
            while (true) {
                if (PCerr < propSim->integParams.tolPC) {
                    break;
                }
                if (PCerr > PCerrPrev && PCIter > 2) {
                    break;
                }
                if (PCIter > PCmaxIter) {
                    break;
                }
                PCerrPrev = PCerr;
                PCIter++;
                for (size_t hIdx = 1; hIdx < nh; hIdx++) {
                    approx_xInteg(xInteg0, accInteg0, xInteg, dt, hVec[hIdx], b,
                                  propSim->integParams.nInteg);
                    accIntegArr[hIdx] = get_state_der(
                        t + hVec[hIdx] * dt, xInteg, propSim);
                    compute_g_and_b(accIntegArr, hIdx, g, b, dim);
                    if (hIdx == 7) {
                        for (size_t i = 0; i < dim; i++) {
                            b6Tilde[i] = b[6][i] - bOld[6][i];
                        }
                        vabs_max(b6Tilde, dim, b6TildeMax);
                        vabs_max(accIntegArr[7], accIntegArr7Max);
                        PCerr = b6TildeMax / accIntegArr7Max;
                        bOld = b;
                    }
                }
            }
            // printf("ks: %d. osd: %d. t: %0.10e, dt: %0.8e", keepStepping, oneStepDone,t, dt);
            // printf(". iterations: %d",PCIter);
            // printf(". predictor_corrector_error: %e",PCerr);
            vabs_max(b[6], b6Max);
            vabs_max(accIntegArr[7], accIntegNextMax);
            relError = b6Max / accIntegNextMax;
            if (propSim->integParams.adaptiveTimestep) {
                if (std::isnormal(relError)) {
                    dtReq = root7(propSim->integParams.tolInteg/relError)*dt;
                } else {
                    dtReq = dt/propSim->integParams.dtChangeFactor;
                }
            } else {
                dtReq = dt;
            }
            // printf(". dtReq: %e. relError: %e", dtReq,relError);
            if (fabs(dtReq) < propSim->integParams.dtMin) {
                dtReq = copysign(propSim->integParams.dtMin, dtReq);
            }
            if (fabs(dtReq/dt) < propSim->integParams.dtChangeFactor) {
                if (propSim->integParams.timestepCounter > 1) {
                    refine_b(b, e, dtReq / dt, dim);
                }
                dt = dtReq;
                std::cout << "Restarting next while loop at time t = " << t << std::endl;
                continue;
            }
            if (dtReq / dt > 1.0 / propSim->integParams.dtChangeFactor) {
                dtReq = dt / propSim->integParams.dtChangeFactor;
            }
            // final cleanup
            // printf(". cleaning up...\n");
            propSim->interpParams.bStack.push_back(bOld);
            propSim->interpParams.accIntegStack.push_back(accInteg0);
            if (propSim->tEval.size() != propSim->xIntegEval.size()) {
                // interpolate(t, dt, xInteg0, accInteg0, b, propSim);
                interpolate_on_the_fly(propSim, t, dt);
            }
            approx_xInteg(xInteg0, accInteg0, xInteg, dt, 1.0, b,
                        propSim->integParams.nInteg);
            accInteg0 = get_state_der(t + dt, xInteg, propSim);
            t += dt;
            propSim->t = t;
            check_and_apply_events(propSim, t, tNextEvent, nextEventIdx,
                                    xInteg);
            propSim->xInteg = xInteg;
            propSim->interpParams.tStack.push_back(t);
            propSim->interpParams.xIntegStack.push_back(xInteg);
            propSim->integParams.timestepCounter++;
            if (propSim->integParams.timestepCounter > 1) {
                refine_b(b, e, dtReq / dt, dim);
            }
            check_ca_or_impact(propSim, t-dt, xInteg0, t, xInteg, keepStepping);
            if ((propSim->integParams.tf > propSim->integParams.t0 &&
                    t >= propSim->integParams.tf) ||
                (propSim->integParams.tf < propSim->integParams.t0 &&
                    t <= propSim->integParams.tf)) {
                keepStepping = 0;
            }
            dt = dtReq;
            if ((propSim->integParams.tf > propSim->integParams.t0 &&
                 t + dt > tNextEvent) ||
                (propSim->integParams.tf < propSim->integParams.t0 &&
                 t + dt < tNextEvent)) {
                dt = tNextEvent - t;
            }
            oneStepDone = 1;
        }
    }
    delete[] g;
    delete[] e;
    delete[] b6Tilde;
}
