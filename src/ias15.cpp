/**
 * @file    ias15.cpp
 * @brief   IAS15 integrator.
 * @author  Rahil Makadia <makadia2@illinois.edu>
 * @details This file implements the integration scheme for GRSS.
 * The integratior used here is based mostly on the IAS15 integrator published
 * by Rein and Spiegel in 2014. Current and previous versions of the integrator
 * are a hybrid amalgamation of the work done in the following sources:
 * 1. The original RADAU code by Everhart (1985)
 * 2. The IAS15 modifications by Rein and Spiegel (2014)
 * 3. A version of IAS15 used in the book Moving Planets Around (2020)
 * 4. The new timestepping criterion introduced by Pham et al. (2024)
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

#include "ias15.h"

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @return real Initial timestep.
 */
real get_initial_timestep(PropSimulation *propSim){
    real dt0 = 1.0L;
    if (propSim->integParams.dt0 != 0.0) {
        dt0 = fabs(propSim->integParams.dt0);
    }
    real span = fabs(propSim->integParams.tf - propSim->integParams.t0);
    if (span < dt0) {
        dt0 = span;
    }
    if (propSim->integParams.tf < propSim->integParams.t0) {
        dt0 *= -1.0;
    }
    return dt0;
}

/**
 * @param[in] b Matrix of interpolation coefficients.
 * @param[in] dim Dimension of the system (number of 2nd derivatives).
 * @param[out] g Matrix of interpolation coefficients.
 */
void update_g_with_b(const std::vector<real> &b, const size_t &dim, std::vector<real> &g) {
    for (size_t i = 0; i < dim; i++) {
        g[0*dim+i] = b[6*dim+i]*dVec[15] + b[5*dim+i]*dVec[10] + b[4*dim+i]*dVec[6] + b[3*dim+i]*dVec[3] + b[2*dim+i]*dVec[1] + b[1*dim+i]*dVec[0] + b[0*dim+i];
        g[1*dim+i] = b[6*dim+i]*dVec[16] + b[5*dim+i]*dVec[11] + b[4*dim+i]*dVec[7] + b[3*dim+i]*dVec[4] + b[2*dim+i]*dVec[2] + b[1*dim+i];
        g[2*dim+i] = b[6*dim+i]*dVec[17] + b[5*dim+i]*dVec[12] + b[4*dim+i]*dVec[8] + b[3*dim+i]*dVec[5] + b[2*dim+i];
        g[3*dim+i] = b[6*dim+i]*dVec[18] + b[5*dim+i]*dVec[13] + b[4*dim+i]*dVec[9] + b[3*dim+i];
        g[4*dim+i] = b[6*dim+i]*dVec[19] + b[5*dim+i]*dVec[14] + b[4*dim+i];
        g[5*dim+i] = b[6*dim+i]*dVec[20] + b[5*dim+i];
        g[6*dim+i] = b[6*dim+i];
    }
}

/**
 * @param[in] AccIntegArr Array of acceleration values at the nodes.
 * @param[in] hIdx Gauss-Radau node index.
 * @param[out] g Matrix of interpolation coefficients.
 * @param[out] bCompCoeffs Vector of compensated summation coefficients for computing the b-matrix.
 * @param[out] b Matrix of interpolation coefficients.
 * @param[in] dim Dimension of the system (number of 2nd derivatives).
 * @param[out] PCerr Error in the predictor-corrector step (only for hIdx = 7)
 */
void compute_g_and_b(const std::vector<std::vector<real>> &AccIntegArr,
                     const size_t &hIdx, std::vector<real> &g, std::vector<real> &bCompCoeffs,
                     std::vector<real> &b, const size_t &dim, real &PCerr) {
    switch (hIdx) {
        case 0:
        {
            throw std::runtime_error("hIdx cannot be 0");
        }
        case 1:
        {
            for (size_t i = 0; i < dim; i++) {
                real dummy = 0;
                real temp = g[0*dim+i];
                real gVal = AccIntegArr[hIdx][i];
                comp_sum(-AccIntegArr[0][i], &gVal, &dummy);
                g[0*dim+i] = gVal/rVec[0];
                comp_sum(g[0*dim+i]-temp, &(b[0*dim+i]), &(bCompCoeffs[0*dim+i]));
            }
            break;
        }
        case 2:
        {
            for (size_t i = 0; i < dim; i++) {
                real dummy = 0;
                real temp = g[1*dim+i];
                real gVal = AccIntegArr[hIdx][i];
                comp_sum(-AccIntegArr[0][i], &gVal, &dummy);
                g[1*dim+i] = (gVal/rVec[1] - g[0*dim+i])/rVec[2];
                temp = g[1*dim+i] - temp;
                comp_sum(temp*cVec[0], &(b[0*dim+i]), &(bCompCoeffs[0*dim+i]));
                comp_sum(temp        , &(b[1*dim+i]), &(bCompCoeffs[1*dim+i]));
            }
            break;
        }
        case 3:
        {
            for (size_t i = 0; i < dim; i++) {
                real dummy = 0;
                real temp = g[2*dim+i];
                real gVal = AccIntegArr[hIdx][i];
                comp_sum(-AccIntegArr[0][i], &gVal, &dummy);
                g[2*dim+i] = ((gVal/rVec[3] - g[0*dim+i])/rVec[4] - g[1*dim+i])/rVec[5];
                temp = g[2*dim+i] - temp;
                comp_sum(temp*cVec[1], &(b[0*dim+i]), &(bCompCoeffs[0*dim+i]));
                comp_sum(temp*cVec[2], &(b[1*dim+i]), &(bCompCoeffs[1*dim+i]));
                comp_sum(temp        , &(b[2*dim+i]), &(bCompCoeffs[2*dim+i]));
            }
            break;
        }
        case 4:
        {
            for (size_t i = 0; i < dim; i++) {
                real dummy = 0;
                real temp = g[3*dim+i];
                real gVal = AccIntegArr[hIdx][i];
                comp_sum(-AccIntegArr[0][i], &gVal, &dummy);
                g[3*dim+i] = (((gVal/rVec[6] - g[0*dim+i])/rVec[7] - g[1*dim+i])/rVec[8] - g[2*dim+i])/rVec[9];
                temp = g[3*dim+i] - temp;
                comp_sum(temp*cVec[3], &(b[0*dim+i]), &(bCompCoeffs[0*dim+i]));
                comp_sum(temp*cVec[4], &(b[1*dim+i]), &(bCompCoeffs[1*dim+i]));
                comp_sum(temp*cVec[5], &(b[2*dim+i]), &(bCompCoeffs[2*dim+i]));
                comp_sum(temp        , &(b[3*dim+i]), &(bCompCoeffs[3*dim+i]));
            }
            break;
        }
        case 5:
        {
            for (size_t i = 0; i < dim; i++) {
                real dummy = 0;
                real temp = g[4*dim+i];
                real gVal = AccIntegArr[hIdx][i];
                comp_sum(-AccIntegArr[0][i], &gVal, &dummy);
                g[4*dim+i] = ((((gVal/rVec[10] - g[0*dim+i])/rVec[11] - g[1*dim+i])/rVec[12] - g[2*dim+i])/rVec[13] - g[3*dim+i])/rVec[14];
                temp = g[4*dim+i] - temp;
                comp_sum(temp*cVec[6], &(b[0*dim+i]), &(bCompCoeffs[0*dim+i]));
                comp_sum(temp*cVec[7], &(b[1*dim+i]), &(bCompCoeffs[1*dim+i]));
                comp_sum(temp*cVec[8], &(b[2*dim+i]), &(bCompCoeffs[2*dim+i]));
                comp_sum(temp*cVec[9], &(b[3*dim+i]), &(bCompCoeffs[3*dim+i]));
                comp_sum(temp        , &(b[4*dim+i]), &(bCompCoeffs[4*dim+i]));
            }
            break;
        }
        case 6:
        {
            for (size_t i = 0; i < dim; i++) {
                real dummy = 0;
                real temp = g[5*dim+i];
                real gVal = AccIntegArr[hIdx][i];
                comp_sum(-AccIntegArr[0][i], &gVal, &dummy);
                g[5*dim+i] = (((((gVal/rVec[15] - g[0*dim+i])/rVec[16] - g[1*dim+i])/rVec[17] - g[2*dim+i])/rVec[18] - g[3*dim+i])/rVec[19] - g[4*dim+i])/rVec[20];
                temp = g[5*dim+i] - temp;
                comp_sum(temp*cVec[10], &(b[0*dim+i]), &(bCompCoeffs[0*dim+i]));
                comp_sum(temp*cVec[11], &(b[1*dim+i]), &(bCompCoeffs[1*dim+i]));
                comp_sum(temp*cVec[12], &(b[2*dim+i]), &(bCompCoeffs[2*dim+i]));
                comp_sum(temp*cVec[13], &(b[3*dim+i]), &(bCompCoeffs[3*dim+i]));
                comp_sum(temp*cVec[14], &(b[4*dim+i]), &(bCompCoeffs[4*dim+i]));
                comp_sum(temp         , &(b[5*dim+i]), &(bCompCoeffs[5*dim+i]));
            }
            break;
        }
        case 7:
        {
            real maxAccI = 0;
            real maxB6Tilde = 0;
            for (size_t i = 0; i < dim; i++) {
                real dummy = 0;
                real temp = g[6*dim+i];
                real gVal = AccIntegArr[hIdx][i];
                comp_sum(-AccIntegArr[0][i], &gVal, &dummy);
                g[6*dim+i] = ((((((gVal/rVec[21] - g[0*dim+i])/rVec[22] - g[1*dim+i])/rVec[23] - g[2*dim+i])/rVec[24] - g[3*dim+i])/rVec[25] - g[4*dim+i])/rVec[26] - g[5*dim+i])/rVec[27];
                temp = g[6*dim+i] - temp;
                comp_sum(temp*cVec[15], &(b[0*dim+i]), &(bCompCoeffs[0*dim+i]));
                comp_sum(temp*cVec[16], &(b[1*dim+i]), &(bCompCoeffs[1*dim+i]));
                comp_sum(temp*cVec[17], &(b[2*dim+i]), &(bCompCoeffs[2*dim+i]));
                comp_sum(temp*cVec[18], &(b[3*dim+i]), &(bCompCoeffs[3*dim+i]));
                comp_sum(temp*cVec[19], &(b[4*dim+i]), &(bCompCoeffs[4*dim+i]));
                comp_sum(temp*cVec[20], &(b[5*dim+i]), &(bCompCoeffs[5*dim+i]));
                comp_sum(temp         , &(b[6*dim+i]), &(bCompCoeffs[6*dim+i]));
                const real accI = fabs(AccIntegArr[hIdx][i]);
                if (std::isnormal(accI) && accI > maxAccI) {
                    maxAccI = accI;
                }
                const real b6Tilde = fabs(temp);
                if (std::isnormal(b6Tilde) && b6Tilde > maxB6Tilde) {
                    maxB6Tilde = b6Tilde;
                }
            }
            PCerr = maxB6Tilde / maxAccI;
            break;
        }
        default:
            throw std::runtime_error("hIdx is out of range");
    }
}

/**
 * @param[inout] b Matrix of interpolation coefficients.
 * @param[inout] e Refinement matrix.
 * @param[in] dtRatio Ratio of the required timestep to the previous timestep.
 * @param[in] dim Dimension of the system (number of 2nd derivatives).
 */
void refine_b(std::vector<real> &b, std::vector<real> &e, const real &dtRatio, const size_t &dim) {
    const real q = dtRatio;
    const real q2 = q * q;
    const real q3 = q2 * q;
    const real q4 = q2 * q2;
    const real q5 = q2 * q3;
    const real q6 = q3 * q3;
    const real q7 = q2 * q5;

    for (size_t i = 0; i < dim; i++){
        const real bDiff0 = b[0*dim+i] - e[0*dim+i];
        const real bDiff1 = b[1*dim+i] - e[1*dim+i];
        const real bDiff2 = b[2*dim+i] - e[2*dim+i];
        const real bDiff3 = b[3*dim+i] - e[3*dim+i];
        const real bDiff4 = b[4*dim+i] - e[4*dim+i];
        const real bDiff5 = b[5*dim+i] - e[5*dim+i];
        const real bDiff6 = b[6*dim+i] - e[6*dim+i];

        e[0*dim+i] = q  * (b[6*dim+i] * 7.0  + b[5*dim+i] * 6.0  + b[4*dim+i] * 5.0  + b[3*dim+i] * 4.0 + b[2*dim+i] * 3.0 + b[1*dim+i] * 2.0 + b[0*dim+i]);
        e[1*dim+i] = q2 * (b[6*dim+i] * 21.0 + b[5*dim+i] * 15.0 + b[4*dim+i] * 10.0 + b[3*dim+i] * 6.0 + b[2*dim+i] * 3.0 + b[1*dim+i]);
        e[2*dim+i] = q3 * (b[6*dim+i] * 35.0 + b[5*dim+i] * 20.0 + b[4*dim+i] * 10.0 + b[3*dim+i] * 4.0 + b[2*dim+i]);
        e[3*dim+i] = q4 * (b[6*dim+i] * 35.0 + b[5*dim+i] * 15.0 + b[4*dim+i] * 5.0  + b[3*dim+i]);
        e[4*dim+i] = q5 * (b[6*dim+i] * 21.0 + b[5*dim+i] * 6.0  + b[4*dim+i]);
        e[5*dim+i] = q6 * (b[6*dim+i] * 7.0  + b[5*dim+i]);
        e[6*dim+i] = q7 * (b[6*dim+i]);

        b[0*dim+i] = e[0*dim+i] + bDiff0;
        b[1*dim+i] = e[1*dim+i] + bDiff1;
        b[2*dim+i] = e[2*dim+i] + bDiff2;
        b[3*dim+i] = e[3*dim+i] + bDiff3;
        b[4*dim+i] = e[4*dim+i] + bDiff4;
        b[5*dim+i] = e[5*dim+i] + bDiff5;
        b[6*dim+i] = e[6*dim+i] + bDiff6;
    }
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] t Current time.
 * @param[in] xInteg Vector of integrated states after the current timestep.
 */
void check_and_apply_impulsive_events(PropSimulation *propSim, const real &t,
                                      std::vector<real> &xInteg) {
    size_t *nextImpEventIdx = &propSim->eventMngr.nextImpEventIdx;
    real *tNextImpEvent = &propSim->eventMngr.tNextImpEvent;
    while (*nextImpEventIdx < propSim->eventMngr.nImpEvents && t == *tNextImpEvent) {
        // apply events for the state just reached by the integrator
        propSim->eventMngr.impulsiveEvents[*nextImpEventIdx].apply_impulsive(propSim, t, xInteg);
        // update next event index and time
        (*nextImpEventIdx)++;
        if (*nextImpEventIdx < propSim->eventMngr.nImpEvents) {
            *tNextImpEvent = propSim->eventMngr.impulsiveEvents[*nextImpEventIdx].t;
        } else {
            *tNextImpEvent = propSim->integParams.tf;
        }
    }
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] t Current time.
 */
void check_continuous_events(PropSimulation *propSim, const real &t) {
    if (!propSim->eventMngr.allConEventDone){
        bool allDone = true;
        bool forwardProp = propSim->integParams.tf > propSim->integParams.t0;
        for (size_t i = 0; i < propSim->eventMngr.nConEvents; i++) {
            bool eventStarted, eventEnded;
            if (forwardProp) {
                eventStarted = t >= propSim->eventMngr.continuousEvents[i].t;
                eventEnded = t >= propSim->integParams.tf;
            } else {
                eventStarted = t <= propSim->integParams.t0;
                eventEnded = t < propSim->eventMngr.continuousEvents[i].t;
            }
            if (eventStarted && !eventEnded) {
                propSim->eventMngr.continuousEvents[i].hasStarted = true;
            } else {
                propSim->eventMngr.continuousEvents[i].hasStarted = false;
            }
            if (!eventEnded) {
                allDone = false;
            }
        }
        if (allDone) {
            propSim->eventMngr.allConEventDone = true;
            propSim->eventMngr.tNextConEvent = propSim->integParams.tf;
        }
    }
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] t Current time.
 * @param[in] xInteg Vector of integrated states after the current timestep.
 */
void check_events(PropSimulation *propSim, const real &t, std::vector<real> &xInteg){
    check_and_apply_impulsive_events(propSim, t, xInteg);
    check_continuous_events(propSim, t);
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] dt Current timestep.
 */
void event_timestep_check(PropSimulation *propSim, real &dt) {
    real tNextEvent = propSim->integParams.tf;
    const bool forwardProp = propSim->integParams.tf > propSim->integParams.t0;
    if (propSim->eventMngr.nImpEvents > 0) {
        const real tNextImpEvent = propSim->eventMngr.tNextImpEvent;
        if (forwardProp) {
            tNextEvent = fmin(tNextEvent, tNextImpEvent);
        } else {
            tNextEvent = fmax(tNextEvent, tNextImpEvent);
        }
    }
    if (!propSim->eventMngr.allConEventDone) {
        for (size_t i = 0; i < propSim->eventMngr.nConEvents; i++) { 
            const real tNextConEvent = propSim->eventMngr.continuousEvents[i].t;
            if (forwardProp && propSim->t < tNextConEvent) {
                tNextEvent = fmin(tNextEvent, tNextConEvent);
            } else if (!forwardProp && propSim->t > tNextConEvent) {
                tNextEvent = fmax(tNextEvent, tNextConEvent);
            }
            // break up the timestep if early in a continuous event
            if (propSim->eventMngr.continuousEvents[i].hasStarted) {
                const real timeConstant = propSim->eventMngr.continuousEvents[i].tau;
                const real timeConstantFac = 5.0L;
                const real numSegments = 100.0L;
                const real dtSegment = timeConstantFac*timeConstant/numSegments;
                if (forwardProp && propSim->t < tNextConEvent+timeConstantFac*timeConstant) {
                    dt = fmin(dt, dtSegment);
                } else if (!forwardProp && propSim->t < tNextConEvent+timeConstantFac*timeConstant) {
                    dt = fmax(dt, -dtSegment);
                }
            }
        }
    }
    if ((forwardProp && propSim->t + dt > tNextEvent) ||
        (!forwardProp && propSim->t + dt < tNextEvent)) {
        dt = tNextEvent - propSim->t;
    }
}

/**
 * @brief Compute the 7th root of a number.
 * 
 * @param[in] num Number to compute the 7th root of.
 * @return real 7th root of the number.
 */
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

/**
 * @brief Get the next timestep based on the new IAS15 adaptive timestep criterion.
 * 
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] dt Current timestep.
 * @param[in] accInteg0 Vector of acceleration values at the current timestep.
 * @param[in] dim Dimension of the system (number of 2nd derivatives).
 * @param[in] b Matrix of interpolation coefficients.
 * @return real Next timestep.
 */
static real get_adaptive_timestep(PropSimulation *propSim, const real &dt,
                                  const std::vector<real> &accInteg0,
                                  const size_t &dim, const std::vector<real> &b) {
    real minTimescale2 = 1e300L;
    size_t startb = 0;
    for (size_t i = 0; i < propSim->integParams.nInteg; i++){
        real y2 = 0.0L;
        real y3 = 0.0L;
        real y4 = 0.0L;
        for (size_t j = 0; j < 3; j++){
            const size_t k = startb+j;
            real temp = accInteg0[k] + b[0*dim+k] + b[1*dim+k] + b[2*dim+k] + b[3*dim+k] + b[4*dim+k] + b[5*dim+k] + b[6*dim+k];
            y2 += temp*temp;
            temp = b[0*dim+k] + 2.0L*b[1*dim+k] + 3.0L*b[2*dim+k] + 4.0L*b[3*dim+k] + 5.0L*b[4*dim+k] + 6.0L*b[5*dim+k] + 7.0L*b[6*dim+k];
            y3 += temp*temp;
            temp = 2.0L*b[1*dim+k] + 6.0L*b[2*dim+k] + 12.0L*b[3*dim+k] + 20.0L*b[4*dim+k] + 30.0L*b[5*dim+k] + 42.0L*b[6*dim+k];
            y4 += temp*temp;
        }
        const real timescale2 = 2.0L*y2/(y3+sqrt(y4*y2));
        if (std::isnormal(timescale2) && timescale2 < minTimescale2){
            minTimescale2 = timescale2;
        }
        startb += propSim->integBodies[i].n2Derivs;
    }
    real dtReq;
    if (std::isnormal(minTimescale2)){
        // Numerical factor below is for relating the epsilon in original IAS15 to the 
        // eta in the new IAS15 stepping from https://arxiv.org/pdf/2401.02849.pdf
        dtReq = sqrt(minTimescale2) * root7(propSim->integParams.tolInteg*5040.0) * dt;
    }else{
        dtReq = dt/propSim->integParams.dtChangeFactor;
    }
    return dtReq;
}

/**
 * @param[inout] propSim PropSimulation object for the integration.
 */
void ias15(PropSimulation *propSim) {
    real t = propSim->t;
    std::vector<real> xInteg = propSim->xInteg;
    if (!std::isfinite(t)) {
        throw std::runtime_error("t is not finite");
    }
    for (size_t i = 0; i < xInteg.size(); i++) {
        if (!std::isfinite(xInteg[i])) {
            throw std::runtime_error("xInteg is not finite");
        }
    }
    const size_t nh = 8;
    const size_t dim = propSim->integParams.n2Derivs;
    real dt = get_initial_timestep(propSim);
    propSim->integParams.timestepCounter = 0;
    std::vector<real> b(7 * dim, 0.0);
    std::vector<real> bCompCoeffs(7 * dim, 0.0);
    std::vector<real> g(7 * dim, 0.0);
    std::vector<real> e(7 * dim, 0.0);
    real dtReq;
    check_events(propSim, t, xInteg);
    event_timestep_check(propSim, dt);
    propSim->interpParams.tStack.push_back(t);
    propSim->interpParams.xIntegStack.push_back(xInteg);
    std::vector<real> accInteg0(dim, 0.0);
    get_state_der(propSim, t, xInteg, accInteg0);
    std::vector<std::vector<real>> accIntegArr(nh, std::vector<real>(dim, 0.0));
    size_t PCmaxIter = 12;
    int keepStepping = 1;
    int oneStepDone = 0;
    if (propSim->integParams.t0 == propSim->integParams.tf) {
        keepStepping = 0;
    }
    std::vector<real> xInteg0(xInteg.size(), 0.0);
    std::vector<real> xIntegCompCoeffs(xInteg.size(), 0.0);
    while (keepStepping) {
        xInteg0 = xInteg;
        oneStepDone = 0;
        update_g_with_b(b, dim, g);
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
                    approx_xInteg(xInteg0, accInteg0, dt, hVec[hIdx], b, dim,
                                  propSim->integBodies, xInteg, xIntegCompCoeffs);
                    get_state_der(propSim, t + hVec[hIdx] * dt, xInteg, accIntegArr[hIdx]);
                    compute_g_and_b(accIntegArr, hIdx, g, bCompCoeffs, b, dim, PCerr);
                }
            }
            if (propSim->integParams.adaptiveTimestep) {
                dtReq = get_adaptive_timestep(propSim, dt, accInteg0, dim, b);
            } else {
                dtReq = dt;
            }
            if (fabs(dtReq) < propSim->integParams.dtMin) {
                dtReq = copysign(propSim->integParams.dtMin, dtReq);
            }
            if (fabs(dtReq/dt) < propSim->integParams.dtChangeFactor) {
                // std::cout << "Restarting next while loop at time t = " << t
                //           << " with dt = " << dtReq << " instead of " << dt
                //           << std::endl;
                dt = dtReq;
                continue;
            }
            // accept step
            if (dtReq / dt > 1.0 / propSim->integParams.dtChangeFactor) {
                dtReq = dt / propSim->integParams.dtChangeFactor;
            }
            propSim->interpParams.t0 = t;
            propSim->interpParams.dt0 = dt;
            propSim->interpParams.xInteg0 = xInteg0;
            propSim->interpParams.accInteg0 = accInteg0;
            propSim->interpParams.b0 = b;
            propSim->interpParams.bStack.push_back(b);
            propSim->interpParams.accIntegStack.push_back(accInteg0);
            if (propSim->tEval.size() != propSim->xIntegEval.size()) {
                interpolate_on_the_fly(propSim, t, dt);
            }
            approx_xInteg(xInteg0, accInteg0, dt, 1.0, b, dim,
                        propSim->integBodies, xInteg, xIntegCompCoeffs);
            t += dt;
            check_events(propSim, t, xInteg);
            get_state_der(propSim, t, xInteg, accInteg0);
            propSim->interpParams.tStack.push_back(t);
            propSim->interpParams.xIntegStack.push_back(xInteg);
            propSim->t = t;
            propSim->xInteg = xInteg;
            propSim->integParams.timestepCounter++;
            refine_b(b, e, dtReq/dt, dim);
            check_ca_or_impact(propSim, t-dt, xInteg0, t, xInteg);
            if ((propSim->integParams.tf > propSim->integParams.t0 &&
                    t >= propSim->integParams.tf) ||
                (propSim->integParams.tf < propSim->integParams.t0 &&
                    t <= propSim->integParams.tf)) {
                keepStepping = 0;
            }
            dt = dtReq;
            event_timestep_check(propSim, dt);
            oneStepDone = 1;
        }
    }
    propSim->interpParams.bStack.push_back(b);
    propSim->interpParams.accIntegStack.push_back(accInteg0);
}
