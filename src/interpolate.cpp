#include "interpolate.h"

void approx_xInteg_math(const std::vector<real> &xInteg0,
                        const std::vector<real> &accInteg0, const real &dt,
                        const real &h, const std::vector<std::vector<real>> &b,
                        const size_t starti, const size_t startb,
                        const size_t &iterStep, std::vector<real> &xIntegNext,
                        std::vector<real> &xIntegCompCoeffs) {
    if (h == 1.0) {
        for (size_t j = 0; j < iterStep; j++) {
            xIntegNext[starti+j]          = xInteg0[starti+j];
            xIntegNext[starti+j+iterStep] = xInteg0[starti+j+iterStep];
            comp_sum(b[6][startb+j]*dt*dt/72.0, &(xIntegNext[starti+j]), &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[5][startb+j]*dt*dt/56.0, &(xIntegNext[starti+j]), &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[4][startb+j]*dt*dt/42.0, &(xIntegNext[starti+j]), &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[3][startb+j]*dt*dt/30.0, &(xIntegNext[starti+j]), &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[2][startb+j]*dt*dt/20.0, &(xIntegNext[starti+j]), &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[1][startb+j]*dt*dt/12.0, &(xIntegNext[starti+j]), &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[0][startb+j]*dt*dt/6.0 , &(xIntegNext[starti+j]), &(xIntegCompCoeffs[starti+j]));
            comp_sum(accInteg0[startb+j]*dt*dt/2.0, &(xIntegNext[starti+j]), &(xIntegCompCoeffs[starti+j]));
            comp_sum(xInteg0[starti+j+iterStep]*dt, &(xIntegNext[starti+j]), &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[6][startb+j]*dt/8.0, &(xIntegNext[starti+j+iterStep]), &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(b[5][startb+j]*dt/7.0, &(xIntegNext[starti+j+iterStep]), &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(b[4][startb+j]*dt/6.0, &(xIntegNext[starti+j+iterStep]), &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(b[3][startb+j]*dt/5.0, &(xIntegNext[starti+j+iterStep]), &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(b[2][startb+j]*dt/4.0, &(xIntegNext[starti+j+iterStep]), &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(b[1][startb+j]*dt/3.0, &(xIntegNext[starti+j+iterStep]), &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(b[0][startb+j]*dt/2.0, &(xIntegNext[starti+j+iterStep]), &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(accInteg0[startb+j]*dt, &(xIntegNext[starti+j+iterStep]), &(xIntegCompCoeffs[starti+j+iterStep]));
        }
    } else {
        for (size_t j = 0; j < iterStep; j++) {
            const real inc0 = ((((((((b[6][startb+j]*7.*h/9. + b[5][startb+j])*3.*h/4. + b[4][startb+j])*5.*h/7. + b[3][startb+j])*2.*h/3. + b[2][startb+j])*3.*h/5. + b[1][startb+j])*h/2.    + b[0][startb+j])*h/3. + accInteg0[startb+j])*dt*h/2. + xInteg0[starti+j+iterStep])*dt*h;;
            const real inc1 = (((((((b[6][startb+j]*7.*h/8.  + b[5][startb+j])*6.*h/7. + b[4][startb+j])*5.*h/6. + b[3][startb+j])*4.*h/5. + b[2][startb+j])*3.*h/4. + b[1][startb+j])*2.*h/3. + b[0][startb+j])*h/2. + accInteg0[startb+j])*dt*h;
            xIntegNext[starti+j]          = xInteg0[starti+j] + inc0 - xIntegCompCoeffs[starti+j];
            xIntegNext[starti+j+iterStep] = xInteg0[starti+j+iterStep] + inc1 - xIntegCompCoeffs[starti+j+iterStep];
        }
    }
}

void approx_xInteg(const std::vector<real> &xInteg0,
                   const std::vector<real> &accInteg0, const real &dt,
                   const real &h, const std::vector<std::vector<real>> &b,
                   const std::vector<IntegBody> &integBodies,
                   std::vector<real> &xIntegNext,
                   std::vector<real> &xIntegCompCoeffs) {
    size_t starti = 0;
    size_t startb = 0;
    for (size_t i = 0; i < integBodies.size(); i++) {
        // do pos/vel first, then STM
        approx_xInteg_math(xInteg0, accInteg0, dt, h, b, starti, startb, 3, xIntegNext, xIntegCompCoeffs);
        starti += 6;
        startb += 3;
        if (integBodies[i].propStm) {
            // do 6x6 STM first, then 6x1 fitted parameters
            approx_xInteg_math(xInteg0, accInteg0, dt, h, b, starti, startb, 18, xIntegNext, xIntegCompCoeffs);
            starti += 36;
            startb += 18;
            if (integBodies[i].stm.size() > 36) {
                // do fitted parameters
                const size_t numParams = (integBodies[i].stm.size()-36)/6;
                for (size_t param = 0; param < numParams; param++) {
                    approx_xInteg_math(xInteg0, accInteg0, dt, h, b, starti, startb, 3, xIntegNext, xIntegCompCoeffs);
                    starti += 6;
                    startb += 3;
                }
            }
        }
    }
}

std::vector<real> propSimulation::interpolate(const real t) {
    std::vector<real> xIntegInterp = std::vector<real>(this->xInteg.size(), 0.0);
    // find the index of the last element in this->interpParams.tStack that is
    // less than t if propagating forward, or the first element that is greater
    // than t if propagating backward
    size_t idx = 0;
    const bool forwardProp = this->integParams.t0 < this->integParams.tf;
    const bool backwardProp = this->integParams.t0 > this->integParams.tf;
    if (forwardProp) {
        if (t+this->tEvalMargin < this->integParams.t0 || t-this->tEvalMargin > this->integParams.tf) {
            throw std::runtime_error("The interpolation time is outside the integration time window");
        }
        while (idx < this->interpParams.tStack.size()-1 &&
               this->interpParams.tStack[idx+1] < t) {
            idx++;
        }
    } else if (backwardProp) {
        if (t-this->tEvalMargin > this->integParams.t0 || t+this->tEvalMargin < this->integParams.tf) {
            throw std::runtime_error("The interpolation time is outside the integration time window");
        }
        while (idx < this->interpParams.tStack.size()-1 &&
               this->interpParams.tStack[idx+1] > t) {
            idx++;
        }
    }
    const real t0 = this->interpParams.tStack[idx];
    real dt;
    if (idx == this->interpParams.tStack.size() - 1) {
        real tRef;
        if (forwardProp) {
            tRef = this->integParams.tf + this->tEvalMargin;
        } else {
            tRef = this->integParams.tf - this->tEvalMargin;
        }
        dt = tRef - t0;
    } else {
        dt = this->interpParams.tStack[idx + 1] - t0;
    }
    const real h = (t - t0) / dt;
    std::vector<real> dummyCompCoeffs = std::vector<real>(this->xInteg.size(), 0.0);
    approx_xInteg(this->interpParams.xIntegStack[idx],
                  this->interpParams.accIntegStack[idx], dt, h,
                  this->interpParams.bStack[idx], this->integBodies,
                  xIntegInterp, dummyCompCoeffs);
    return xIntegInterp;
}

void interpolate_on_the_fly(propSimulation *propSim, const real &t, const real &dt) {
    const real tNext = t + dt;
    size_t &interpIdx = propSim->interpIdx;
    const bool forwardProp = propSim->integParams.t0 < propSim->integParams.tf;
    const bool backwardProp = propSim->integParams.t0 > propSim->integParams.tf;
    bool interpIdxInWindow;
    get_interpIdxInWindow(propSim, t, tNext, forwardProp,
                          backwardProp, interpIdxInWindow);
    while (interpIdx < propSim->tEval.size() && interpIdxInWindow) {
        real tInterpGeom;
        if (propSim->tEvalUTC) {
            SpiceDouble etMinusUtc;
            real secPastJ2000Utc;
            mjd_to_et(propSim->tEval[interpIdx], secPastJ2000Utc);
            deltet_c(secPastJ2000Utc, "UTC", &etMinusUtc);
            tInterpGeom = et_to_mjd(secPastJ2000Utc + etMinusUtc);
        } else {
            tInterpGeom = propSim->tEval[interpIdx];
        }
        std::vector<real> xInterpGeom(propSim->xInteg.size(), 0.0);
        evaluate_one_interpolation(propSim, t, dt, tInterpGeom, xInterpGeom);
        if (propSim->evalApparentState) {
            std::vector<real> lightTime(propSim->integParams.nInteg, 0.0);
            std::vector<real> xInterpApparent(propSim->xInteg.size(), 0.0);
            get_lightTime_and_xRelative(propSim, interpIdx, t, dt, tInterpGeom, xInterpGeom,
                                        lightTime, xInterpApparent);
            propSim->lightTimeEval.push_back(lightTime);
            propSim->xIntegEval.push_back(xInterpApparent);
            if (propSim->evalMeasurements) {
                get_measurement(propSim, interpIdx, t, dt, tInterpGeom, xInterpGeom, xInterpApparent);
            }
        } else {
            propSim->xIntegEval.push_back(xInterpGeom);
        }
        interpIdx++;
        get_interpIdxInWindow(propSim, t, tNext, forwardProp,
                              backwardProp, interpIdxInWindow);
    }
}

void evaluate_one_interpolation(const propSimulation *propSim, const real &t,
                                const real &dt, const real &tInterp,
                                std::vector<real> &xInterp) {
    const real h = (tInterp - t) / dt;
    const size_t idx = propSim->interpParams.bStack.size() - 1;
    std::vector<real> dummyCompCoeffs = std::vector<real>(propSim->xInteg.size(), 0.0);
    approx_xInteg(propSim->interpParams.xIntegStack[idx],
                  propSim->interpParams.accIntegStack[idx], dt, h,
                  propSim->interpParams.bStack[idx], propSim->integBodies,
                  xInterp, dummyCompCoeffs);
}

void get_interpIdxInWindow(const propSimulation *propSim,
                           const real &tWindowStart, const real &tNext,
                           const bool &forwardProp,
                           const bool &backwardProp,
                           bool &interpIdxInWindow) {
    interpIdxInWindow = false;
    const size_t interpIdx = propSim->interpIdx;
    bool fwInWindow, fwInWindowMarginStart, fwInWindowMarginEnd;
    bool bwInWindow, bwInWindowMarginStart, bwInWindowMarginEnd;
    fwInWindow =
        (forwardProp && propSim->tEval[interpIdx] >= tWindowStart &&
         propSim->tEval[interpIdx] <= tNext);
    fwInWindowMarginStart =
        (forwardProp &&
         propSim->tEval[interpIdx] <= propSim->integParams.t0 &&
         propSim->tEval[interpIdx] + propSim->tEvalMargin >=
             propSim->integParams.t0 &&
         propSim->tEval[interpIdx] + propSim->tEvalMargin >= tWindowStart &&
         propSim->tEval[interpIdx] + propSim->tEvalMargin <= tNext);
    fwInWindowMarginEnd =
        (forwardProp &&
         propSim->tEval[interpIdx] >= propSim->integParams.tf &&
         propSim->tEval[interpIdx] - propSim->tEvalMargin <=
             propSim->integParams.tf &&
         propSim->tEval[interpIdx] - propSim->tEvalMargin >= tWindowStart &&
         propSim->tEval[interpIdx] - propSim->tEvalMargin <= tNext);
    bwInWindow =
        (backwardProp && propSim->tEval[interpIdx] <= tWindowStart &&
         propSim->tEval[interpIdx] >= tNext);
    bwInWindowMarginStart =
        (backwardProp &&
         propSim->tEval[interpIdx] >= propSim->integParams.t0 &&
         propSim->tEval[interpIdx] - propSim->tEvalMargin <=
             propSim->integParams.t0 &&
         propSim->tEval[interpIdx] - propSim->tEvalMargin <= tWindowStart &&
         propSim->tEval[interpIdx] - propSim->tEvalMargin >= tNext);
    bwInWindowMarginEnd =
        (backwardProp &&
         propSim->tEval[interpIdx] <= propSim->integParams.tf &&
         propSim->tEval[interpIdx] + propSim->tEvalMargin >=
             propSim->integParams.tf &&
         propSim->tEval[interpIdx] + propSim->tEvalMargin <= tWindowStart &&
         propSim->tEval[interpIdx] + propSim->tEvalMargin >= tNext);
    interpIdxInWindow = fwInWindow || fwInWindowMarginStart ||
        fwInWindowMarginEnd || bwInWindow || bwInWindowMarginStart ||
        bwInWindowMarginEnd;
}

void get_lightTime_and_xRelative(propSimulation *propSim,
                                 const size_t interpIdx, const real &t,
                                 const real &dt, const real tInterpGeom,
                                 const std::vector<real> &xInterpGeom,
                                 std::vector<real> &lightTime,
                                 std::vector<real> &xInterpApparent) {
    size_t numStates = xInterpGeom.size();
    std::vector<real> xObserver = propSim->xObserver[interpIdx];
    bool bouncePointAtLeadingEdge = false;
    if (propSim->observerInfo[interpIdx].size() == 9 ||
        propSim->observerInfo[interpIdx].size() == 10) {
        bouncePointAtLeadingEdge = propSim->observerInfo[interpIdx][8] == 1.0;
    }
    size_t starti = 0;
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        real lightTimeTemp;
        std::vector<real> xInterpApparentTemp(numStates, 0.0);
        std::vector<real> xInterpApparentBary(6, 0.0);
        get_lightTimeOneBody(propSim, i, tInterpGeom, xInterpGeom, xObserver,
                             bouncePointAtLeadingEdge, t, dt,
                             lightTimeTemp);
        evaluate_one_interpolation(propSim, t, dt, tInterpGeom-lightTimeTemp, xInterpApparentTemp);
        lightTime[i] = lightTimeTemp;
        for (size_t j = 0; j < 6; j++) {
            xInterpApparentBary[j] = xInterpApparentTemp[starti + j];
        }
        get_glb_correction(propSim, tInterpGeom, xInterpApparentBary);
        for (size_t j = 0; j < 6; j++) {
            xInterpApparent[starti + j] = xInterpApparentBary[j] - xObserver[j];
        }
        starti += 2*propSim->integBodies[i].n2Derivs;
    }
}

void get_lightTimeOneBody(propSimulation *propSim, const size_t &i,
                          const real tInterpGeom, std::vector<real> xInterpGeom,
                          std::vector<real> xObserver,
                          const bool bouncePointAtLeadingEdge, const real &t,
                          const real &dt, real &lightTimeOneBody) {
    size_t numStates = xInterpGeom.size();
    std::vector<real> xInterpApparentFull(numStates, 0.0);
    std::vector<real> xInterpApparent(6, 0.0);
    std::vector<real> xRelativeOneBody(6, 0.0);
    real distRelativeOneBody;

    size_t starti = 0;
    for (size_t j = 0; j < i; j++) {
        starti += 2*propSim->integBodies[j].n2Derivs;
    }
    for (size_t j = 0; j < 6; j++) {
        xRelativeOneBody[j] = xInterpGeom[starti + j] - xObserver[j];
    }
    vnorm({xRelativeOneBody[0], xRelativeOneBody[1], xRelativeOneBody[2]},
          distRelativeOneBody);
    if (bouncePointAtLeadingEdge) {
        distRelativeOneBody -= propSim->integBodies[i].radius;
    }
    lightTimeOneBody = distRelativeOneBody / propSim->consts.clight;
    if (propSim->convergedLightTime) {
        real lightTimeTol = 1e-10 / 86400.0L;
        real lightTimeOneBodyPrev = 0.0L;
        real deltaLightTimeRelativistic;
        size_t maxIter = 20;
        size_t iter = 0;
        // keep iterating until max iterations or light time tolerance is met
        while (iter < maxIter &&
               fabs(lightTimeOneBody - lightTimeOneBodyPrev) > lightTimeTol) {
            evaluate_one_interpolation(propSim, t, dt, tInterpGeom - lightTimeOneBody,
                                       xInterpApparentFull);
            for (size_t j = 0; j < 6; j++) {
                xRelativeOneBody[j] =
                    xInterpApparentFull[starti + j] - xObserver[j];
                xInterpApparent[j] = xInterpApparentFull[starti + j];
            }
            vnorm(
                {xRelativeOneBody[0], xRelativeOneBody[1], xRelativeOneBody[2]},
                distRelativeOneBody);
            if (bouncePointAtLeadingEdge) {
                distRelativeOneBody -= propSim->integBodies[i].radius;
            }
            lightTimeOneBodyPrev = lightTimeOneBody;
            get_delta_delay_relativistic(
                propSim, tInterpGeom - lightTimeOneBody, xInterpApparent,
                deltaLightTimeRelativistic);
            lightTimeOneBody = distRelativeOneBody / propSim->consts.clight +
                deltaLightTimeRelativistic;
            iter++;
        }
        if (iter >= maxIter) {
            std::cout
                << "Warning: Downleg light time did not converge for body "
                << propSim->integBodies[i].name << " at time " << tInterpGeom
                << ", change from previous iteration was "
                << fabs(lightTimeOneBody - lightTimeOneBodyPrev) << std::endl;
        }
    }
}

void get_glb_correction(propSimulation *propSim, const real &tInterpGeom,
                        std::vector<real> &xInterpApparentBary) {
    double sunState[9];
    double earthState[9];
    get_spk_state(10, tInterpGeom, propSim->ephem, sunState);
    get_spk_state(399, tInterpGeom, propSim->ephem, earthState);

    std::vector<real> sunEarthPos = {earthState[0] - sunState[0],
                                     earthState[1] - sunState[1],
                                     earthState[2] - sunState[2]};
    real sunEarthDist;
    vnorm(sunEarthPos, sunEarthDist);
    std::vector<real> sunTargetPos = {xInterpApparentBary[0] - sunState[0],
                                      xInterpApparentBary[1] - sunState[1],
                                      xInterpApparentBary[2] - sunState[2]};
    real sunTargetDist;
    vnorm(sunTargetPos, sunTargetDist);
    std::vector<real> earthTargetPos = {xInterpApparentBary[0] - earthState[0],
                                        xInterpApparentBary[1] - earthState[1],
                                        xInterpApparentBary[2] - earthState[2]};
    real earthTargetDist;
    vnorm(earthTargetPos, earthTargetDist);

    real G = propSim->consts.G;
    real sunGM = 0;
    for (size_t i = 0; i < propSim->integParams.nSpice; i++) {
        if (propSim->spiceBodies[i].spiceId == 10) {
            sunGM = G * propSim->spiceBodies[i].mass;
        }
    }
    if (sunGM == 0) {
        throw std::runtime_error(
            "Sun GM not found in get_delta_delay_relativistic");
    }
    real c = propSim->consts.clight;

    // from section 7.2.4 in the Explanatory Supplement to the Astronomical
    // Almanac, 3rd edition
    std::vector<real> e(3, 0.0);
    vunit(sunEarthPos, e);
    std::vector<real> q(3, 0.0);
    vunit(sunTargetPos, q);
    std::vector<real> p(3, 0.0);
    vunit(earthTargetPos, p);
    std::vector<real> deltaP1Targ(3, 0.0);
    std::vector<real> deltaP1Star(3, 0.0);
    std::vector<real> p1(3, 0.0);

    real pDotQ, eDotP, qDotE;
    vdot(p, q, pDotQ);
    vdot(e, p, eDotP);
    vdot(q, e, qDotE);

    real g1 = 2 * sunGM / c / c / sunEarthDist;
    real g2 = 1.0L + qDotE;

    deltaP1Targ[0] = g1 * (pDotQ * e[0] - eDotP * q[0]) / g2;
    deltaP1Targ[1] = g1 * (pDotQ * e[1] - eDotP * q[1]) / g2;
    deltaP1Targ[2] = g1 * (pDotQ * e[2] - eDotP * q[2]) / g2;

    deltaP1Star[0] = g1 * (e[0] - eDotP * p[0]) / (1 + eDotP);
    deltaP1Star[1] = g1 * (e[1] - eDotP * p[1]) / (1 + eDotP);
    deltaP1Star[2] = g1 * (e[2] - eDotP * p[2]) / (1 + eDotP);

    p1[0] = p[0] - deltaP1Star[0] + deltaP1Targ[0];
    p1[1] = p[1] - deltaP1Star[1] + deltaP1Targ[1];
    p1[2] = p[2] - deltaP1Star[2] + deltaP1Targ[2];

    earthTargetPos[0] = earthTargetDist * p1[0];
    earthTargetPos[1] = earthTargetDist * p1[1];
    earthTargetPos[2] = earthTargetDist * p1[2];

    xInterpApparentBary[0] = earthState[0] + earthTargetPos[0];
    xInterpApparentBary[1] = earthState[1] + earthTargetPos[1];
    xInterpApparentBary[2] = earthState[2] + earthTargetPos[2];
}

void get_measurement(propSimulation *propSim, const size_t &interpIdx,
                     const real &t, const real &dt, const real tInterpGeom,
                     const std::vector<real> &xInterpGeom,
                     const std::vector<real> &xInterpApparent) {
    std::vector<real> opticalMeasurement(2*propSim->integParams.nInteg,
                                       std::numeric_limits<real>::quiet_NaN());
    std::vector<real> opticalPartials(12*propSim->integParams.nInteg,
                                       std::numeric_limits<real>::quiet_NaN());
    std::vector<real> radarMeasurement(propSim->integParams.nInteg,
                                       std::numeric_limits<real>::quiet_NaN());
    std::vector<real> radarPartials(6*propSim->integParams.nInteg,
                                       std::numeric_limits<real>::quiet_NaN());
    switch (propSim->radarObserver[interpIdx]) {
    case 0:
        get_optical_measurement(propSim, xInterpApparent, opticalMeasurement,
                                opticalPartials);
        break;
    case 1: case 2:
        get_radar_measurement(propSim, interpIdx, t, dt, tInterpGeom,
                              xInterpGeom, radarMeasurement, radarPartials);
        break;
    default:
        throw std::runtime_error(
            "get_measurement: radarObserver flag must be 0, 1, or 2");
        break;
    }
    propSim->opticalObs.push_back(opticalMeasurement);
    propSim->opticalPartials.push_back(opticalPartials);
    propSim->radarObs.push_back(radarMeasurement);
    propSim->radarPartials.push_back(radarPartials);
}

void get_optical_measurement(propSimulation *propSim,
                             const std::vector<real> &xInterpApparent,
                             std::vector<real> &opticalMeasurement,
                             std::vector<real> &opticalPartials) {
    size_t starti = 0;
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        std::vector<real> xInterpApparentOneBody(6, 0.0);
        for (size_t j = 0; j < 6; j++) {
            xInterpApparentOneBody[j] = xInterpApparent[starti + j];
        }
        real dist;
        vnorm({xInterpApparentOneBody[0], xInterpApparentOneBody[1],
                xInterpApparentOneBody[2]}, dist);
        real r_asc = atan2(xInterpApparentOneBody[1],
                           xInterpApparentOneBody[0]);
        if (r_asc < 0) {
            r_asc = r_asc + 2*PI;
        }
        const real x2py2 = xInterpApparentOneBody[0]*xInterpApparentOneBody[0] +
            xInterpApparentOneBody[1]*xInterpApparentOneBody[1];
        const real dradx = -xInterpApparentOneBody[1]/x2py2;
        const real drady = xInterpApparentOneBody[0]/x2py2;
        real dec = asin(xInterpApparentOneBody[2]/dist);
        const real r = sqrt(x2py2 + xInterpApparentOneBody[2]*xInterpApparentOneBody[2]);
        const real ddecdx = -xInterpApparentOneBody[0]*xInterpApparentOneBody[2]/r/r/sqrt(x2py2);
        const real ddecdy = -xInterpApparentOneBody[1]*xInterpApparentOneBody[2]/r/r/sqrt(x2py2);
        const real ddecdz = sqrt(x2py2)/r/r;
        const real conv = 180.0L/PI*3600.0L; // radians -> arcsec
        r_asc *= conv;
        dec *= conv;
        opticalMeasurement[2*i] = r_asc;
        opticalMeasurement[2*i+1] = dec;
        std::fill(opticalPartials.begin()+12*i, opticalPartials.begin()+12*(i+1), 0.0);
        opticalPartials[12*i] = dradx*conv;
        opticalPartials[12*i+1] = drady*conv;
        opticalPartials[12*i+6] = ddecdx*conv;
        opticalPartials[12*i+7] = ddecdy*conv;
        opticalPartials[12*i+8] = ddecdz*conv;
        starti += 2*propSim->integBodies[i].n2Derivs;
    }
}

void get_radar_measurement(propSimulation *propSim, const size_t &interpIdx,
                           const real &t, const real &dt,
                           const real tInterpGeom,
                           const std::vector<real> &xInterpGeom,
                           std::vector<real> &radarMeasurement,
                           std::vector<real> &radarPartials) {
    if (propSim->observerInfo[interpIdx].size() != 9 &&
        propSim->observerInfo[interpIdx].size() != 10) {
        throw std::runtime_error(
            "Error: observerInfo must be a 9 or 10 "
            "element vector for radar measurements");
    }
    real receiveTimeTDB = tInterpGeom;
    real transmitTimeTDB;
    std::vector<real> xTrgtBaryRcv = xInterpGeom;
    std::vector<real> xObsBaryRcv(6, 0.0);
    std::vector<real> xTrgtBaryBounce(6, 0.0);
    std::vector<real> xObsBaryTx(6, 0.0);
    real transmitFreq = 0.0;
    if (propSim->radarObserver[interpIdx] == 2) {
        transmitFreq = propSim->observerInfo[interpIdx][9];
    }
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        std::fill(radarPartials.begin()+6*i, radarPartials.begin()+6*(i+1), 0.0);
        real delayMeasurement;
        get_delay_measurement(propSim, interpIdx, t, dt, i, tInterpGeom,
                              xInterpGeom, receiveTimeTDB, transmitTimeTDB,
                              xObsBaryRcv, xTrgtBaryBounce, xObsBaryTx,
                              delayMeasurement, radarPartials);
        if (propSim->radarObserver[interpIdx] == 1) {
            radarMeasurement[i] = delayMeasurement;
        } else if (propSim->radarObserver[interpIdx] == 2) {
            real dopplerMeasurement;
            get_doppler_measurement(propSim, i, receiveTimeTDB, transmitTimeTDB,
                                    xObsBaryRcv, xTrgtBaryBounce, xObsBaryTx,
                                    transmitFreq, dopplerMeasurement, radarPartials);
            radarMeasurement[i] = dopplerMeasurement;
        }
    } 
}

void get_delay_measurement(propSimulation *propSim, const size_t &interpIdx,
                           const real &t, const real &dt, const size_t &i,
                           const real tInterpGeom,
                           const std::vector<real> &xInterpGeom,
                           const real &receiveTimeTDB, real &transmitTimeTDB,
                           std::vector<real> &xObsBaryRcv,
                           std::vector<real> &xTrgtBaryBounce,
                           std::vector<real> &xObsBaryTx, real &delayMeasurement,
                           std::vector<real> &delayPartials) {
    size_t numStates = xInterpGeom.size();
    std::vector<real> receiverInfo = {propSim->observerInfo[interpIdx][0],
                                      propSim->observerInfo[interpIdx][1],
                                      propSim->observerInfo[interpIdx][2],
                                      propSim->observerInfo[interpIdx][3]};
    get_observer_state(receiveTimeTDB, receiverInfo, propSim, false,
                       xObsBaryRcv);
    std::vector<real> transmitterInfo = {propSim->observerInfo[interpIdx][4],
                                         propSim->observerInfo[interpIdx][5],
                                         propSim->observerInfo[interpIdx][6],
                                         propSim->observerInfo[interpIdx][7]};
    bool bouncePointAtLeadingEdge = propSim->observerInfo[interpIdx][8] == 1.0;
    real delayDownleg;
    real bounceTimeTDB;
    real delayUpleg;
    std::vector<real> xTrgtBaryBounceAllBody(numStates, 0.0);
    // downleg delay is the already evaluated light time
    delayDownleg = propSim->lightTimeEval[interpIdx][i];
    bounceTimeTDB = receiveTimeTDB - delayDownleg;
    evaluate_one_interpolation(propSim, t, dt, bounceTimeTDB,
                                    xTrgtBaryBounceAllBody);
    size_t starti = 0;
    for (size_t j = 0; j < i; j++) {
        starti += 2*propSim->integBodies[j].n2Derivs;
    }
    for (size_t j = 0; j < 6; j++) {
        xTrgtBaryBounce[j] = xTrgtBaryBounceAllBody[starti + j];
    }
    // iterate to get upleg delay
    delayUpleg = delayDownleg;
    if (propSim->convergedLightTime) {
        real distRelativeUpleg;
        real lightTimeTol = 1e-10 / 86400.0L;
        real delayUplegPrev = 0.0L;
        real deltaDelayUplegRelativistic;
        size_t maxIter = 20;
        size_t iter = 0;
        while (iter < maxIter &&
                fabs(delayUpleg - delayUplegPrev) > lightTimeTol) {
            get_observer_state(bounceTimeTDB - delayUpleg, transmitterInfo,
                                propSim, false, xObsBaryTx);
            std::vector<real> xRelativeOneBody(6, 0.0);
            for (size_t j = 0; j < 6; j++) {
                xRelativeOneBody[j] = xTrgtBaryBounce[j] - xObsBaryTx[j];
            }
            vnorm({xRelativeOneBody[0], xRelativeOneBody[1],
                    xRelativeOneBody[2]},
                    distRelativeUpleg);
            if (bouncePointAtLeadingEdge) {
                distRelativeUpleg -= propSim->integBodies[i].radius;
            }
            delayUplegPrev = delayUpleg;
            get_delta_delay_relativistic(
                propSim, bounceTimeTDB - delayUpleg, xTrgtBaryBounce,
                deltaDelayUplegRelativistic);
            delayUpleg = distRelativeUpleg / propSim->consts.clight +
                deltaDelayUplegRelativistic;
            iter++;
        }
        if (iter >= maxIter) {
            std::cout
                << "Warning: Upleg light time did not converge for body "
                << propSim->integBodies[i].name << " at time "
                << tInterpGeom << ", change from previous iteration was "
                << fabs(delayUpleg - delayUplegPrev) << std::endl;
        }
    }
    transmitTimeTDB = bounceTimeTDB - delayUpleg;
    get_observer_state(transmitTimeTDB, transmitterInfo, propSim,
                        false, xObsBaryTx);
    // get delay measurement
    delayMeasurement = (delayDownleg + delayUpleg) * 86400.0L *
        1e6;  // days -> seconds -> microseconds
    if (propSim->tEvalUTC) {
        SpiceDouble etMinusUtcReceiveTime;
        SpiceDouble etMinusUtcTransmitTime;
        real receiveTimeET;
        real transmitTimeET;
        mjd_to_et(receiveTimeTDB, receiveTimeET);
        mjd_to_et(transmitTimeTDB, transmitTimeET);
        deltet_c(receiveTimeET, "ET", &etMinusUtcReceiveTime);
        deltet_c(transmitTimeET, "ET", &etMinusUtcTransmitTime);
        delayMeasurement +=
            (etMinusUtcTransmitTime - etMinusUtcReceiveTime) * 1e6;
    }

    std::vector<real> uplegState(6, 0.0);
    std::vector<real> downlegState(6, 0.0);
    for (size_t j = 0; j < 6; j++) {
        uplegState[j] = xTrgtBaryBounce[j] - xObsBaryTx[j];
        downlegState[j] = xTrgtBaryBounce[j] - xObsBaryRcv[j];
    }
    const real uplegDist = sqrt(uplegState[0]*uplegState[0] +
                                uplegState[1]*uplegState[1] +
                                uplegState[2]*uplegState[2]);
    const real downlegDist = sqrt(downlegState[0]*downlegState[0] +
                                  downlegState[1]*downlegState[1] +
                                  downlegState[2]*downlegState[2]);
    for (size_t j = 0; j < 3; j++){
        delayPartials[6*i+j] = 1/propSim->consts.clight * (uplegState[j]/uplegDist + downlegState[j]/downlegDist);
        delayPartials[6*i+j] *= 86400.0L * 1.0e6;  // days -> seconds -> microseconds
        delayPartials[6*i+3+j] = 0.0L;
    }
}

void get_delta_delay_relativistic(propSimulation *propSim,
                                  const real &tForSpice,
                                  const std::vector<real> &targetState,
                                  real &deltaDelayRelativistic) {
    // from Standish (1990),
    // https://ui.adsabs.harvard.edu/abs/1990A&A...233..252S
    double sunState[9];
    double earthState[9];
    get_spk_state(10, tForSpice, propSim->ephem, sunState);
    get_spk_state(399, tForSpice, propSim->ephem, earthState);

    std::vector<real> sunEarthPos = {earthState[0] - sunState[0],
                                     earthState[1] - sunState[1],
                                     earthState[2] - sunState[2]};
    real sunEarthDist;
    vnorm(sunEarthPos, sunEarthDist);
    std::vector<real> sunTargetPos = {targetState[0] - sunState[0],
                                      targetState[1] - sunState[1],
                                      targetState[2] - sunState[2]};
    real sunTargetDist;
    vnorm(sunTargetPos, sunTargetDist);
    std::vector<real> earthTargetPos = {targetState[0] - earthState[0],
                                        targetState[1] - earthState[1],
                                        targetState[2] - earthState[2]};
    real earthTargetDist;
    vnorm(earthTargetPos, earthTargetDist);

    real G = propSim->consts.G;
    real sunGM = 0;
    for (size_t i = 0; i < propSim->integParams.nSpice; i++) {
        if (propSim->spiceBodies[i].spiceId == 10) {
            sunGM = G * propSim->spiceBodies[i].mass;
        }
    }
    if (sunGM == 0) {
        throw std::runtime_error(
            "Sun GM not found in get_delta_delay_relativistic");
    }
    real c = propSim->consts.clight;
    real gamma = 1.0L;  // PPN parameter

    deltaDelayRelativistic = (1 + gamma) * sunGM * pow(c, -3) *
        log((sunEarthDist + sunTargetDist + earthTargetDist) /
            (sunEarthDist + sunTargetDist - earthTargetDist));
}

void get_doppler_measurement(propSimulation *propSim, const size_t &i,
                             const real receiveTimeTDB,
                             const real transmitTimeTDB,
                             const std::vector<real> xObsBaryRcv,
                             const std::vector<real> xTrgtBaryBounce,
                             const std::vector<real> xObsBaryTx,
                             const real transmitFreq, real &dopplerMeasurement,
                             std::vector<real> &dopplerPartials) {
    // based on "Mathematical Formulation of the Double-Precision Orbit
    // Determination Program (DPODP)" by T.D.Moyer (1971),
    // https://ntrs.nasa.gov/citations/19710017134
    std::vector<real> pos3(3), vel3(3);
    std::vector<real> pos2(3), vel2(3);
    std::vector<real> pos1(3), vel1(3);
    pos3 = {xObsBaryRcv[0], xObsBaryRcv[1], xObsBaryRcv[2]};
    vel3 = {xObsBaryRcv[3], xObsBaryRcv[4], xObsBaryRcv[5]};
    pos2 = {xTrgtBaryBounce[0], xTrgtBaryBounce[1], xTrgtBaryBounce[2]};
    vel2 = {xTrgtBaryBounce[3], xTrgtBaryBounce[4], xTrgtBaryBounce[5]};
    pos1 = {xObsBaryTx[0], xObsBaryTx[1], xObsBaryTx[2]};
    vel1 = {xObsBaryTx[3], xObsBaryTx[4], xObsBaryTx[5]};

    real r3, r2, r1;
    vnorm(pos3, r3);
    vnorm(pos2, r2);
    vnorm(pos1, r1);
    real v3, v1;
    vnorm(vel3, v3);
    vnorm(vel1, v1);

    real rdot1, rdot2, rdot3;
    vdot(pos1, vel1, rdot1);
    rdot1 /= r1;
    vdot(pos2, vel2, rdot2);
    rdot2 /= r2;
    vdot(pos3, vel3, rdot3);
    rdot3 /= r3;

    double xSun3[9];
    double xSun1[9];
    get_spk_state(10, receiveTimeTDB, propSim->ephem, xSun3);
    get_spk_state(10, transmitTimeTDB, propSim->ephem, xSun1);

    std::vector<real> posHelio3(3), velHelio3(3), posHelio1(3), velHelio1(3);
    posHelio3[0] = pos3[0] - xSun3[0];
    posHelio3[1] = pos3[1] - xSun3[1];
    posHelio3[2] = pos3[2] - xSun3[2];
    velHelio3[0] = vel3[0] - xSun3[3];
    velHelio3[1] = vel3[1] - xSun3[4];
    velHelio3[2] = vel3[2] - xSun3[5];
    posHelio1[0] = pos1[0] - xSun1[0];
    posHelio1[1] = pos1[1] - xSun1[1];
    posHelio1[2] = pos1[2] - xSun1[2];
    velHelio1[0] = vel1[0] - xSun1[3];
    velHelio1[1] = vel1[1] - xSun1[4];
    velHelio1[2] = vel1[2] - xSun1[5];
    real rHelio3, vHelio3, rHelio1, vHelio1;
    vnorm(posHelio3, rHelio3);
    vnorm(velHelio3, vHelio3);
    vnorm(posHelio1, rHelio1);
    vnorm(velHelio1, vHelio1);

    std::vector<real> pos23(3), vel23(3);
    std::vector<real> pos12(3), vel12(3);
    vsub(pos3, pos2, pos23);
    vsub(vel3, vel2, vel23);
    vsub(pos2, pos1, pos12);
    vsub(vel2, vel1, vel12);

    real r23, r12;
    vnorm(pos23, r23);
    vnorm(pos12, r12);

    real rdot23, rdot12;
    vdot(pos23, vel23, rdot23);
    rdot23 /= r23;
    vdot(pos12, vel12, rdot12);
    rdot12 /= r12;

    real pdot23, pdot12;
    vdot(pos23, vel2, pdot23);
    pdot23 /= r23;
    vdot(pos12, vel1, pdot12);
    pdot12 /= r12;

    real G = propSim->consts.G;
    real sunGM = 0;
    for (size_t i = 0; i < propSim->integParams.nSpice; i++) {
        if (propSim->spiceBodies[i].spiceId == 10) {
            sunGM = G * propSim->spiceBodies[i].mass;
        }
    }
    if (sunGM == 0) {
        throw std::runtime_error("Sun GM not found in get_doppler_measurement");
    }
    real c = propSim->consts.clight;

    real fac1, term1;
    fac1 = -transmitFreq / c;
    term1 = fac1 * (rdot12 + rdot23);

    real fac2, term2, term2a, term2b, term2c, term2d, term2e;
    fac2 = -transmitFreq / c / c;
    term2a = rdot12 * pdot12;
    term2b = rdot23 * pdot23;
    term2c = -rdot12 * rdot23;
    term2d = sunGM * (1 / rHelio1 - 1 / rHelio3);
    term2e = 0.5L * (v1 * v1 - v3 * v3);
    term2 = fac2 * (term2a + term2b + term2c + term2d + term2e);

    real fac3, term3, term3a, term3b, term3c, term3d, term3e, gamma, eps12,
        eps23;
    fac3 = -transmitFreq / c / c / c;
    term3a = rdot12 * pdot12 * pdot12;
    term3b = rdot23 * pdot23 * pdot23;
    term3c = -rdot12 * rdot23 * (pdot12 + pdot23);
    term3d = (rdot12 + rdot23) *
        (sunGM * (1 / rHelio1 - 1 / rHelio3) + 0.5L * (v1 * v1 - v3 * v3));
    gamma = 1.0L;
    eps12 = (rdot1 + rdot2 - rdot12) / (r1 + r2 - r12) -
        (rdot1 + rdot2 + rdot12) / (r1 + r2 + r12);
    eps23 = (rdot2 + rdot3 - rdot23) / (r2 + r3 - r23) -
        (rdot2 + rdot3 + rdot23) / (r2 + r3 + r23);
    term3e = (1 + gamma) * sunGM * (eps12 + eps23);
    term3 = fac3 * (term3a + term3b + term3c + term3d + term3e);

    real dopplerShift;
    dopplerShift = term1 + term2 + term3;
    dopplerMeasurement = dopplerShift;

    std::vector<real> uplegState(6, 0.0);
    std::vector<real> downlegState(6, 0.0);
    for (size_t j = 0; j < 6; j++) {
        uplegState[j] = xTrgtBaryBounce[j] - xObsBaryTx[j];
        downlegState[j] = xTrgtBaryBounce[j] - xObsBaryRcv[j];
    }
    const real uplegDist = sqrt(uplegState[0]*uplegState[0] +
                                uplegState[1]*uplegState[1] +
                                uplegState[2]*uplegState[2]);
    const real downlegDist = sqrt(downlegState[0]*downlegState[0] +
                                  downlegState[1]*downlegState[1] +
                                  downlegState[2]*downlegState[2]);
    const real uplegRhoDot = (uplegState[0]*uplegState[3] +
                              uplegState[1]*uplegState[4] +
                              uplegState[2]*uplegState[5]) / uplegDist;
    const real downlegRhoDot = (downlegState[0]*downlegState[3] +
                                downlegState[1]*downlegState[4] +
                                downlegState[2]*downlegState[5]) / downlegDist;
    // dopplerPartials already has the delay partials since delay measurements gets called first
    for (size_t j = 0; j < 3; j++) {
        dopplerPartials[6*i+3+j] = -transmitFreq * dopplerPartials[6*i+j] / 1.0e6/86400.0L; // convert partial microsec -> sec -
        dopplerPartials[6*i+j] = transmitFreq *
            ((uplegRhoDot * uplegState[j] / uplegDist - uplegState[3 + j]) / uplegDist +
             (downlegRhoDot * downlegState[j] / downlegDist - downlegState[3 + j]) / downlegDist) / c;
    }
}
