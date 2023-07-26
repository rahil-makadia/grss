#include "interpolate.h"

void interpolate(const real &t, const real &dt,
                 const std::vector<real> &xInteg0,
                 const std::vector<real> &accInteg0,
                 const std::vector<std::vector<real>> &b,
                 propSimulation *propSim) {
    size_t nh = 8;
    std::vector<real> tVecForInterp(nh, 0.0);
    std::vector<std::vector<real>> xIntegForInterp(
        nh, std::vector<real>(xInteg0.size(), 0.0));
    static std::vector<real> tVecForInterpPrev(tVecForInterp.size(), 0.0);
    static std::vector<std::vector<real>> xIntegForInterpPrev(
        xIntegForInterp.size(),
        std::vector<real>(xIntegForInterp[0].size(), 0.0));

    tVecForInterp[0] = t;
    xIntegForInterp[0] = xInteg0;
    for (size_t hIdx = 1; hIdx < nh; hIdx++) {
        tVecForInterp[hIdx] = t + hVec[hIdx] * dt;
        approx_xInteg(xInteg0, accInteg0, xIntegForInterp[hIdx], dt, hVec[hIdx],
                      b, propSim->integParams.nInteg);
    }
    one_timestep_interpolation(t + dt, tVecForInterp, xIntegForInterp,
                               tVecForInterpPrev, xIntegForInterpPrev, propSim);
    tVecForInterpPrev = tVecForInterp;
    xIntegForInterpPrev = xIntegForInterp;
}

void get_coeffs(const std::vector<real> &tVecForInterp,
                const std::vector<std::vector<real>> &xIntegForInterp,
                std::vector<std::vector<real>> &coeffs) {
    size_t tLen = tVecForInterp.size();
    size_t numStates = xIntegForInterp[0].size();
    std::vector<std::vector<std::vector<real>>> c(
        numStates,
        std::vector<std::vector<real>>(tLen, std::vector<real>(tLen, 0.0)));
    for (size_t i = 0; i < numStates; i++) {
        for (size_t j = 0; j < tLen; j++) {
            c[i][j][0] = xIntegForInterp[j][i];
        }
        for (size_t j = 1; j < tLen; j++) {
            for (size_t k = 0; k < tLen - j; k++) {
                c[i][k][j] = (c[i][k + 1][j - 1] - c[i][k][j - 1]) /
                    (tVecForInterp[k + j] - tVecForInterp[k]);
            }
        }
    }

    for (size_t i = 0; i < numStates; i++) {
        for (size_t j = 0; j < tLen; j++) {
            coeffs[i][j] = c[i][0][j];
        }
    }
}

void evaluate_one_interpolation(
    const propSimulation *propSim, const real &tInterp,
    const std::vector<real> &tVecForInterp,
    const std::vector<std::vector<real>> &coeffs,
    const std::vector<real> &tVecForInterpPrev,
    const std::vector<std::vector<real>> &coeffsPrev,
    std::vector<real> &xInterp) {
    bool forwardIntegrate = tVecForInterp[0] < tVecForInterp[1];
    bool backwardIntegrate = tVecForInterp[0] > tVecForInterp[1];
    size_t tLen = tVecForInterp.size();
    size_t numStates = xInterp.size();
    size_t n = tLen - 1;

    std::vector<std::vector<real>> coeffsToUse;
    std::vector<real> tVecForInterpToUse;
    if ((forwardIntegrate &&
         (tInterp < tVecForInterp[0] && tInterp != propSim->integParams.t0)) ||
        (backwardIntegrate &&
         (tInterp > tVecForInterp[0] && tInterp != propSim->integParams.t0))) {
        coeffsToUse = coeffsPrev;
        tVecForInterpToUse = tVecForInterpPrev;
    } else {
        coeffsToUse = coeffs;
        tVecForInterpToUse = tVecForInterp;
    }

    for (size_t i = 0; i < numStates; i++) {
        xInterp[i] = coeffsToUse[i][n];
    }
    for (size_t i = 1; i < n + 1; i++) {
        for (size_t j = 0; j < numStates; j++) {
            xInterp[j] = coeffsToUse[j][n - i] +
                (tInterp - tVecForInterpToUse[n - i]) * xInterp[j];
        }
    }
}

void get_interpIdxInWindow(const propSimulation *propSim,
                           const real &tWindowStart, const real &tNext,
                           const bool &forwardIntegrate,
                           const bool &backwardIntegrate,
                           bool &interpIdxInWindow) {
    interpIdxInWindow = false;
    const size_t interpIdx = propSim->interpIdx;
    bool fwInWindow, fwInWindowMarginStart, fwInWindowMarginEnd;
    bool bwInWindow, bwInWindowMarginStart, bwInWindowMarginEnd;
    fwInWindow =
        (forwardIntegrate && propSim->tEval[interpIdx] >= tWindowStart &&
         propSim->tEval[interpIdx] <= tNext);
    fwInWindowMarginStart =
        (forwardIntegrate &&
         propSim->tEval[interpIdx] <= propSim->integParams.t0 &&
         propSim->tEval[interpIdx] + propSim->tEvalMargin >=
             propSim->integParams.t0 &&
         propSim->tEval[interpIdx] + propSim->tEvalMargin >= tWindowStart &&
         propSim->tEval[interpIdx] + propSim->tEvalMargin <= tNext);
    fwInWindowMarginEnd =
        (forwardIntegrate &&
         propSim->tEval[interpIdx] >= propSim->integParams.tf &&
         propSim->tEval[interpIdx] - propSim->tEvalMargin <=
             propSim->integParams.tf &&
         propSim->tEval[interpIdx] - propSim->tEvalMargin >= tWindowStart &&
         propSim->tEval[interpIdx] - propSim->tEvalMargin <= tNext);
    bwInWindow =
        (backwardIntegrate && propSim->tEval[interpIdx] <= tWindowStart &&
         propSim->tEval[interpIdx] >= tNext);
    bwInWindowMarginStart =
        (backwardIntegrate &&
         propSim->tEval[interpIdx] >= propSim->integParams.t0 &&
         propSim->tEval[interpIdx] - propSim->tEvalMargin <=
             propSim->integParams.t0 &&
         propSim->tEval[interpIdx] - propSim->tEvalMargin <= tWindowStart &&
         propSim->tEval[interpIdx] - propSim->tEvalMargin >= tNext);
    bwInWindowMarginEnd =
        (backwardIntegrate &&
         propSim->tEval[interpIdx] <= propSim->integParams.tf &&
         propSim->tEval[interpIdx] + propSim->tEvalMargin >=
             propSim->integParams.tf &&
         propSim->tEval[interpIdx] + propSim->tEvalMargin <= tWindowStart &&
         propSim->tEval[interpIdx] + propSim->tEvalMargin >= tNext);
    interpIdxInWindow = fwInWindow || fwInWindowMarginStart ||
        fwInWindowMarginEnd || bwInWindow || bwInWindowMarginStart ||
        bwInWindowMarginEnd;
}

void one_timestep_interpolation(
    const real &tNext, const std::vector<real> &tVecForInterp,
    const std::vector<std::vector<real>> &xIntegForInterp,
    std::vector<real> &tVecForInterpPrev,
    std::vector<std::vector<real>> &xIntegForInterpPrev,
    propSimulation *propSim) {
    size_t tLen = tVecForInterp.size();
    size_t numStates = xIntegForInterp[0].size();
    std::vector<std::vector<real>> coeffs(numStates,
                                          std::vector<real>(tLen, 0.0));
    static std::vector<std::vector<real>> coeffsPrev(
        numStates, std::vector<real>(tLen, 0.0));

    get_coeffs(tVecForInterp, xIntegForInterp, coeffs);
    // get_coeffs(tVecForInterpPrev, xIntegForInterpPrev, coeffsPrev);
    size_t &interpIdx = propSim->interpIdx;
    if (interpIdx == 0) {
        coeffsPrev = coeffs;
        tVecForInterpPrev = tVecForInterp;
    }
    const bool forwardIntegrate = tVecForInterp[0] < tVecForInterp[tLen - 1];
    const bool backwardIntegrate = tVecForInterp[0] > tVecForInterp[tLen - 1];
    bool interpIdxInWindow;
    get_interpIdxInWindow(propSim, tVecForInterp[0], tNext, forwardIntegrate,
                          backwardIntegrate, interpIdxInWindow);
    while (interpIdx < propSim->tEval.size() && interpIdxInWindow) {
        real tInterpGeom;
        if (!propSim->evalApparentState) {
            for (size_t i = 0; i < tLen; i++) {
                tInterpGeom = propSim->tEval[interpIdx];
                if (tInterpGeom == tVecForInterp[i]) {
                    propSim->xIntegEval.push_back(xIntegForInterp[i]);
                    interpIdx++;
                    continue;
                }
            }
        }
        if (propSim->tEvalUTC) {
            SpiceDouble etMinusUtc;
            real secPastJ2000Utc;
            mjd_to_et(propSim->tEval[interpIdx], secPastJ2000Utc);
            deltet_c(secPastJ2000Utc, "UTC", &etMinusUtc);
            tInterpGeom = et_to_mjd(secPastJ2000Utc + etMinusUtc);
        } else {
            tInterpGeom = propSim->tEval[interpIdx];
        }
        std::vector<real> xInterpGeom(numStates, 0.0);
        evaluate_one_interpolation(propSim, tInterpGeom, tVecForInterp, coeffs,
                                   tVecForInterpPrev, coeffsPrev, xInterpGeom);
        if (propSim->evalApparentState) {
            std::vector<real> lightTime(propSim->integParams.nInteg, 0.0);
            std::vector<real> xInterpApparent(numStates, 0.0);
            std::vector<real> radarMeasurement(
                propSim->integParams.nInteg,
                std::numeric_limits<real>::quiet_NaN());
            get_lightTime_and_xRelative(interpIdx, tInterpGeom, xInterpGeom,
                                        tVecForInterp, coeffs,
                                        tVecForInterpPrev, coeffsPrev, propSim,
                                        lightTime, xInterpApparent);
            if (propSim->radarObserver[interpIdx] == 1 ||
                propSim->radarObserver[interpIdx] == 2) {
                radarMeasurement =
                    std::vector<real>(propSim->integParams.nInteg, 0.0);
                get_radar_measurement(interpIdx, tInterpGeom, xInterpGeom,
                                      tVecForInterp, coeffs, tVecForInterpPrev,
                                      coeffsPrev, propSim, radarMeasurement);
            }
            propSim->lightTimeEval.push_back(lightTime);
            propSim->xIntegEval.push_back(xInterpApparent);
            propSim->radarObsEval.push_back(radarMeasurement);
        } else {
            propSim->xIntegEval.push_back(xInterpGeom);
        }
        interpIdx++;
        get_interpIdxInWindow(propSim, tVecForInterp[0], tNext,
                              forwardIntegrate, backwardIntegrate,
                              interpIdxInWindow);
    }
    coeffsPrev = coeffs;
    if (interpIdx == propSim->tEval.size()) {
        interpIdx = 0;
        tVecForInterpPrev.clear();
        xIntegForInterpPrev.clear();
        coeffsPrev.clear();
    }
}

void get_lightTime_and_xRelative(
    const size_t interpIdx, const real tInterpGeom,
    const std::vector<real> &xInterpGeom,
    const std::vector<real> &tVecForInterp,
    const std::vector<std::vector<real>> &coeffs,
    const std::vector<real> &tVecForInterpPrev,
    const std::vector<std::vector<real>> &coeffsPrev,
    const propSimulation *propSim, std::vector<real> &lightTime,
    std::vector<real> &xInterpApparent) {
    size_t numStates = xInterpGeom.size();
    std::vector<real> xObserver = propSim->xObserver[interpIdx];
    bool bouncePointAtLeadingEdge = false;
    if (propSim->observerInfo[interpIdx].size() == 9 ||
        propSim->observerInfo[interpIdx].size() == 10) {
        bouncePointAtLeadingEdge = propSim->observerInfo[interpIdx][8] == 1.0;
    }
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        real lightTimeTemp;
        std::vector<real> xInterpApparentTemp(numStates, 0.0);
        std::vector<real> xInterpApparentBary(6, 0.0);
        get_lightTimeOneBody(i, tInterpGeom, xInterpGeom, xObserver,
                             bouncePointAtLeadingEdge, tVecForInterp, coeffs,
                             tVecForInterpPrev, coeffsPrev, propSim,
                             lightTimeTemp);
        evaluate_one_interpolation(propSim, tInterpGeom - lightTimeTemp,
                                   tVecForInterp, coeffs, tVecForInterpPrev,
                                   coeffsPrev, xInterpApparentTemp);
        lightTime[i] = lightTimeTemp;
        for (size_t j = 0; j < 6; j++) {
            xInterpApparentBary[j] = xInterpApparentTemp[6 * i + j];
        }
        get_glb_correction(propSim, tInterpGeom, xInterpApparentBary);
        for (size_t j = 0; j < 6; j++) {
            xInterpApparent[6 * i + j] = xInterpApparentBary[j] - xObserver[j];
        }
    }
}

void get_lightTimeOneBody(const size_t &i, const real tInterpGeom,
                          std::vector<real> xInterpGeom,
                          std::vector<real> xObserver,
                          const bool bouncePointAtLeadingEdge,
                          const std::vector<real> &tVecForInterp,
                          const std::vector<std::vector<real>> &coeffs,
                          const std::vector<real> &tVecForInterpPrev,
                          const std::vector<std::vector<real>> &coeffsPrev,
                          const propSimulation *propSim,
                          real &lightTimeOneBody) {
    size_t numStates = xInterpGeom.size();
    std::vector<real> xInterpApparentFull(numStates, 0.0);
    std::vector<real> xInterpApparent(6, 0.0);
    std::vector<real> xRelativeOneBody(6, 0.0);
    real distRelativeOneBody;

    for (size_t j = 0; j < 6; j++) {
        xRelativeOneBody[j] = xInterpGeom[6 * i + j] - xObserver[j];
    }
    vnorm({xRelativeOneBody[0], xRelativeOneBody[1], xRelativeOneBody[2]},
          distRelativeOneBody);
    if (bouncePointAtLeadingEdge) {
        distRelativeOneBody -= propSim->forceParams.radii[i];
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
            evaluate_one_interpolation(propSim, tInterpGeom - lightTimeOneBody,
                                       tVecForInterp, coeffs, tVecForInterpPrev,
                                       coeffsPrev, xInterpApparentFull);
            for (size_t j = 0; j < 6; j++) {
                xRelativeOneBody[j] =
                    xInterpApparentFull[6 * i + j] - xObserver[j];
                xInterpApparent[j] = xInterpApparentFull[6 * i + j];
            }
            vnorm(
                {xRelativeOneBody[0], xRelativeOneBody[1], xRelativeOneBody[2]},
                distRelativeOneBody);
            if (bouncePointAtLeadingEdge) {
                distRelativeOneBody -= propSim->forceParams.radii[i];
            }
            lightTimeOneBodyPrev = lightTimeOneBody;
            get_delta_delay_relativistic(
                propSim, tInterpGeom - lightTimeOneBody, xInterpApparent,
                propSim->consts, deltaLightTimeRelativistic);
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

void get_glb_correction(const propSimulation *propSim, const real &tInterpGeom,
                        std::vector<real> &xInterpApparentBary) {
    Constants consts = propSim->consts;
    SpiceDouble sunState[6];
    SpiceDouble sunLightTime;
    get_spice_state_lt(10, tInterpGeom, consts, sunState, sunLightTime);
    SpiceDouble earthState[6];
    SpiceDouble earthLightTime;
    get_spice_state_lt(399, tInterpGeom, consts, earthState, earthLightTime);

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

    real G = consts.G;
    real sunGM = 0;
    for (size_t i = propSim->integParams.nInteg;
         i < propSim->integParams.nTotal; i++) {
        if (propSim->forceParams.spiceIdList[i] == 10) {
            sunGM = G * propSim->forceParams.masses[i];
        }
    }
    if (sunGM == 0) {
        throw std::runtime_error(
            "Sun GM not found in get_delta_delay_relativistic");
    }
    real c = consts.clight;

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

void get_radar_measurement(const size_t interpIdx, const real tInterpGeom,
                           const std::vector<real> &xInterpGeom,
                           const std::vector<real> &tVecForInterp,
                           const std::vector<std::vector<real>> &coeffs,
                           const std::vector<real> &tVecForInterpPrev,
                           const std::vector<std::vector<real>> &coeffsPrev,
                           const propSimulation *propSim,
                           std::vector<real> &radarMeasurement) {
    if (propSim->observerInfo[interpIdx].size() != 9 &&
        propSim->observerInfo[interpIdx].size() != 10) {
        throw std::runtime_error(
            "Error: observerInfo must be a 9 or 10 "
            "element vector for radar measurements");
    }
    size_t numStates = xInterpGeom.size();
    real receiveTimeTDB = tInterpGeom;
    std::vector<real> xTrgtBaryRcv = xInterpGeom;
    std::vector<real> xObsBaryRcv(6, 0.0);
    std::vector<real> receiverInfo = {propSim->observerInfo[interpIdx][0],
                                      propSim->observerInfo[interpIdx][1],
                                      propSim->observerInfo[interpIdx][2],
                                      propSim->observerInfo[interpIdx][3]};
    get_observer_state(receiveTimeTDB, receiverInfo, propSim->consts, false,
                       xObsBaryRcv);

    std::vector<real> xObsBaryTx(6, 0.0);
    std::vector<real> transmitterInfo = {propSim->observerInfo[interpIdx][4],
                                         propSim->observerInfo[interpIdx][5],
                                         propSim->observerInfo[interpIdx][6],
                                         propSim->observerInfo[interpIdx][7]};

    bool bouncePointAtLeadingEdge = propSim->observerInfo[interpIdx][8] == 1.0;
    real transmitFreq = 0.0;
    if (propSim->radarObserver[interpIdx] == 2) {
        transmitFreq = propSim->observerInfo[interpIdx][9];
    }
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        real delayDownleg;
        real bounceTimeTDB;
        real delayUpleg;
        real transmitTimeTDB;
        real delayMeasurement;
        std::vector<real> xTrgtBaryBounceAllBody(numStates, 0.0);
        std::vector<real> xTrgtBaryBounce(6, 0.0);
        // get downleg delay
        get_lightTimeOneBody(i, receiveTimeTDB, xTrgtBaryRcv, xObsBaryRcv,
                             bouncePointAtLeadingEdge, tVecForInterp, coeffs,
                             tVecForInterpPrev, coeffsPrev, propSim,
                             delayDownleg);
        bounceTimeTDB = receiveTimeTDB - delayDownleg;
        evaluate_one_interpolation(propSim, bounceTimeTDB, tVecForInterp,
                                   coeffs, tVecForInterpPrev, coeffsPrev,
                                   xTrgtBaryBounceAllBody);
        for (size_t j = 0; j < 6; j++) {
            xTrgtBaryBounce[j] = xTrgtBaryBounceAllBody[6 * i + j];
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
                                   propSim->consts, false, xObsBaryTx);
                std::vector<real> xRelativeOneBody(6, 0.0);
                for (size_t j = 0; j < 6; j++) {
                    xRelativeOneBody[j] = xTrgtBaryBounce[j] - xObsBaryTx[j];
                }
                vnorm({xRelativeOneBody[0], xRelativeOneBody[1],
                       xRelativeOneBody[2]},
                      distRelativeUpleg);
                if (bouncePointAtLeadingEdge) {
                    distRelativeUpleg -= propSim->forceParams.radii[i];
                }
                delayUplegPrev = delayUpleg;
                get_delta_delay_relativistic(
                    propSim, bounceTimeTDB - delayUpleg, xTrgtBaryBounce,
                    propSim->consts, deltaDelayUplegRelativistic);
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
        get_observer_state(transmitTimeTDB, transmitterInfo, propSim->consts,
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
        if (propSim->radarObserver[interpIdx] == 1) {
            radarMeasurement[i] = delayMeasurement;
        } else if (propSim->radarObserver[interpIdx] == 2) {
            real dopplerMeasurement;
            get_doppler_measurement(propSim, receiveTimeTDB, transmitTimeTDB,
                                    xObsBaryRcv, xTrgtBaryBounce, xObsBaryTx,
                                    transmitFreq, dopplerMeasurement);
            radarMeasurement[i] = dopplerMeasurement;
        }
    }
}

void get_delta_delay_relativistic(const propSimulation *propSim,
                                  const real &tForSpice,
                                  const std::vector<real> &targetState,
                                  const Constants &consts,
                                  real &deltaDelayRelativistic) {
    // from Standish (1990),
    // https://ui.adsabs.harvard.edu/abs/1990A&A...233..252S
    SpiceDouble sunState[6];
    SpiceDouble sunLightTime;
    get_spice_state_lt(10, tForSpice, consts, sunState, sunLightTime);
    SpiceDouble earthState[6];
    SpiceDouble earthLightTime;
    get_spice_state_lt(399, tForSpice, consts, earthState, earthLightTime);

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

    real G = consts.G;
    real sunGM = 0;
    for (size_t i = propSim->integParams.nInteg;
         i < propSim->integParams.nTotal; i++) {
        if (propSim->forceParams.spiceIdList[i] == 10) {
            sunGM = G * propSim->forceParams.masses[i];
        }
    }
    if (sunGM == 0) {
        throw std::runtime_error(
            "Sun GM not found in get_delta_delay_relativistic");
    }
    real c = consts.clight;
    real gamma = 1.0L;  // PPN parameter

    deltaDelayRelativistic = (1 + gamma) * sunGM * pow(c, -3) *
        log((sunEarthDist + sunTargetDist + earthTargetDist) /
            (sunEarthDist + sunTargetDist - earthTargetDist));
}

void get_doppler_measurement(
    const propSimulation *propSim, const real receiveTimeTDB,
    const real transmitTimeTDB, const std::vector<real> xObsBaryRcv,
    const std::vector<real> xTrgtBaryBounce, const std::vector<real> xObsBaryTx,
    const real transmitFreq, real &dopplerMeasurement) {
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

    SpiceDouble lt;
    SpiceDouble xSun3[6];
    get_spice_state_lt(10, receiveTimeTDB, propSim->consts, xSun3, lt);
    SpiceDouble xSun1[6];
    get_spice_state_lt(10, transmitTimeTDB, propSim->consts, xSun1, lt);

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
    for (size_t i = propSim->integParams.nInteg;
         i < propSim->integParams.nTotal; i++) {
        if (propSim->forceParams.spiceIdList[i] == 10) {
            sunGM = G * propSim->forceParams.masses[i];
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
}
