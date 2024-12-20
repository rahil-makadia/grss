/**
 * @file    observe.cpp
 * @brief   Source file for the observable computations.
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

#include "observe.h"

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] interpIdx Index of the next interpolation time.
 * @param[in] tInterpGeom Time to interpolate to.
 * @param[out] xInterpApparentBary Apparent state vector of the target body.
 */
void get_glb_correction(PropSimulation *propSim, const size_t &interpIdx,
                        const real &tInterpGeom,
                        std::vector<real> &xInterpApparentBary) {
    const real G = propSim->consts.G;
    real sunGM = 0;
    real earthGM = 0;
    for (size_t i = 0; i < propSim->integParams.nSpice; i++) {
        if (propSim->spiceBodies[i].spiceId == 10) {
            sunGM = G * propSim->spiceBodies[i].mass;
        }
        if (propSim->spiceBodies[i].spiceId == 399) {
            earthGM = G * propSim->spiceBodies[i].mass;
        }
    }
    if (sunGM == 0 || earthGM == 0) {
        throw std::runtime_error(
            "Sun GM or Earth GM not found in get_delta_delay_relativistic");
    }
    const real c = propSim->consts.clight;

    double sunState[9];
    double earthState[9];
    get_spk_state(10, tInterpGeom, propSim->spkEphem, sunState);
    get_spk_state(399, tInterpGeom, propSim->spkEphem, earthState);
    std::vector<real> observerState = propSim->xObserver[interpIdx];
    std::vector<real> observerTargetPos = {xInterpApparentBary[0] - observerState[0],
                                            xInterpApparentBary[1] - observerState[1],
                                            xInterpApparentBary[2] - observerState[2]};
    real observerTargetDist;
    vnorm(observerTargetPos, observerTargetDist);
    std::vector<real> p(3, 0.0);
    vunit(observerTargetPos, p);

    std::vector<real> p1(p);
    size_t numBendingBodies = 1; // sun (1) or sun+earth (2)
    for (size_t i = 0; i < numBendingBodies; i++){
        double centralBodyState[6];
        real centralBodyGM;
        switch (i) {
        case 0:
            for (size_t j = 0; j < 6; j++) {
                centralBodyState[j] = sunState[j];
            }
            centralBodyGM = sunGM;
            break;
        case 1:
            for (size_t j = 0; j < 6; j++) {
                centralBodyState[j] = earthState[j];
            }
            centralBodyGM = earthGM;
            break;
        default:
            throw std::runtime_error(
                "get_glb_correction: central body index must be 0 or 1");
            break;
        }

        std::vector<real> centralBodyObserverPos = {observerState[0] - centralBodyState[0],
                                                    observerState[1] - centralBodyState[1],
                                                    observerState[2] - centralBodyState[2]};
        real centralBodyObserverDist;
        vnorm(centralBodyObserverPos, centralBodyObserverDist);
        std::vector<real> centralBodyTargetPos = {xInterpApparentBary[0] - centralBodyState[0],
                                                    xInterpApparentBary[1] - centralBodyState[1],
                                                    xInterpApparentBary[2] - centralBodyState[2]};
        real centralBodyTargetDist;
        vnorm(centralBodyTargetPos, centralBodyTargetDist);
        
        // from section 7.2.4 in the Explanatory Supplement to the Astronomical
        // Almanac, 3rd edition
        std::vector<real> e(3, 0.0);
        vunit(centralBodyObserverPos, e);
        std::vector<real> q(3, 0.0);
        vunit(centralBodyTargetPos, q);
        
        real pDotQ, eDotP, qDotE;
        vdot(p, q, pDotQ);
        vdot(e, p, eDotP);
        vdot(q, e, qDotE);

        const real g1 = 2 * centralBodyGM / c / c / centralBodyObserverDist;
        const real g2 = 1.0L + qDotE;

        std::vector<real> deltaP1Targ(3, 0.0);
        deltaP1Targ[0] += g1 * (pDotQ * e[0] - eDotP * q[0]) / g2;
        deltaP1Targ[1] += g1 * (pDotQ * e[1] - eDotP * q[1]) / g2;
        deltaP1Targ[2] += g1 * (pDotQ * e[2] - eDotP * q[2]) / g2;

        std::vector<real> deltaP1Star(3, 0.0);
        deltaP1Star[0] += g1 * (e[0] - eDotP * p[0]) / (1 + eDotP);
        deltaP1Star[1] += g1 * (e[1] - eDotP * p[1]) / (1 + eDotP);
        deltaP1Star[2] += g1 * (e[2] - eDotP * p[2]) / (1 + eDotP);

        // do absolute correction first
        p1[0] += deltaP1Targ[0];
        p1[1] += deltaP1Targ[1];
        p1[2] += deltaP1Targ[2];
        // if not gaia obs, do full relative correction
        if (propSim->obsType[interpIdx] != 3) {
            p1[0] -= deltaP1Star[0];
            p1[1] -= deltaP1Star[1];
            p1[2] -= deltaP1Star[2];
        }
        // re-normalize
        const real p1Norm = sqrt(p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2]);
        for (size_t j = 0; j < 3; j++) {
            p1[j] /= p1Norm;
        }
    }

    observerTargetPos[0] = observerTargetDist * p1[0];
    observerTargetPos[1] = observerTargetDist * p1[1];
    observerTargetPos[2] = observerTargetDist * p1[2];

    xInterpApparentBary[0] = observerState[0] + observerTargetPos[0];
    xInterpApparentBary[1] = observerState[1] + observerTargetPos[1];
    xInterpApparentBary[2] = observerState[2] + observerTargetPos[2];
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] interpIdx Index of the next interpolation time.
 * @param[in] tInterpGeom Time to interpolate to.
 * @param[in] xInterpGeom Geometric state vector of the target body at the interpolation time.
 * @param[in] xInterpApparent Apparent state vector of the target body at the interpolation time.
 */
void get_measurement(PropSimulation *propSim, const size_t &interpIdx,
                     const real tInterpGeom,
                     const std::vector<real> &xInterpGeom,
                     const std::vector<real> &xInterpApparent) {
    std::vector<real> opticalMeasurement(2*propSim->integParams.nInteg,
                                       std::numeric_limits<real>::quiet_NaN());
    std::vector<real> opticalMeasurementDot(2*propSim->integParams.nInteg,
                                       std::numeric_limits<real>::quiet_NaN());
    std::vector<real> opticalPartials(12*propSim->integParams.nInteg,
                                       std::numeric_limits<real>::quiet_NaN());
    std::vector<real> photocenterCorr(2*propSim->integParams.nInteg,
                                       std::numeric_limits<real>::quiet_NaN());
    std::vector<real> radarMeasurement(propSim->integParams.nInteg,
                                       std::numeric_limits<real>::quiet_NaN());
    std::vector<real> radarPartials(6*propSim->integParams.nInteg,
                                       std::numeric_limits<real>::quiet_NaN());
    switch (propSim->obsType[interpIdx]) {
    case 0: case 3:
        get_optical_measurement(propSim, xInterpApparent, opticalMeasurement,
                                opticalMeasurementDot, opticalPartials);
        get_photocenter_correction(propSim, interpIdx, tInterpGeom,
                                    xInterpApparent, photocenterCorr);
        break;
    case 1: case 2:
        get_radar_measurement(propSim, interpIdx, tInterpGeom,
                                xInterpGeom, radarMeasurement, radarPartials);
        break;
    default:
        throw std::runtime_error(
            "get_measurement: obsType flag must be 0, 1, 2, or 3");
        break;
    }
    propSim->opticalObs.push_back(opticalMeasurement);
    propSim->opticalObsDot.push_back(opticalMeasurementDot);
    propSim->opticalPartials.push_back(opticalPartials);
    propSim->opticalObsCorr.push_back(photocenterCorr);
    propSim->radarObs.push_back(radarMeasurement);
    propSim->radarPartials.push_back(radarPartials);
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] xInterpApparent Apparent state vector of the target body.
 * @param[out] opticalMeasurement Optical measurement (RA and Dec).
 * @param[out] opticalMeasurementDot Time derivative of the optical measurement.
 * @param[out] opticalPartials Partials of the optical measurement.
 */
void get_optical_measurement(PropSimulation *propSim,
                             const std::vector<real> &xInterpApparent,
                             std::vector<real> &opticalMeasurement,
                             std::vector<real> &opticalMeasurementDot,
                             std::vector<real> &opticalPartials) {
    size_t starti = 0;
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        real rho[3], rhoHat[3], rhoDot[3];
        for (size_t j = 0; j < 3; j++) {
            rho[j] = xInterpApparent[starti + j];
            rhoDot[j] = xInterpApparent[starti + 3 + j];
        }
        real dist;
        vnorm(rho, 3, dist);
        for (size_t j = 0; j < 3; j++) {
            rhoHat[j] = rho[j]/dist;
        }
        real r_asc = atan2(rhoHat[1], rhoHat[0]);
        if (r_asc < 0) {
            r_asc = r_asc + 2*PI;
        }
        const real x2py2 = rho[0]*rho[0] + rho[1]*rho[1];
        const real xydist = sqrt(x2py2);
        const real dradx = -rho[1]/x2py2;
        const real drady = rho[0]/x2py2;
        real dec = asin(rhoHat[2]);
        const real dist2 = dist*dist;
        const real ddecdx = -rho[0]*rho[2]/dist2/xydist;
        const real ddecdy = -rho[1]*rho[2]/dist2/xydist;
        const real ddecdz = xydist/dist2;
        const real conv = 180.0L/PI*3600.0L; // radians -> arcsec
        opticalMeasurement[2*i] = r_asc*conv;
        opticalMeasurement[2*i+1] = dec*conv;
        const real rdot = (rhoDot[0]*rho[0] + rhoDot[1]*rho[1] + rhoDot[2]*rho[2])/dist;
        const real raDot = (rhoDot[0]*rho[1] - rhoDot[1]*rho[0])/-x2py2;
        const real decDot = (rhoDot[2] - rdot*sin(dec))/xydist;
        opticalMeasurementDot[2*i] = raDot*conv;
        opticalMeasurementDot[2*i+1] = decDot*conv;
        std::fill(opticalPartials.begin()+12*i, opticalPartials.begin()+12*(i+1), 0.0);
        opticalPartials[12*i] = dradx*conv;
        opticalPartials[12*i+1] = drady*conv;
        opticalPartials[12*i+6] = ddecdx*conv;
        opticalPartials[12*i+7] = ddecdy*conv;
        opticalPartials[12*i+8] = ddecdz*conv;
        starti += 2*propSim->integBodies[i].n2Derivs;
    }
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] interpIdx Index of the next interpolation time.
 * @param[in] tInterpGeom Time to interpolate to.
 * @param[in] xInterpApparent Apparent state vector of the target body.
 * @param[out] photocenterCorr Photocenter-barycenter correction to the optical measurement.
 */
void get_photocenter_correction(PropSimulation *propSim, const size_t &interpIdx,
                                const real &tInterpGeom,
                                const std::vector<real> &xInterpApparent,
                                std::vector<real> &photocenterCorr){
    const size_t obsType = propSim->obsType[interpIdx];
    if (obsType == 1 || obsType == 2){
        return;
    }
    if (obsType == 3) {
        // // polynomial coefficients from JPL paper (fuentes-munoz et al. 2024)
        // std::vector<real> polyCoeffs = {
        //     -0.02384,
        //     0.05579,
        //     0.329,
        //     0};
        // polynomial coefficients from GRSS fit (backyard/random/grss_corr.ipynb)
        std::vector<real> polyCoeffs = {
            -0.02352667223191772,
            0.05403589930067203,
            0.3318343597581827,
            -0.001233060221632019
        };
        size_t starti = 0;
        std::vector<real> xObserver = propSim->xObserver[interpIdx];
        for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
            const real radius = propSim->integBodies[i].radius;
            if (radius == 0.0){
                photocenterCorr[2*i] = 0.0;
                photocenterCorr[2*i+1] = 0.0;
                continue;
            }
            std::vector<real> xInterpApparentBaryOneBody(6, 0.0);
            for (size_t j = 0; j < 6; j++) {
                xInterpApparentBaryOneBody[j] = xInterpApparent[starti+j] + xObserver[j];
            }
            const real tForSpice = tInterpGeom - propSim->lightTimeEval[interpIdx][i];
            double sunState[9];
            get_spk_state(10, tForSpice, propSim->spkEphem, sunState);
            real rhoVec[3], rhoHat[3], rho, rVec[3], rHat[3];
            for (size_t j = 0; j < 3; j++) {
                rhoVec[j] = xInterpApparent[starti+j];
                rVec[j] = xInterpApparentBaryOneBody[j] - sunState[j];
            }
            vnorm(rhoVec, 3, rho);
            for (size_t j = 0; j < 3; j++) {
                rhoHat[j] = rhoVec[j]/rho;
            }
            vunit(rVec, 3, rHat);
            real rHatDotRhoHat;
            vdot(rHat, rhoHat, 3, rHatDotRhoHat);
            const real alpha = acos(rHatDotRhoHat);
            real fval = 0.0;
            for (size_t j = 0; j < polyCoeffs.size(); j++) {
                fval = fval*alpha + polyCoeffs[j];
            }
            fval *= radius;
            real tHat[3], tVec[3];
            for (size_t j = 0; j < 3; j++) {
                tHat[j] = (rHatDotRhoHat*rhoHat[j] - rHat[j])/sin(alpha);
                tVec[j] = fval/rho*tHat[j];
            }
            const real conv = 180.0L/PI*3600.0L; // radians -> arcsec
            const real raHatGetter[3] = {0, 0, 1};
            real raDir[3], raHat[3], decHat[3];
            vcross(raHatGetter, rhoHat, raDir);
            vunit(raDir, 3, raHat);
            vcross(rhoHat, raHat, decHat);
            real raCosDecCorrection, decCorrection;
            vdot(tVec, raHat, 3, raCosDecCorrection);
            vdot(tVec, decHat, 3, decCorrection);
            photocenterCorr[2*i] = raCosDecCorrection*conv;
            photocenterCorr[2*i+1] = decCorrection*conv;
            starti += 2*propSim->integBodies[i].n2Derivs;
        }
    } else if (obsType == 0) {
        std::fill(photocenterCorr.begin(), photocenterCorr.end(), 0.0);
    } else {
        throw std::runtime_error(
            "get_photocenter_correction: obsType flag must be 0, 1, 2, or 3");
    }
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] interpIdx Index of the next interpolation time.
 * @param[in] tInterpGeom Time to interpolate to.
 * @param[in] xInterpGeom Geometric state vector of the target body at the interpolation time.
 * @param[out] radarMeasurement Radar measurement (range/Doppler).
 * @param[out] radarPartials Partials of the radar measurement.
 */
void get_radar_measurement(PropSimulation *propSim, const size_t &interpIdx,
                           const real tInterpGeom,
                           const std::vector<real> &xInterpGeom,
                           std::vector<real> &radarMeasurement,
                           std::vector<real> &radarPartials) {
    if (propSim->obsType[interpIdx] != 1 && propSim->obsType[interpIdx] != 2) {
        throw std::runtime_error(
            "get_radar_measurement: obsType must be 1 or 2 for radar "
            "measurements");
    }
    real receiveTimeTDB = tInterpGeom;
    real transmitTimeTDB;
    std::vector<real> xTrgtBaryRcv = xInterpGeom;
    std::vector<real> xObsBaryRcv(6, 0.0);
    std::vector<real> xTrgtBaryBounce(6, 0.0);
    std::vector<real> xObsBaryTx(6, 0.0);
    real transmitFreq = 0.0;
    if (propSim->obsType[interpIdx] == 2) {
        transmitFreq = propSim->observerInfo[interpIdx][9]*1.0e6L;
    }
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        std::fill(radarPartials.begin()+6*i, radarPartials.begin()+6*(i+1), 0.0);
        real delayMeasurement;
        get_delay_measurement(propSim, interpIdx, i, tInterpGeom,
                              xInterpGeom, receiveTimeTDB, transmitTimeTDB,
                              xObsBaryRcv, xTrgtBaryBounce, xObsBaryTx,
                              delayMeasurement, radarPartials);
        if (propSim->obsType[interpIdx] == 1) {
            radarMeasurement[i] = delayMeasurement;
        } else if (propSim->obsType[interpIdx] == 2) {
            real dopplerMeasurement;
            get_doppler_measurement(propSim, i, receiveTimeTDB, transmitTimeTDB,
                                    xObsBaryRcv, xTrgtBaryBounce, xObsBaryTx,
                                    transmitFreq, dopplerMeasurement, radarPartials);
            radarMeasurement[i] = dopplerMeasurement;
        }
    } 
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] interpIdx Index of the next interpolation time.
 * @param[in] i Index of the target body.
 * @param[in] tInterpGeom Time to interpolate to.
 * @param[in] xInterpGeom Geometric state vector of the target body at the interpolation time.
 * @param[in] receiveTimeTDB Time of reception of the radar signal (TDB).
 * @param[out] transmitTimeTDB Time of transmission of the radar signal (TDB).
 * @param[out] xObsBaryRcv Barycentric state vector of the observer at the reception time.
 * @param[out] xTrgtBaryBounce Barycentric state vector of the target body at the bounce time.
 * @param[out] xObsBaryTx Barycentric state vector of the observer at the transmission time.
 * @param[out] delayMeasurement Radar delay measurement.
 * @param[out] delayPartials Partials of the radar delay measurement.
 */
void get_delay_measurement(PropSimulation *propSim, const size_t &interpIdx,
                           const size_t &i, const real tInterpGeom,
                           const std::vector<real> &xInterpGeom,
                           const real &receiveTimeTDB, real &transmitTimeTDB,
                           std::vector<real> &xObsBaryRcv,
                           std::vector<real> &xTrgtBaryBounce,
                           std::vector<real> &xObsBaryTx,
                           real &delayMeasurement,
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
    bool bouncePointAtCenterOfMass = propSim->observerInfo[interpIdx][8] == 1.0;
    real delayDownleg;
    real bounceTimeTDB;
    real delayUpleg;
    std::vector<real> xTrgtBaryBounceAllBody(numStates, 0.0);
    // downleg delay is the already evaluated light time
    delayDownleg = propSim->lightTimeEval[interpIdx][i];
    bounceTimeTDB = receiveTimeTDB - delayDownleg;
    evaluate_one_interpolation(propSim, bounceTimeTDB,
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
            if (!bouncePointAtCenterOfMass) {
                distRelativeUpleg -= propSim->integBodies[i].radius;
            }
            delayUplegPrev = delayUpleg;
            delayUpleg = distRelativeUpleg / propSim->consts.clight;
            iter++;
        }
        real deltaDelayUplegRelativistic;
        get_delta_delay_relativistic(
            propSim, bounceTimeTDB - delayUpleg, xTrgtBaryBounce,
            deltaDelayUplegRelativistic);
        delayUpleg += deltaDelayUplegRelativistic;
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
        const real etMinusUtcReceiveTime = delta_et_tdb(receiveTimeTDB);
        const real etMinusUtcTransmitTime = delta_et_tdb(transmitTimeTDB);
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

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] tForSpice Time for querying SPICE Ephemeris for the Sun and Earth positions.
 * @param[in] targetState State vector of the target body.
 * @param[out] deltaDelayRelativistic Relativistic delay correction (Shapiro delay).
 */
void get_delta_delay_relativistic(PropSimulation *propSim,
                                  const real &tForSpice,
                                  const std::vector<real> &targetState,
                                  real &deltaDelayRelativistic) {
    // from Standish (1990),
    // https://ui.adsabs.harvard.edu/abs/1990A&A...233..252S
    double sunState[9];
    double earthState[9];
    get_spk_state(10, tForSpice, propSim->spkEphem, sunState);
    get_spk_state(399, tForSpice, propSim->spkEphem, earthState);

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

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] i Index of the target body.
 * @param[in] receiveTimeTDB Time of reception of the radar signal (TDB).
 * @param[in] transmitTimeTDB Time of transmission of the radar signal (TDB).
 * @param[in] xObsBaryRcv Barycentric state vector of the observer at the reception time.
 * @param[in] xTrgtBaryBounce Barycentric state vector of the target body at the bounce time.
 * @param[in] xObsBaryTx Barycentric state vector of the observer at the transmission time.
 * @param[in] transmitFreq Frequency of the transmitted radar signal.
 * @param[out] dopplerMeasurement Doppler measurement.
 * @param[out] dopplerPartials Partials of the Doppler measurement.
 */
void get_doppler_measurement(PropSimulation *propSim, const size_t &i,
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
    get_spk_state(10, receiveTimeTDB, propSim->spkEphem, xSun3);
    get_spk_state(10, transmitTimeTDB, propSim->spkEphem, xSun1);

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
