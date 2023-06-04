#include "interpolate.h"

void interpolate(const real &t, const real &dt, const std::vector<real> &xInteg0, const std::vector<real> &accInteg0, const std::vector< std::vector<real> > &b, propSimulation *propSim){
    size_t nh = 8;
    std::vector<real> tVecForInterp (nh, 0.0);
    std::vector< std::vector<real> > xIntegForInterp(nh, std::vector<real>(xInteg0.size(), 0.0));
    static std::vector<real> tVecForInterpPrev(tVecForInterp.size(), 0.0);
    static std::vector< std::vector<real> > xIntegForInterpPrev(xIntegForInterp.size(), std::vector<real>(xIntegForInterp[0].size(), 0.0));
    
    tVecForInterp[0] = t;
    xIntegForInterp[0] = xInteg0;
    for (size_t hIdx = 1; hIdx < nh; hIdx++) {
        tVecForInterp[hIdx] = t + hVec[hIdx]*dt;
        approx_xInteg(xInteg0, accInteg0, xIntegForInterp[hIdx], dt, hVec[hIdx], b, propSim->integParams.nInteg);
    }
    // tVecForInterp.push_back(t+dt);
    // xIntegForInterp.push_back(xInteg);

    one_timestep_interpolation(t+dt, tVecForInterp, xIntegForInterp, tVecForInterpPrev, xIntegForInterpPrev, propSim);
    tVecForInterpPrev = tVecForInterp;
    xIntegForInterpPrev = xIntegForInterp;
}

void get_coeffs(const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &xIntegForInterp, std::vector< std::vector<real> > &coeffs){
    size_t tLen = tVecForInterp.size();
    size_t numStates = xIntegForInterp[0].size();
    std::vector< std::vector< std::vector<real> > > c(numStates, std::vector< std::vector<real> >(tLen, std::vector<real>(tLen, 0.0)));
    for (size_t i = 0; i < numStates; i++){
        for (size_t j = 0; j < tLen; j++){
            c[i][j][0] = xIntegForInterp[j][i];
        }
        for (size_t j = 1; j < tLen; j++){
            for (size_t k = 0; k < tLen-j; k++){
                c[i][k][j] = (c[i][k+1][j-1] - c[i][k][j-1])/(tVecForInterp[k+j]-tVecForInterp[k]);
            }
        }
    }

    for (size_t i = 0; i < numStates; i++){
        for (size_t j = 0; j < tLen; j++){
            coeffs[i][j] = c[i][0][j];
        }
    }
}

void evaluate_one_interpolation(const real &tInterp, const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &coeffs, std::vector<real> &xInterp){
    size_t tLen = tVecForInterp.size();
    size_t numStates = xInterp.size();
    size_t n = tLen-1;
    for (size_t i = 0; i < numStates; i++){
        xInterp[i] = coeffs[i][n];
    }
    for (size_t i = 1; i < n+1; i++){
        for (size_t j = 0; j < numStates; j++){
            xInterp[j] = coeffs[j][n-i] + (tInterp-tVecForInterp[n-i])*xInterp[j];
        }
    }
}

void one_timestep_interpolation(const real &tNext, const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &xIntegForInterp, std::vector<real> &tVecForInterpPrev, std::vector< std::vector<real> > &xIntegForInterpPrev, propSimulation *propSim){
    size_t tLen = tVecForInterp.size();
    size_t numStates = xIntegForInterp[0].size();
    std::vector< std::vector<real> > coeffs(numStates, std::vector<real>(tLen, 0.0));
    static std::vector< std::vector<real> > coeffsPrev(numStates, std::vector<real>(tLen, 0.0));

    get_coeffs(tVecForInterp, xIntegForInterp, coeffs);
    // get_coeffs(tVecForInterpPrev, xIntegForInterpPrev, coeffsPrev);
    static size_t interpIdx = 0;
    bool forwardIntegrate = tVecForInterp[0] < tVecForInterp[tLen-1];
    bool backwardIntegrate = tVecForInterp[0] > tVecForInterp[tLen-1];
    while ( interpIdx < propSim->tEval.size()
            &&( (forwardIntegrate && (propSim->tEval[interpIdx] == tVecForInterp[0] || (propSim->tEval[interpIdx] > tVecForInterp[0] && propSim->tEval[interpIdx] <= tNext)))
            ||  (forwardIntegrate && propSim->tEval[interpIdx] <= propSim->integParams.t0 && propSim->tEval[interpIdx]+propSim->tEvalMargin >= propSim->integParams.t0) || (forwardIntegrate && propSim->tEval[interpIdx] >= propSim->integParams.tf && propSim->tEval[interpIdx]-propSim->tEvalMargin <= propSim->integParams.tf)
            ||  (backwardIntegrate && (propSim->tEval[interpIdx] == tVecForInterp[0] || (propSim->tEval[interpIdx] < tVecForInterp[0] && propSim->tEval[interpIdx] >= tNext)))
            ||  (backwardIntegrate && propSim->tEval[interpIdx] >= propSim->integParams.t0 && propSim->tEval[interpIdx]-propSim->tEvalMargin <= propSim->integParams.t0) || (backwardIntegrate && propSim->tEval[interpIdx] <= propSim->integParams.tf && propSim->tEval[interpIdx]+propSim->tEvalMargin >= propSim->integParams.tf) )
        ){
        real tInterpGeom;
        if (!propSim->evalApparentState){
            for (size_t i = 0; i < tLen; i++){
                tInterpGeom = propSim->tEval[interpIdx];
                if (tInterpGeom == tVecForInterp[i]){
                    propSim->xIntegEval.push_back(xIntegForInterp[i]);
                    interpIdx++;
                    continue;
                }
            }
        }
        if (propSim->tEvalUTC){
            SpiceDouble etMinusUtc;
            real secPastJ2000Utc;
            mjd_to_et(propSim->tEval[interpIdx], secPastJ2000Utc);
            deltet_c (secPastJ2000Utc, "UTC", &etMinusUtc);
            tInterpGeom = et_to_mjd(secPastJ2000Utc + etMinusUtc);
        } else {
            tInterpGeom = propSim->tEval[interpIdx];
        }
        std::vector<real> xInterpGeom(numStates, 0.0);
        evaluate_one_interpolation(tInterpGeom, tVecForInterp, coeffs, xInterpGeom);
        if (propSim->evalApparentState){
            std::vector<real> lightTime(propSim->integParams.nInteg, 0.0);
            std::vector<real> xInterpApparent(numStates, 0.0);
            std::vector<real> radarMeasurement(propSim->integParams.nInteg, std::numeric_limits<real>::quiet_NaN());
            get_lightTime_and_xRelative(interpIdx, tInterpGeom, xInterpGeom, tVecForInterp, coeffs, tVecForInterpPrev, coeffsPrev, propSim, lightTime, xInterpApparent);
            if (propSim->radarObserver[interpIdx] == 1){
                radarMeasurement = std::vector<real>(propSim->integParams.nInteg, 0.0);
                get_radar_measurement(interpIdx, tInterpGeom, xInterpGeom, tVecForInterp, coeffs, tVecForInterpPrev, coeffsPrev, propSim, radarMeasurement);
            }
            else if (propSim->radarObserver[interpIdx] == 2){
                throw std::runtime_error("Doppler measurements not implemented yet. Will be added to get_radar_measurement()");
            }
            propSim->lightTimeEval.push_back(lightTime);
            propSim->xIntegEval.push_back(xInterpApparent);
            // std::cout << "tEval = " << propSim->tEval[interpIdx] << std::endl;
            propSim->radarObsEval.push_back(radarMeasurement);
        } else {
            propSim->xIntegEval.push_back(xInterpGeom);
        }
        interpIdx++;
    }
    coeffsPrev = coeffs;
    if (interpIdx == propSim->tEval.size()){
        interpIdx = 0;
        tVecForInterpPrev.clear();
        xIntegForInterpPrev.clear();
        coeffsPrev.clear();
    }
}

void get_lightTime_and_xRelative(const size_t interpIdx, const real tInterpGeom, const std::vector<real> &xInterpGeom, const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &coeffs, const std::vector<real> &tVecForInterpPrev, const std::vector< std::vector<real> > &coeffsPrev, const propSimulation *propSim, std::vector<real> &lightTime, std::vector<real> &xInterpApparent){
    size_t numStates = xInterpGeom.size();
    std::vector<real> xObserver = propSim->xObserver[interpIdx];
    bool bouncePointAtLeadingEdge = false;
    bool forwardIntegrate = tVecForInterp[0] < tVecForInterp[1];
    bool backwardIntegrate = tVecForInterp[0] > tVecForInterp[1];
    for (size_t i = 0; i < propSim->integParams.nInteg; i++){
        real lightTimeTemp;
        std::vector<real> xInterpApparentTemp(numStates, 0.0);
        get_lightTimeOneBody(i, tInterpGeom, xInterpGeom, xObserver, bouncePointAtLeadingEdge, tVecForInterp, coeffs, tVecForInterpPrev, coeffsPrev, propSim, lightTimeTemp);
        if ( (forwardIntegrate && (tInterpGeom-lightTimeTemp < tVecForInterp[0] && tInterpGeom != propSim->integParams.t0))
            || (backwardIntegrate && (tInterpGeom-lightTimeTemp > tVecForInterp[0] && tInterpGeom != propSim->integParams.t0)) ){
            evaluate_one_interpolation(tInterpGeom-lightTimeTemp, tVecForInterpPrev, coeffsPrev, xInterpApparentTemp);
        } else {
            evaluate_one_interpolation(tInterpGeom-lightTimeTemp, tVecForInterp, coeffs, xInterpApparentTemp);
        }
        lightTime[i] = lightTimeTemp;
        for (size_t j = 0; j < 6; j++){
            xInterpApparent[6*i+j] = xInterpApparentTemp[6*i+j] - xObserver[j];
        }
    }
}

void get_lightTimeOneBody(const size_t &i, const real tInterpGeom, std::vector<real> xInterpGeom, std::vector<real> xObserver, const bool bouncePointAtLeadingEdge, const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &coeffs, const std::vector<real> &tVecForInterpPrev, const std::vector< std::vector<real> > &coeffsPrev, const propSimulation *propSim, real &lightTimeOneBody){
    size_t numStates = xInterpGeom.size();
    std::vector<real> xInterpApparentTemp(numStates, 0.0);
    std::vector<real> xRelativeOneBody(6, 0.0);
    real distRelativeOneBody;

    for (size_t j = 0; j < 6; j++){
        xRelativeOneBody[j] = xInterpGeom[6*i+j] - xObserver[j];
    }
    vnorm({xRelativeOneBody[0], xRelativeOneBody[1], xRelativeOneBody[2]}, distRelativeOneBody);
    if (bouncePointAtLeadingEdge){
        distRelativeOneBody -= propSim->forceParams.radii[i];
    }
    lightTimeOneBody = distRelativeOneBody/propSim->consts.clight;
    bool forwardIntegrate = tVecForInterp[0] < tVecForInterp[1];
    bool backwardIntegrate = tVecForInterp[0] > tVecForInterp[1];
    if (propSim->convergedLightTime){
        real lightTimeTol = 1e-16/86400.0L;
        real lightTimeOneBodyPrev = 0.0L;
        size_t maxIter = 20;
        size_t iter = 0;
        // keep iterating until max iterations or light time tolerance is met
        while (iter < maxIter && fabs(lightTimeOneBody - lightTimeOneBodyPrev) > lightTimeTol){
            if ( (forwardIntegrate && (tInterpGeom-lightTimeOneBody < tVecForInterp[0] && tInterpGeom != propSim->integParams.t0))
                || (backwardIntegrate && (tInterpGeom-lightTimeOneBody > tVecForInterp[0] && tInterpGeom != propSim->integParams.t0)) ){
                evaluate_one_interpolation(tInterpGeom-lightTimeOneBody, tVecForInterpPrev, coeffsPrev, xInterpApparentTemp);
            } else {
                evaluate_one_interpolation(tInterpGeom-lightTimeOneBody, tVecForInterp, coeffs, xInterpApparentTemp);
            }
            for (size_t j = 0; j < 6; j++){
                xRelativeOneBody[j] = xInterpApparentTemp[6*i+j] - xObserver[j];
            }
            vnorm({xRelativeOneBody[0], xRelativeOneBody[1], xRelativeOneBody[2]}, distRelativeOneBody);
            if (bouncePointAtLeadingEdge){
                distRelativeOneBody -= propSim->forceParams.radii[i];
            }
            lightTimeOneBodyPrev = lightTimeOneBody;
            lightTimeOneBody = distRelativeOneBody/propSim->consts.clight;
            iter++;
        }
    }
}

void get_radar_measurement(const size_t interpIdx, const real tInterpGeom, const std::vector<real> &xInterpGeom, const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &coeffs, const std::vector<real> &tVecForInterpPrev, const std::vector< std::vector<real> > &coeffsPrev, const propSimulation *propSim, std::vector<real> &delayMeasurement){
    if (propSim->observerInfo[interpIdx].size() !=9){
        throw std::runtime_error("Error: observerInfo must be a 9-element vector for delay measurement.");
    }
    size_t numStates = xInterpGeom.size();
    real receiveTimeTDB = tInterpGeom;
    std::vector<real> xTargetBarycentricReceiveTime = xInterpGeom;
    std::vector<real> xObserverBarycentricReceiveTime(6, 0.0);
    std::vector<real> receiverInfo = {propSim->observerInfo[interpIdx][0], propSim->observerInfo[interpIdx][1], propSim->observerInfo[interpIdx][2], propSim->observerInfo[interpIdx][3]};
    get_observer_state(receiveTimeTDB, receiverInfo, propSim->consts, false, xObserverBarycentricReceiveTime);

    std::vector<real> xObserverBarycentricTransmitTime(6, 0.0);
    std::vector<real> transmitterInfo = {propSim->observerInfo[interpIdx][4], propSim->observerInfo[interpIdx][5], propSim->observerInfo[interpIdx][6], propSim->observerInfo[interpIdx][7]};

    bool bouncePointAtLeadingEdge = propSim->observerInfo[interpIdx][8] == 1.0;
    bool forwardIntegrate = tVecForInterp[0] < tVecForInterp[1];
    bool backwardIntegrate = tVecForInterp[0] > tVecForInterp[1];
    for (size_t i = 0; i < propSim->integParams.nInteg; i++){
        real delayDownleg;
        real bounceTimeTDB;
        real delayUpleg;
        real transmitTimeTDB;
        std::vector<real> xTargetBarycentricBounceTimeAllBody(numStates, 0.0);
        std::vector<real> xTargetBarycentricBounceTime(6, 0.0);
        // get downleg delay
        get_lightTimeOneBody(i, receiveTimeTDB, xTargetBarycentricReceiveTime, xObserverBarycentricReceiveTime, bouncePointAtLeadingEdge, tVecForInterp, coeffs, tVecForInterpPrev, coeffsPrev, propSim, delayDownleg);
        bounceTimeTDB = receiveTimeTDB-delayDownleg;
        if ( (forwardIntegrate && (bounceTimeTDB < tVecForInterp[0] && receiveTimeTDB != propSim->integParams.t0))
            || (backwardIntegrate && (bounceTimeTDB > tVecForInterp[0] && receiveTimeTDB != propSim->integParams.t0)) ){
            evaluate_one_interpolation(bounceTimeTDB, tVecForInterpPrev, coeffsPrev, xTargetBarycentricBounceTimeAllBody);
        } else {
            evaluate_one_interpolation(bounceTimeTDB, tVecForInterp, coeffs, xTargetBarycentricBounceTimeAllBody);
        }
        for (size_t j = 0; j < 6; j++){
            xTargetBarycentricBounceTime[j] = xTargetBarycentricBounceTimeAllBody[6*i+j];
        }
        real downlegDeltaDelayRelativistic;
        get_delta_delay_relativistic(propSim, receiveTimeTDB, xTargetBarycentricBounceTime, propSim->consts, downlegDeltaDelayRelativistic);
        // std::cout<< "downlegDeltaDelayRelativistic (receiveTime)= " << downlegDeltaDelayRelativistic*86400.0L*1e6 << " microseconds" << std::endl;
        delayDownleg += downlegDeltaDelayRelativistic;
        // iterate to get upleg delay
        delayUpleg = delayDownleg;
        if (propSim->convergedLightTime){
            real distRelativeUpleg;
            real lightTimeTol = 1e-16/86400.0L;
            real delayUplegPrev = 0.0L;
            size_t maxIter = 20;
            size_t iter = 0;
            while (iter < maxIter && fabs(delayUpleg - delayUplegPrev) > lightTimeTol){
                transmitTimeTDB = bounceTimeTDB-delayUpleg;
                get_observer_state(transmitTimeTDB, transmitterInfo, propSim->consts, false, xObserverBarycentricTransmitTime);
                vnorm({xTargetBarycentricBounceTime[0]-xObserverBarycentricTransmitTime[0], xTargetBarycentricBounceTime[1]-xObserverBarycentricTransmitTime[1], xTargetBarycentricBounceTime[2]-xObserverBarycentricTransmitTime[2]}, distRelativeUpleg);
                if (bouncePointAtLeadingEdge){
                    distRelativeUpleg -= propSim->forceParams.radii[i];
                }
                delayUplegPrev = delayUpleg;
                delayUpleg = distRelativeUpleg/propSim->consts.clight;
                iter++;
            }
        }
        transmitTimeTDB = bounceTimeTDB-delayUpleg;
        get_observer_state(transmitTimeTDB, transmitterInfo, propSim->consts, false, xObserverBarycentricTransmitTime);
        real uplegDeltaDelayRelativistic;
        get_delta_delay_relativistic(propSim, transmitTimeTDB, xTargetBarycentricBounceTime, propSim->consts, uplegDeltaDelayRelativistic);
        // std::cout<< "uplegDeltaDelayRelativistic (transmitTime)= " << uplegDeltaDelayRelativistic*86400.0L*1e6 << " microseconds" << std::endl;
        delayUpleg += uplegDeltaDelayRelativistic;
        // get delay measurement
        delayMeasurement[i] = (delayDownleg + delayUpleg)*86400.0L*1e6; // days -> seconds -> microseconds
        if (propSim->tEvalUTC){
            SpiceDouble etMinusUtcReceiveTime;
            SpiceDouble etMinusUtcTransmitTime; 
            real receiveTimeET;
            real transmitTimeET;
            mjd_to_et(receiveTimeTDB, receiveTimeET);
            mjd_to_et(transmitTimeTDB, transmitTimeET);
            deltet_c(receiveTimeET, "ET", &etMinusUtcReceiveTime);
            deltet_c(transmitTimeET, "ET", &etMinusUtcTransmitTime);
            delayMeasurement[i] += (etMinusUtcTransmitTime-etMinusUtcReceiveTime)*1e6;
        }
    }
}

void get_delta_delay_relativistic(const propSimulation *propSim, const real &tForSpice, const std::vector<real> &targetState, const Constants &consts, real &deltaDelayRelativistic){
    // from Standish (1990), https://ui.adsabs.harvard.edu/abs/1990A&A...233..252S
    SpiceDouble sunState[6];
    SpiceDouble sunLightTime;
    get_spice_state_lt(10, tForSpice, consts, sunState, sunLightTime);
    SpiceDouble earthState[6];
    SpiceDouble earthLightTime;
    get_spice_state_lt(399, tForSpice, consts, earthState, earthLightTime);

    std::vector<real> sunEarthPos = {earthState[0]-sunState[0], earthState[1]-sunState[1], earthState[2]-sunState[2]};
    real sunEarthDist;
    vnorm(sunEarthPos, sunEarthDist);
    std::vector<real> sunTargetPos = {targetState[0]-sunState[0], targetState[1]-sunState[1], targetState[2]-sunState[2]};
    real sunTargetDist;
    vnorm(sunTargetPos, sunTargetDist);
    std::vector<real> earthTargetPos = {targetState[0]-earthState[0], targetState[1]-earthState[1], targetState[2]-earthState[2]};
    real earthTargetDist;
    vnorm(earthTargetPos, earthTargetDist);

    real G = 6.6743e-11L/(149597870700.0L*149597870700.0L*149597870700.0L)*86400.0L*86400.0L; // default kg au^3 / day^2
    real sunGM = 0;
    for (size_t i=propSim->integParams.nInteg; i<propSim->integParams.nTotal; i++){
        if (propSim->forceParams.spiceIdList[i]==10){
            sunGM = G*propSim->forceParams.masses[i];
        }
    }
    if (sunGM==0){
        throw std::runtime_error("Sun GM not found in get_delta_delay_relativistic");
    }
    real c = consts.clight;
    real gamma = 1.0L; // PPN parameter

    deltaDelayRelativistic = (1+gamma)*sunGM*pow(c, -3)*log((sunEarthDist+sunTargetDist+earthTargetDist)/(sunEarthDist+sunTargetDist-earthTargetDist));
}
