/**
 * @file    interpolate.cpp
 * @brief   Source file for integrator interpolation functions.
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

#include "interpolate.h"

/**
 * @param[in] xInteg0 Initial state vector.
 * @param[in] accInteg0 Initial acceleration vector.
 * @param[in] dt Integration time step.
 * @param[in] h Fraction of the time step to use for the approximation.
 * @param[in] b Interpolation coefficients for the Gauss-Radau polynomial.
 * @param[in] dim Dimension of the system (number of 2nd derivatives).
 * @param[in] starti Starting index of the state vector to approximate.
 * @param[in] startb Starting index of the interpolation coefficients.
 * @param[in] iterStep Number of derivatives to evaluate.
 * @param[out] xIntegNext Approximated state vector.
 * @param[out] xIntegCompCoeffs Compensation coefficients.
 */
void approx_xInteg_math(const std::vector<real> &xInteg0,
                        const std::vector<real> &accInteg0, const real &dt,
                        const real &h, const std::vector<real> &b, const size_t &dim,
                        const size_t starti, const size_t startb,
                        const size_t &iterStep, std::vector<real> &xIntegNext,
                        std::vector<real> &xIntegCompCoeffs) {
    if (h == 1.0) {
        for (size_t j = 0; j < iterStep; j++) {
            xIntegNext[starti+j]          = xInteg0[starti+j];
            xIntegNext[starti+j+iterStep] = xInteg0[starti+j+iterStep];
            comp_sum(b[6*dim+startb+j]*dt*dt/72.0, &(xIntegNext[starti+j]), &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[5*dim+startb+j]*dt*dt/56.0, &(xIntegNext[starti+j]), &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[4*dim+startb+j]*dt*dt/42.0, &(xIntegNext[starti+j]), &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[3*dim+startb+j]*dt*dt/30.0, &(xIntegNext[starti+j]), &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[2*dim+startb+j]*dt*dt/20.0, &(xIntegNext[starti+j]), &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[1*dim+startb+j]*dt*dt/12.0, &(xIntegNext[starti+j]), &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[0*dim+startb+j]*dt*dt/6.0 , &(xIntegNext[starti+j]), &(xIntegCompCoeffs[starti+j]));
            comp_sum(accInteg0[startb+j]*dt*dt/2.0, &(xIntegNext[starti+j]), &(xIntegCompCoeffs[starti+j]));
            comp_sum(xInteg0[starti+j+iterStep]*dt, &(xIntegNext[starti+j]), &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[6*dim+startb+j]*dt/8.0, &(xIntegNext[starti+j+iterStep]), &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(b[5*dim+startb+j]*dt/7.0, &(xIntegNext[starti+j+iterStep]), &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(b[4*dim+startb+j]*dt/6.0, &(xIntegNext[starti+j+iterStep]), &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(b[3*dim+startb+j]*dt/5.0, &(xIntegNext[starti+j+iterStep]), &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(b[2*dim+startb+j]*dt/4.0, &(xIntegNext[starti+j+iterStep]), &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(b[1*dim+startb+j]*dt/3.0, &(xIntegNext[starti+j+iterStep]), &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(b[0*dim+startb+j]*dt/2.0, &(xIntegNext[starti+j+iterStep]), &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(accInteg0[startb+j]*dt, &(xIntegNext[starti+j+iterStep]), &(xIntegCompCoeffs[starti+j+iterStep]));
        }
    } else {
        for (size_t j = 0; j < iterStep; j++) {
            const real inc0 = ((((((((b[6*dim+startb+j]*7.*h/9. + b[5*dim+startb+j])*3.*h/4. + b[4*dim+startb+j])*5.*h/7. + b[3*dim+startb+j])*2.*h/3. + b[2*dim+startb+j])*3.*h/5. + b[1*dim+startb+j])*h/2.    + b[0*dim+startb+j])*h/3. + accInteg0[startb+j])*dt*h/2. + xInteg0[starti+j+iterStep])*dt*h;;
            const real inc1 = (((((((b[6*dim+startb+j]*7.*h/8.  + b[5*dim+startb+j])*6.*h/7. + b[4*dim+startb+j])*5.*h/6. + b[3*dim+startb+j])*4.*h/5. + b[2*dim+startb+j])*3.*h/4. + b[1*dim+startb+j])*2.*h/3. + b[0*dim+startb+j])*h/2. + accInteg0[startb+j])*dt*h;
            xIntegNext[starti+j]          = xInteg0[starti+j] + inc0 - xIntegCompCoeffs[starti+j];
            xIntegNext[starti+j+iterStep] = xInteg0[starti+j+iterStep] + inc1 - xIntegCompCoeffs[starti+j+iterStep];
        }
    }
}

/**
 * @param[in] xInteg0 Initial state vector.
 * @param[in] accInteg0 Initial acceleration vector.
 * @param[in] dt Integration time step.
 * @param[in] h Fraction of the time step to use for the approximation.
 * @param[in] b Interpolation coefficients for the Gauss-Radau polynomial.
 * @param[in] dim Dimension of the system (number of 2nd derivatives).
 * @param[in] integBodies List of integrated bodies in the PropSimulation.
 * @param[out] xIntegNext Approximated state vector.
 * @param[out] xIntegCompCoeffs Compensation coefficients.
 */
void approx_xInteg(const std::vector<real> &xInteg0,
                   const std::vector<real> &accInteg0, const real &dt,
                   const real &h, const std::vector<real> &b, const size_t &dim,
                   const std::vector<IntegBody> &integBodies,
                   std::vector<real> &xIntegNext,
                   std::vector<real> &xIntegCompCoeffs) {
    size_t starti = 0;
    size_t startb = 0;
    for (size_t i = 0; i < integBodies.size(); i++) {
        // do pos/vel first, then STM
        approx_xInteg_math(xInteg0, accInteg0, dt, h, b, dim, starti, startb, 3, xIntegNext, xIntegCompCoeffs);
        starti += 6;
        startb += 3;
        if (integBodies[i].propStm) {
            // do 6x6 STM first, then 6x1 fitted parameters
            approx_xInteg_math(xInteg0, accInteg0, dt, h, b, dim, starti, startb, 18, xIntegNext, xIntegCompCoeffs);
            starti += 36;
            startb += 18;
            if (integBodies[i].stm.size() > 36) {
                // do fitted parameters
                const size_t numParams = (integBodies[i].stm.size()-36)/6;
                for (size_t param = 0; param < numParams; param++) {
                    approx_xInteg_math(xInteg0, accInteg0, dt, h, b, dim, starti, startb, 3, xIntegNext, xIntegCompCoeffs);
                    starti += 6;
                    startb += 3;
                }
            }
        }
    }
}

/**
 * @param[in] t Time to interpolate the integrated state from the PropSimulation at.
 * @return std::vector<real> The interpolated geometric barycentric state at time t.
 */
std::vector<real> PropSimulation::interpolate(const real t) {
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
        while (idx < this->interpParams.bStack.size()-1 &&
               this->interpParams.tStack[idx+1] < t) {
            idx++;
        }
    } else if (backwardProp) {
        if (t-this->tEvalMargin > this->integParams.t0 || t+this->tEvalMargin < this->integParams.tf) {
            throw std::runtime_error("The interpolation time is outside the integration time window");
        }
        while (idx < this->interpParams.bStack.size()-1 &&
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
    const size_t dim = this->interpParams.accIntegStack[idx].size();
    std::vector<real> dummyCompCoeffs = std::vector<real>(this->xInteg.size(), 0.0);
    approx_xInteg(this->interpParams.xIntegStack[idx],
                  this->interpParams.accIntegStack[idx], dt, h,
                  this->interpParams.bStack[idx], dim, this->integBodies,
                  xIntegInterp, dummyCompCoeffs);
    return xIntegInterp;
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] t Time at the beginning of the time step.
 * @param[in] dt Completed time step size.
 */
void interpolate_on_the_fly(PropSimulation *propSim, const real &t, const real &dt) {
    const real tNext = t + dt;
    size_t &interpIdx = propSim->interpParams.interpIdx;
    const bool forwardProp = propSim->integParams.t0 < propSim->integParams.tf;
    const bool backwardProp = propSim->integParams.t0 > propSim->integParams.tf;
    bool interpIdxInWindow;
    get_interpIdxInWindow(propSim, t, tNext, forwardProp,
                          backwardProp, interpIdxInWindow);
    while (interpIdx < propSim->tEval.size() && interpIdxInWindow) {
        real tInterpGeom = propSim->tEval[interpIdx];
        if (propSim->tEvalUTC) {
            tInterpGeom += delta_et_utc(tInterpGeom)/86400.0;
        }
        std::vector<real> xInterpGeom(propSim->xInteg.size(), 0.0);
        evaluate_one_interpolation(propSim, tInterpGeom, xInterpGeom);
        if (propSim->evalApparentState) {
            std::vector<real> lightTime(propSim->integParams.nInteg, 0.0);
            std::vector<real> xInterpApparent(propSim->xInteg.size(), 0.0);
            get_lightTime_and_xRelative(propSim, interpIdx, tInterpGeom, xInterpGeom,
                                        lightTime, xInterpApparent);
            propSim->lightTimeEval.push_back(lightTime);
            propSim->xIntegEval.push_back(xInterpApparent);
            if (propSim->evalMeasurements) {
                // // THIS IS NOT NECESSARY, IT WAS ADDED FOR TESTING
                // if (propSim->obsType[interpIdx] == 3) {  // gaia
                //     // check stellar aberration necessity for gaia
                //     apply_stellar_aberration(propSim, interpIdx, xInterpApparent);
                // }
                get_measurement(propSim, interpIdx, tInterpGeom, xInterpGeom,
                                xInterpApparent);
            }
        } else {
            propSim->xIntegEval.push_back(xInterpGeom);
        }
        interpIdx++;
        get_interpIdxInWindow(propSim, t, tNext, forwardProp,
                              backwardProp, interpIdxInWindow);
    }
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] tInterp Time to interpolate to.
 * @param[out] xInterp Interpolated state vector.
 */
void evaluate_one_interpolation(const PropSimulation *propSim, const real &tInterp,
                                std::vector<real> &xInterp) {
    const real h = (tInterp - propSim->interpParams.t0) / propSim->interpParams.dt0;
    std::vector<real> dummyCompCoeffs = std::vector<real>(propSim->xInteg.size(), 0.0);
    approx_xInteg(propSim->interpParams.xInteg0,
                  propSim->interpParams.accInteg0, propSim->interpParams.dt0, h,
                  propSim->interpParams.b0, propSim->integParams.n2Derivs,
                  propSim->integBodies, xInterp, dummyCompCoeffs);
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] tWindowStart Start time of the interpolation window.
 * @param[in] tNext Time of the next interpolation.
 * @param[in] forwardProp Flag to indicate forward integration.
 * @param[in] backwardProp Flag to indicate backward integration.
 * @param[out] interpIdxInWindow Flag to indicate whether the next interpolation
 * index is within the window.
 */
void get_interpIdxInWindow(const PropSimulation *propSim,
                           const real &tWindowStart, const real &tNext,
                           const bool &forwardProp,
                           const bool &backwardProp,
                           bool &interpIdxInWindow) {
    interpIdxInWindow = false;
    const size_t interpIdx = propSim->interpParams.interpIdx;
    real tCheck = propSim->tEval[interpIdx];
    if (propSim->tEvalUTC) {
        tCheck += delta_et_utc(tCheck)/86400.0;
    }
    bool fwInWindow, fwInWindowMarginStart, fwInWindowMarginEnd;
    bool bwInWindow, bwInWindowMarginStart, bwInWindowMarginEnd;
    fwInWindow = (forwardProp && tCheck >= tWindowStart && tCheck < tNext);
    fwInWindowMarginStart =
        (forwardProp && tCheck <= propSim->integParams.t0 &&
         tCheck + propSim->tEvalMargin >= propSim->integParams.t0 &&
         tCheck + propSim->tEvalMargin >= tWindowStart);
    fwInWindowMarginEnd =
        (forwardProp && tCheck >= propSim->integParams.tf &&
         tCheck - propSim->tEvalMargin <= propSim->integParams.tf &&
         tCheck - propSim->tEvalMargin <= tNext);
    bwInWindow = (backwardProp && tCheck <= tWindowStart && tCheck > tNext);
    bwInWindowMarginStart =
        (backwardProp && tCheck >= propSim->integParams.t0 &&
         tCheck - propSim->tEvalMargin <= propSim->integParams.t0 &&
         tCheck - propSim->tEvalMargin <= tWindowStart);
    bwInWindowMarginEnd =
        (backwardProp && tCheck <= propSim->integParams.tf &&
         tCheck + propSim->tEvalMargin >= propSim->integParams.tf &&
         tCheck + propSim->tEvalMargin >= tNext);
    interpIdxInWindow = fwInWindow || fwInWindowMarginStart ||
        fwInWindowMarginEnd || bwInWindow || bwInWindowMarginStart ||
        bwInWindowMarginEnd;
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] interpIdx Index of the next interpolation time.
 * @param[in] tInterpGeom Time to interpolate to.
 * @param[in] xInterpGeom Geometric state vector of the target body.
 * @param[out] lightTime Light time to the target body.
 * @param[out] xInterpApparent Apparent state vector of the target body.
 */
void get_lightTime_and_xRelative(PropSimulation *propSim,
                                 const size_t interpIdx, const real tInterpGeom,
                                 const std::vector<real> &xInterpGeom,
                                 std::vector<real> &lightTime,
                                 std::vector<real> &xInterpApparent) {
    size_t numStates = xInterpGeom.size();
    std::vector<real> xObserver = propSim->xObserver[interpIdx];
    bool bouncePointAtCenterOfMass = true;
    if (propSim->obsType[interpIdx] == 1 ||  // delay
        propSim->obsType[interpIdx] == 2) {  // doppler
        bouncePointAtCenterOfMass = propSim->observerInfo[interpIdx][8] == 1.0;
    }
    size_t starti = 0;
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        real lightTimeTemp;
        std::vector<real> xInterpApparentBary(numStates, 0.0);
        std::vector<real> xInterpApparentOneBody(2*propSim->integBodies[i].n2Derivs, 0.0);
        get_lightTimeOneBody(propSim, i, tInterpGeom, xInterpGeom, xObserver,
                             bouncePointAtCenterOfMass, lightTimeTemp);
        evaluate_one_interpolation(propSim, tInterpGeom-lightTimeTemp, xInterpApparentBary);
        lightTime[i] = lightTimeTemp;
        for (size_t j = 0; j < 2*propSim->integBodies[i].n2Derivs; j++) {
            xInterpApparentOneBody[j] = xInterpApparentBary[starti + j];
        }
        get_glb_correction(propSim, interpIdx, tInterpGeom, xInterpApparentOneBody);
        for (size_t j = 0; j < 6; j++) {
            xInterpApparent[starti + j] = xInterpApparentOneBody[j] - xObserver[j];
        }
        for (size_t j = 6; j < 2*propSim->integBodies[i].n2Derivs; j++) {
            xInterpApparent[starti + j] = xInterpApparentOneBody[j];
        }
        starti += 2*propSim->integBodies[i].n2Derivs;
    }
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] i Index of the target body.
 * @param[in] tInterpGeom Time to interpolate to.
 * @param[in] xInterpGeom Geometric state vector of the target body.
 * @param[in] xObserver State vector of the observer.
 * @param[in] bouncePointAtCenterOfMass Flag to indicate whether the bounce point
 * is at the center of mass (as opposed to leading edge).
 * @param[out] lightTimeOneBody Light time to the target body.
 */
void get_lightTimeOneBody(PropSimulation *propSim, const size_t &i,
                          const real tInterpGeom, std::vector<real> xInterpGeom,
                          std::vector<real> xObserver,
                          const bool bouncePointAtCenterOfMass,
                          real &lightTimeOneBody) {
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
    if (!bouncePointAtCenterOfMass) {
        distRelativeOneBody -= propSim->integBodies[i].radius;
    }
    lightTimeOneBody = distRelativeOneBody / propSim->consts.clight;
    if (propSim->convergedLightTime) {
        real lightTimeTol = 1e-10 / 86400.0L;
        real lightTimeOneBodyPrev = 0.0L;
        size_t maxIter = 20;
        size_t iter = 0;
        // keep iterating until max iterations or light time tolerance is met
        while (iter < maxIter &&
               fabs(lightTimeOneBody - lightTimeOneBodyPrev) > lightTimeTol) {
            evaluate_one_interpolation(propSim, tInterpGeom - lightTimeOneBody,
                                       xInterpApparentFull);
            for (size_t j = 0; j < 6; j++) {
                xRelativeOneBody[j] =
                    xInterpApparentFull[starti + j] - xObserver[j];
                xInterpApparent[j] = xInterpApparentFull[starti + j];
            }
            vnorm(
                {xRelativeOneBody[0], xRelativeOneBody[1], xRelativeOneBody[2]},
                distRelativeOneBody);
            if (!bouncePointAtCenterOfMass) {
                distRelativeOneBody -= propSim->integBodies[i].radius;
            }
            lightTimeOneBodyPrev = lightTimeOneBody;
            lightTimeOneBody = distRelativeOneBody / propSim->consts.clight;
            iter++;
        }
        real deltaLightTimeRelativistic;
        get_delta_delay_relativistic(
            propSim, tInterpGeom, xInterpApparent,
            deltaLightTimeRelativistic);
        lightTimeOneBody += deltaLightTimeRelativistic;
        if (iter >= maxIter) {
            std::cout
                << "Warning: Downleg light time did not converge for body "
                << propSim->integBodies[i].name << " at time " << tInterpGeom
                << ", change from previous iteration was "
                << fabs(lightTimeOneBody - lightTimeOneBodyPrev) << std::endl;
        }
    }
}

/**
 * @param[in] propSim PropSimulation object for the integration.
 * @param[in] interpIdx Index of the next interpolation time.
 * @param[inout] xInterpApparent Apparent state vector of the target body.
 */
void apply_stellar_aberration(PropSimulation *propSim, const size_t interpIdx,
                              std::vector<real> &xInterpApparent) {
    // ref: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/abcorr.html#If%20correction%20for%20stellar%20aberration%20is%20requested
    const std::vector<real> xObserver = propSim->xObserver[interpIdx];
    const std::vector<real> velObserver = {xObserver[3], xObserver[4], xObserver[5]};
    real vObserver;
    vnorm(velObserver, vObserver);

    size_t starti = 0;
    for (size_t i = 0; i < propSim->integParams.nInteg; i++) {
        // // START LOGGING to compare with SPICE zzstelab.f
        // std::cout.precision(16);
        // std::cout << 
        // "VOBS(1) = " << velObserver[0] << std::endl <<
        // "VOBS(2) = " << velObserver[1] << std::endl <<
        // "VOBS(3) = " << velObserver[2] << std::endl << std::endl <<
        // "STARG(1) = " << xInterpApparent[starti] << std::endl <<
        // "STARG(2) = " << xInterpApparent[starti+1] << std::endl <<
        // "STARG(3) = " << xInterpApparent[starti+2] << std::endl <<
        // "STARG(4) = " << xInterpApparent[starti+3] << std::endl <<
        // "STARG(5) = " << xInterpApparent[starti+4] << std::endl <<
        // "STARG(6) = " << xInterpApparent[starti+5] << std::endl << std::endl;
        const std::vector<real> posApparent = {xInterpApparent[starti],
                                                xInterpApparent[starti+1],
                                                xInterpApparent[starti+2]};
        real rApparent;
        vnorm(posApparent, rApparent);
        std::vector<real> posApparentHat(3, 0.0);
        vunit(posApparent, posApparentHat);
        real velObserverDotPosApparentHat;
        vdot(velObserver, posApparentHat, velObserverDotPosApparentHat);
        std::vector<real> vp(3, 0.0);
        for (size_t j = 0; j < 3; j++) {
            vp[j] = velObserver[j] - velObserverDotPosApparentHat * posApparentHat[j];
        }
        real vpNorm;
        vnorm(vp, vpNorm);
        const real sinPhi = vpNorm / propSim->consts.clight;
        const real cosPhi = fmax(0.0, sqrt(1.0 - sinPhi*sinPhi));
        // // for testing whether i'm doing the rotation in the wrong direction (I wasn't)
        // const real phi = asin(vpNorm / propSim->consts.clight);
        // const real sinPhi = sin(phi);
        // const real cosPhi = cos(phi);
        std::vector<real> vpHat(3, 0.0);
        vunit(vp, vpHat);
        for (size_t j = 0; j < 3; j++) {
            xInterpApparent[starti+j] = rApparent * (sinPhi * vpHat[j] +
                                                cosPhi * posApparentHat[j]);
        }
        // std::cout << "VP(1) = " << vp[0] << std::endl <<
        //     "VP(2) = " << vp[1] << std::endl <<
        //     "VP(3) = " << vp[2] << std::endl <<
        //     "VPMAG = " << vpNorm << std::endl <<
        //     "PTGMAG = " << rApparent << std::endl <<
        //     "VPHAT(1) = " << vpHat[0] << std::endl <<
        //     "VPHAT(2) = " << vpHat[1] << std::endl <<
        //     "VPHAT(3) = " << vpHat[2] << std::endl <<
        //     "RHAT(1) = " << posApparentHat[0] << std::endl <<
        //     "RHAT(2) = " << posApparentHat[1] << std::endl <<
        //     "RHAT(3) = " << posApparentHat[2] << std::endl <<
        //     "C = " << cosPhi << std::endl <<
        //     "S = " << sinPhi << std::endl <<
        //     "STARG(1) = " << xInterpApparent[starti] << std::endl <<
        //     "STARG(2) = " << xInterpApparent[starti+1] << std::endl <<
        //     "STARG(3) = " << xInterpApparent[starti+2] << std::endl <<
        //     "STARG(4) = " << xInterpApparent[starti+3] << std::endl <<
        //     "STARG(5) = " << xInterpApparent[starti+4] << std::endl <<
        //     "STARG(6) = " << xInterpApparent[starti+5] << std::endl << std::endl;
        // throw std::runtime_error(
        //     "apply_stellar_aberration: Not implemented yet");
        starti += 2*propSim->integBodies[i].n2Derivs;
    }
}
