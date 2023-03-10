#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include "simulation.h"
#include "gr15.h"

void interpolate(const real &t, const real &dt, const std::vector<real> &xInteg0, const std::vector<real> &accInteg0, const std::vector< std::vector<real> > &b, const std::vector<real> &hVec, propSimulation &propSim);
void get_coeffs(const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &xIntegForInterp, std::vector< std::vector<real> > &coeffs);
void evaluate_one_interpolation(const real &tInterp, const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &coeffs, std::vector<real> &xInterp);
void one_timestep_interpolation(const real &tNext, const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &xIntegForInterp, const std::vector<real> &tVecForInterpPrev, const std::vector< std::vector<real> > &xIntegForInterpPrev, propSimulation &propSim);
void get_lightTime_and_xRelative(const size_t interpIdx, const real tInterpGeom, const std::vector<real> &xInterpGeom, const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &coeffs, const std::vector<real> &tVecForInterpPrev, const std::vector< std::vector<real> > &coeffsPrev, const propSimulation &propSim, std::vector<real> &lightTime, std::vector<real> &xInterpApparent);
void get_lightTimeOneBody(const size_t &i, const real tInterpGeom, std::vector<real> xInterpGeom, std::vector<real> xObserver, const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &coeffs, const std::vector<real> &tVecForInterpPrev, const std::vector< std::vector<real> > &coeffsPrev, const propSimulation &propSim, real &lightTimeOneBody);
void get_radar_measurement(const size_t interpIdx, const real tInterpGeom, const std::vector<real> &xInterpGeom, const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &coeffs, const std::vector<real> &tVecForInterpPrev, const std::vector< std::vector<real> > &coeffsPrev, const propSimulation &propSim, std::vector<real> &delayMeasurement);
void get_delta_delay_relativistic(const real &t, const std::vector<real> &targetState, const Constants &consts, real &deltaDelayRelativistic);

#endif
