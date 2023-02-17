#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include "simulation.h"
#include "gr15.h"

void interpolate(const real &t, const real &dt, const std::vector<real> &xInteg0, const std::vector<real> &accInteg0, const std::vector< std::vector<real> > &b, const std::vector<real> &hVec, Simulation &sim);
void get_coeffs(const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &xIntegForInterp, std::vector< std::vector<real> > &coeffs);
void evaluate_one_interpolation(const real &tInterp, const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &coeffs, std::vector<real> &xInterp);
void one_timestep_interpolation(const real &tNext, const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &xIntegForInterp, const std::vector<real> &tVecForInterpPrev, const std::vector< std::vector<real> > &xIntegForInterpPrev, Simulation &sim);

#endif