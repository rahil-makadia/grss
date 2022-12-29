#ifndef INTERPOLATE_H
#define INTERPOLATE_H
#include "utilities.h"
#include "simulation.h"

void interpolate(const real &tNext, const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &xIntegForInterp, Simulation &sim);

#endif