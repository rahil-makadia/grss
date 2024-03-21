#ifndef PARALLEL_H
#define PARALLEL_H

#include "gr15.h"
#include <omp.h>

std::vector<PropSimulation> propSim_parallel_omp(
    const PropSimulation refSim, const bool isCometary,
    const std::vector<std::vector<real> > &allBodies);

#endif
