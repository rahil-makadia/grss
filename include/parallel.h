#ifndef PARALLEL_H
#define PARALLEL_H

#include "gr15.h"
#include <omp.h>

std::vector<propSimulation> propSim_parallel_omp(
    const propSimulation refSim, const bool isCometary,
    const std::vector<std::vector<real> > &allBodies);

#endif
