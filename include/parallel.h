#ifndef PARALLEL_H
#define PARALLEL_H

#include "gr15.h"
#include <omp.h>

/**
 * @brief Propagate an array of bodies in parallel via OpenMP, using a reference
 * simulation as a template.
 */
void propSim_parallel_omp(const PropSimulation refSim, const bool isCometary,
                          const std::vector<std::vector<real> > &allBodies);

#endif
