/**
 * @file    parallel.h
 * @brief   Header file for parallel propagation of bodies using OpenMP.
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

#ifndef PARALLEL_H
#define PARALLEL_H

#include "ias15.h"
#include <omp.h>

/**
 * @brief Propagate an array of bodies in parallel via OpenMP, using a reference
 * simulation as a template.
 */
void propSim_parallel_omp(const PropSimulation refSim, const bool isCometary,
                          const std::vector<std::vector<real> > &allBodies,
                          const int &maxThreads = 128);

#endif
