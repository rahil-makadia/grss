/**
 * @file    approach.h
 * @brief   Header file for close approach and impact detection.
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

#ifndef CLOSEAPPROACH_H
#define CLOSEAPPROACH_H

#include "interpolate.h"

/**
 * @brief Check for a close approach or impact between two bodies.
 */
void check_ca_or_impact(PropSimulation *propSim, const real &tOld,
                        const std::vector<real> xIntegOld, const real &t,
                        const std::vector<real> xInteg);

/**
 * @brief Compute relative radial velocity to check for a close approach.
 */
void ca_rdot_calc(PropSimulation *propSim, const size_t &i, const size_t &j,
                  const real &t, real &rDot);

/**
 * @brief Compute relative distance to check for an impact.
 */
void impact_r_calc(PropSimulation *propSim, const size_t &i, const size_t &j,
                  const real &t, real &r);

/**
 * @brief Compute relative state of a body at a given time.
 */
std::vector<real> get_rel_state(PropSimulation *propSim, const size_t &i,
                                const size_t &j, const real &t);

/**
 * @brief Compute partials of B-plane parameters.
 */
void get_bplane_partials(PropSimulation *propSim, CloseApproachParameters *ca,
                         const real &mu, const real &radius);

/**
 * @brief Compute the time of close approach or impact using Brent's method
 *       for root finding in a bracketed interval.
 */
void get_ca_or_impact_time(PropSimulation *propSim, const size_t &i,
                           const size_t &j, const real &x1, const real &x2,
                           real &tCA,
                           void (*zero_func)(PropSimulation *, const size_t &,
                                             const size_t &, const real &,
                                             real &));

#endif
