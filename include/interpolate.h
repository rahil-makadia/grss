/**
 * @file    interpolate.h
 * @brief   Header file for integrator interpolation functions.
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

#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include "observe.h"

/**
 * @brief Compute the sum of two real numbers with a compensated summation.
 * 
 * @param[in] num Number to add to the sum.
 * @param[inout] sum Sum of the numbers.
 * @param[inout] compCoeff Compensation coefficient.
 */
static inline void comp_sum(real num, real *sum, real *compCoeff) {
    const real y = num - *compCoeff;
    const real t = *sum + y;
    *compCoeff = (t - *sum) - y;
    *sum = t;
}

/**
 * @brief Evaluate the Gauss-Radau polynomial.
 */
void approx_xInteg_math(const std::vector<real> &xInteg0,
                        const std::vector<real> &accInteg0, const real &dt,
                        const real &h, const std::vector<real> &b, const size_t &dim,
                        const size_t starti, const size_t startb,
                        const size_t &iterStep, std::vector<real> &xIntegNext,
                        std::vector<real> &xIntegCompCoeffs);

/**
 * @brief Use the Gauss-Radau polynomial to approximate the integral of the
 * acceleration.
 */
void approx_xInteg(const std::vector<real> &xInteg0,
                   const std::vector<real> &accInteg0, const real &dt,
                   const real &h, const std::vector<real> &b, const size_t &dim,
                   const std::vector<IntegBody> &integBodies,
                   std::vector<real> &xIntegNext,
                   std::vector<real> &xIntegCompCoeffs);

/**
 * @brief Interpolate the integrator state for the time step that was just
 * completed.
 */
void interpolate_on_the_fly(PropSimulation *propSim, const real &t, const real &dt);

/**
 * @brief Determine whether the next interpolation index is within the window
 * of the time step that was just completed.
 */
void get_interpIdxInWindow(const PropSimulation *propSim,
                           const real &tWindowStart, const real &tNext,
                           const bool &forwardProp,
                           const bool &backwardProp,
                           bool &interpIdxInWindow);

/**
 * @brief Compute the light time and apparent position of the target body.
 */
void get_lightTime_and_xRelative(PropSimulation *propSim,
                                 const size_t interpIdx, const real tInterpGeom,
                                 const std::vector<real> &xInterpGeom,
                                 std::vector<real> &lightTime,
                                 std::vector<real> &xInterpApparent);

/**
 * @brief Compute the light time to the target body.
 */
void get_lightTimeOneBody(PropSimulation *propSim, const size_t &i,
                          const real tInterpGeom, std::vector<real> xInterpGeom,
                          std::vector<real> xObserver,
                          const bool bouncePointAtCenterOfMass, real &lightTimeOneBody);

/**
 * @brief Apply stellar aberration correction to the apparent position of the target
 * body.
 */
void apply_stellar_aberration(PropSimulation *propSim, const size_t interpIdx,
                             std::vector<real> &xInterpApparent);

#endif
