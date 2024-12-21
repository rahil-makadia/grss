/**
 * @file    stm.h
 * @brief   Header file for the STM functions.
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

#ifndef STM_H
#define STM_H

#include "simulation.h"

/**
 * @brief Structure to hold the STM submatrices and their derivatives.
 * 
 * @param B The B submatrix (dpos_dpos).
 * @param Bdot The Bdot submatrix (dpos_dvel).
 * @param C The C submatrix (dvel_dpos).
 * @param Cdot The Cdot submatrix (dvel_dvel).
 * @param D The D submatrix (dpos_dpar).
 * @param Ddot The Ddot submatrix (dvel_dpar).
 * @param dfdpos The partial derivative of the dynamics with respect to position.
 * @param dfdvel The partial derivative of the dynamics with respect to velocity.
 * @param dfdpar The partial derivative of the dynamics with respect to parameters.
 */
struct STMParameters {
    real *B = nullptr;
    real *Bdot = nullptr;
    real *C = nullptr;
    real *Cdot = nullptr;
    real *D = nullptr;
    real *Ddot = nullptr;
    real *dfdpos = nullptr;
    real *dfdvel = nullptr;
    real *dfdpar = nullptr;
};

/**
 * @brief Unpack the STM submatrices.
 */
void bcd_and_dot(const std::vector<real> &stm, real *B, real *Bdot, real *C,
                 real *Cdot, real *D, real *Ddot);

/**
 * @brief Compute the derivatives of the STM submatrices.
 */
void bcd_2dot(STMParameters &stmParams, size_t numParams,
              size_t stmStarti, std::vector<real> &accInteg);

/**
 * @brief Compute the derivatives of newtonian gravity with respect to position.
 */
void stm_newton(STMParameters &stmParams, const real &gm, const real &dx,
                const real &dy, const real &dz);

/**
 * @brief Compute the derivatives of the PPN relativistic correction with respect to position and velocity.
 */
void stm_ppn_simple(STMParameters &stmParams, const real &gm, const real &c,
                    const real &beta, const real &gamma, const real &dx,
                    const real &dy, const real &dz, const real &dvx,
                    const real &dvy, const real &dvz);

/**
 * @brief Compute the derivatives of the J2 zonal harmonic with respect to position.
 */
void stm_J2(STMParameters &stmParams, const real &gm, const real &J2,
            const real &dxBody, const real &dyBody, const real &dzBody,
            const real &radius, const real &sinRA, const real &cosRA,
            const real &sinDec, const real &cosDec,
            const real &smoothing_threshold);

/**
 * @brief Compute the derivatives of the non-gravitational acceleration with respect to position, velocity, and parameters.
 */
void stm_nongrav(STMParameters &stmParams, const real &g,
                 const NongravParameters &ngParams, const real &dx,
                 const real &dy, const real &dz, const real &dvx,
                 const real &dvy, const real &dvz, real *rVec, real *nVec);

/**
 * @brief Compute the derivatives of a continuous event with respect to position, velocity, and parameters.
 */
void stm_continuous_event(STMParameters &stmParams,
                          const PropSimulation *propSim, const size_t &eventIdx,
                          const real &tPastEvent, const real &postFac);
#endif
