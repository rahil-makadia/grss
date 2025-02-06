/**
 * @file    elements.h
 * @brief   Header file for the orbital elements.
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

#ifndef ELEMENTS_H
#define ELEMENTS_H

#include "spk.h"

/**
 * @brief Solve the Kepler equation for elliptic orbits.
 */
void kepler_solve_elliptic(const real &M, const real &e, real &E,
                           const real &tol, const int &max_iter);

/**
 * @brief Solve the Kepler equation for hyperbolic orbits.
 */
void kepler_solve_hyperbolic(const real &M, const real &e, real &EHyp,
                             const real &tol, const int &max_iter);

/**
 * @brief Solve the Kepler equation.
 */
void kepler_solve(const real &epochMjD, const std::vector<real> &cometaryState,
                  const real &GM, real &M, real &E, real &nu,
                  const real &tol = 1.0e-12L, const int &max_iter = 250);

/**
 * @brief Convert Cometary elements to Keplerian elements.
 */
void cometary_to_keplerian(const real &epochMjd,
                           const std::vector<real> &cometaryState,
                           std::vector<real> &keplerianState,
                           const real GM = 2.959122082855911e-4L);

/**
 * @brief Convert Keplerian elements to Cometary elements.
 */
void keplerian_to_cometary(const real &epochMjd,
                           const std::vector<real> &keplerianState,
                           std::vector<real> &cometaryState,
                           const real GM = 2.959122082855911e-4L);

/**
 * @brief Convert Keplerian elements to Cartesian state.
 */
void keplerian_to_cartesian(const std::vector<real> &keplerianState,
                            std::vector<real> &cartesianState,
                            const real GM = 2.959122082855911e-4L);

/**
 * @brief Convert Cartesian state to Keplerian elements.
 */
void cartesian_to_keplerian(const std::vector<real> &cartesianState,
                            std::vector<real> &keplerianState,
                            const real GM = 2.959122082855911e-4L);

/**
 * @brief Convert Cometary elements to Cartesian state.
 */
void cometary_to_cartesian(const real &epochMjd,
                           const std::vector<real> &cometaryState,
                           std::vector<real> &cartesianState,
                           const real GM = 2.959122082855911e-4L);

/**
 * @brief Convert Cartesian state to Cometary elements.
 */
void cartesian_to_cometary(const real &epochMjd,
                           const std::vector<real> &cartesianState,
                           std::vector<real> &cometaryState,
                           const real GM = 2.959122082855911e-4L);

/**
 * @brief Get the partial derivatives of the Cometary/Keplerian elements with respect to
 * the Cartesian state.
 */
void get_elements_partials(const real &epochMjd, const std::vector<real> &elems,
                           const std::string conversion,
                           std::vector<std::vector<real>> &partials,
                           const real GM = 2.959122082855911e-4L);

/**
 * @brief Get the partial derivatives of the Cartesian state with respect to the
 * Cometary/Keplerian elements.
 */
void get_cartesian_partials(const real &epochMjd,
                            const std::vector<real> &state,
                            const std::string &conversion,
                            std::vector<std::vector<real>> &partials,
                            const real GM = 2.959122082855911e-4L);

#endif
