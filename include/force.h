/**
 * @file    force.h
 * @brief   Header file for the force model.
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

#ifndef FORCE_H
#define FORCE_H

#include "stm.h"

/**
 * @brief Calculate the state derivative of the system.
 */
void get_state_der(PropSimulation *propSim, const real &t,
                   const std::vector<real> &xInteg,
                   std::vector<real> &accInteg);

#endif
