/**
 * @file    grss.h
 * @brief   Top level header file for GRSS
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

#ifndef GRSS_H
#define GRSS_H

#include "parallel.h"

void PropSimulation::integrate() {
    // raise error if xInteg is empty (no integration bodies in simulation)
    if (this->integParams.nInteg < 1) {
        throw std::runtime_error(
            "ERROR: There are no integration bodies in the "
            "simulation. Need at least one body with a full state vector to "
            "integrate.");
    }

    // integrate the system
    this->map_ephemeris();
    this->preprocess();
    ias15(this);
    if (!this->unsafePersistentMemoryMap) {
        this->unmap_ephemeris();
    }

    // flip vectors if integration is backwards in time
    if (this->integParams.t0 > this->integParams.tf) {
        std::reverse(this->eventMngr.impulsiveEvents.begin(), this->eventMngr.impulsiveEvents.end());
        std::reverse(this->eventMngr.continuousEvents.begin(), this->eventMngr.continuousEvents.end());
        std::reverse(this->xObserver.begin(), this->xObserver.end());
        std::reverse(this->observerInfo.begin(), this->observerInfo.end());
        std::reverse(this->tEval.begin(), this->tEval.end());
        std::reverse(this->obsType.begin(), this->obsType.end());
        std::reverse(this->lightTimeEval.begin(), this->lightTimeEval.end());
        std::reverse(this->xIntegEval.begin(), this->xIntegEval.end());
        std::reverse(this->opticalObs.begin(), this->opticalObs.end());
        std::reverse(this->opticalObsDot.begin(), this->opticalObsDot.end());
        std::reverse(this->opticalPartials.begin(), this->opticalPartials.end());
        std::reverse(this->opticalObsCorr.begin(), this->opticalObsCorr.end());
        std::reverse(this->radarObs.begin(), this->radarObs.end());
        std::reverse(this->radarPartials.begin(), this->radarPartials.end());
    }
}

#endif
