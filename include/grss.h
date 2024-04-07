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
    gr15(this);
    this->unmap_ephemeris();

    // flip vectors if integration is backwards in time
    if (this->integParams.t0 > this->integParams.tf) {
        std::reverse(this->events.begin(), this->events.end());
        std::reverse(this->xObserver.begin(), this->xObserver.end());
        std::reverse(this->observerInfo.begin(), this->observerInfo.end());
        std::reverse(this->tEval.begin(), this->tEval.end());
        std::reverse(this->radarObserver.begin(), this->radarObserver.end());
        std::reverse(this->lightTimeEval.begin(), this->lightTimeEval.end());
        std::reverse(this->xIntegEval.begin(), this->xIntegEval.end());
        std::reverse(this->opticalObs.begin(), this->opticalObs.end());
        std::reverse(this->opticalPartials.begin(), this->opticalPartials.end());
        std::reverse(this->radarObs.begin(), this->radarObs.end());
        std::reverse(this->radarPartials.begin(), this->radarPartials.end());
    }
}

#endif
