#ifndef GRSS_H
#define GRSS_H

#include "parallel.h"

void propSimulation::integrate() {
    // check that xInteg has a valid size. If not, raise error
    if (this->integParams.nInteg < 1) {
        throw std::runtime_error(
            "ERROR: There are no integration bodies in the "
            "simulation. Need at least one body with a full state vector to "
            "integrate.");
    }

    this->preprocess();
    furnsh_c(this->DEkernelPath.c_str());
    gr15(this);
    unload_c(this->DEkernelPath.c_str());

    // if backwards integration
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
