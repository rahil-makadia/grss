#ifndef GRSS_H
#define GRSS_H

#include "gr15.h"

void Simulation::integrate(){
    // check that xInteg has a valid size. If not, raise error
    if (this->integParams.nInteg < 1){
        throw std::runtime_error("\n\ngrss.h: ERROR: There are no integration bodies in the simulation. Need at least one body with a full state vector to integrate.\n");
    }

    gr15(this->t, this->xInteg, *this);
}
#endif
