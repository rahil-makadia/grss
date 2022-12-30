#ifndef GRSS_H
#define GRSS_H

#include "gr15.h"

void Simulation::integrate(){
    gr15(this->t, this->xInteg, *this);
}
#endif
