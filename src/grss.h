#include "utilities.h"
#include "force.h"
#include "simulation.h"
#include "gr15.h"

void Simulation::integrate(){
    gr15(this->t, this->xInteg, *this);
}
