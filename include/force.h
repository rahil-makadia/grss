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
