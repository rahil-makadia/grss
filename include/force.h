#ifndef FORCE_H
#define FORCE_H

#include "stm.h"

/**
 * @brief Calculate the state derivative of the system.
 */
std::vector<real> get_state_der(const real &t, const std::vector<real> &xInteg,
                                PropSimulation *propSim);

#endif
