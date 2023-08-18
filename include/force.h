#ifndef FORCE_H
#define FORCE_H
#include "simulation.h"

std::vector<real> get_state_der(const real &t, const std::vector<real> &xInteg,
                                propSimulation *propSim);
void force_newton(const propSimulation *propSim, real *accInteg);
void force_ppn_simple(const propSimulation *propSim, real *accInteg);
void force_ppn_eih(const propSimulation *propSim, real *accInteg);
void force_J2(const propSimulation *propSim, real *accInteg);
void force_nongrav(const propSimulation *propSim, real *accInteg);
void force_thruster(const propSimulation *propSim, real *accInteg);

#endif
