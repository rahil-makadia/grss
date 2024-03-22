#ifndef FORCE_H
#define FORCE_H

#include "stm.h"

std::vector<real> get_state_der(const real &t, const std::vector<real> &xInteg,
                                PropSimulation *propSim);
void force_newton(const PropSimulation *propSim, std::vector<real> &accInteg,
                  std::vector<STMParameters> &allSTMs);
void force_ppn_simple(const PropSimulation *propSim,
                      std::vector<real> &accInteg,
                      std::vector<STMParameters> &allSTMs);
void force_ppn_eih(const PropSimulation *propSim, std::vector<real> &accInteg,
                   std::vector<STMParameters> &allSTMs);
void force_J2(const PropSimulation *propSim, std::vector<real> &accInteg,
              std::vector<STMParameters> &allSTMs);
void force_nongrav(const PropSimulation *propSim, std::vector<real> &accInteg,
                   std::vector<STMParameters> &allSTMs);
void force_thruster(const PropSimulation *propSim, std::vector<real> &accInteg);

#endif
