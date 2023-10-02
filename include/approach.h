#ifndef CLOSEAPPROACH_H
#define CLOSEAPPROACH_H

#include "interpolate.h"

void check_ca_or_impact(propSimulation *propSim, const real &tOld,
                        const std::vector<real> xIntegOld, const real &t,
                        const std::vector<real> xInteg, int &keepStepping);
void ca_rdot_calc(propSimulation *propSim, const size_t &i, const size_t &j,
                  const real &t, real &rDot);
void get_ca_state(propSimulation *propSim, const size_t &i, const size_t &j,
                  const real &t, real xRelCA[6]);
void get_ca_time(propSimulation *propSim, const size_t &i, const size_t &j,
                 const real &x1, const real &x2, real &tCA);
#endif
