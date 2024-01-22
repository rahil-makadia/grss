#ifndef CLOSEAPPROACH_H
#define CLOSEAPPROACH_H

#include "interpolate.h"

void check_ca_or_impact(propSimulation *propSim, const real &tOld,
                        const std::vector<real> xIntegOld, const real &t,
                        const std::vector<real> xInteg);
void ca_rdot_calc(propSimulation *propSim, const size_t &i, const size_t &j,
                  const real &t, real &rDot);
void impact_r_calc(propSimulation *propSim, const size_t &i, const size_t &j,
                  const real &t, real &r);
void get_rel_state(propSimulation *propSim, const size_t &i, const size_t &j,
                   const real &t, real xRel[6]);
void get_ca_or_impact_time(propSimulation *propSim, const size_t &i,
                           const size_t &j, const real &x1, const real &x2,
                           real &tCA,
                           void (*zero_func)(propSimulation *, const size_t &,
                                             const size_t &, const real &,
                                             real &));
#endif
