#ifndef CLOSEAPPROACH_H
#define CLOSEAPPROACH_H

#include "interpolate.h"

/**
 * @brief Check for a close approach or impact between two bodies.
 */
void check_ca_or_impact(PropSimulation *propSim, const real &tOld,
                        const real &t);

/**
 * @brief Compute relative radial velocity to check for a close approach.
 */
void ca_rdot_calc(PropSimulation *propSim, const size_t &i, const size_t &j,
                  const real &t, real &r, real &rDot);

/**
 * @brief Compute relative distance to check for an impact.
 */
void impact_r_calc(PropSimulation *propSim, const size_t &i, const size_t &j,
                  const real &t, real &r, real &dummy);

/**
 * @brief Compute relative state of a body at a given time.
 */
void get_rel_state(PropSimulation *propSim, const size_t &i, const size_t &j,
                   const real &t, real xRel[6]);

/**
 * @brief Compute the time of close approach or impact using Brent's method
 *       for root finding in a bracketed interval.
 */
void get_ca_or_impact_time(PropSimulation *propSim, const size_t &i,
                           const size_t &j, const real &x1, const real &x2,
                           real &tCA,
                           void (*zero_func)(PropSimulation *, const size_t &,
                                             const size_t &, const real &,
                                             real &, real &));

#endif
