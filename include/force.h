#ifndef FORCE_H
#define FORCE_H

#include "stm.h"

/**
 * @brief Calculate the state derivative of the system.
 */
std::vector<real> get_state_der(const real &t, const std::vector<real> &xInteg,
                                PropSimulation *propSim);

/**
 * @brief Compute the acceleration of the system due to newtonian gravity.
 */
void force_newton(const PropSimulation *propSim, std::vector<real> &accInteg,
                  std::vector<STMParameters> &allSTMs);

/**
 * @brief Compute the acceleration of the system due to the PPN relativistic correction (simple heliocentric model).
 */
void force_ppn_simple(const PropSimulation *propSim,
                      std::vector<real> &accInteg,
                      std::vector<STMParameters> &allSTMs);

/**
 * @brief Compute the acceleration of the system due to the PPN relativistic correction (full Einstein-Infeld-Hoffmann model).
 */
void force_ppn_eih(const PropSimulation *propSim, std::vector<real> &accInteg,
                   std::vector<STMParameters> &allSTMs);

/**
 * @brief Compute the acceleration of the system due to the J2 zonal harmonic.
 */
void force_J2(const PropSimulation *propSim,
              std::vector<real> &accInteg, std::vector<STMParameters> &allSTMs);

/**
 * @brief Compute the acceleration of the system due to the nongravitational forces.
 */
void force_nongrav(const PropSimulation *propSim, std::vector<real> &accInteg,
                   std::vector<STMParameters> &allSTMs);

/**
 * @brief Compute the acceleration of the system due to a thruster in the velocity direction.
 */
void force_thruster(const PropSimulation *propSim, std::vector<real> &accInteg);

#endif
