#ifndef OBSERVE_H
#define OBSERVE_H

#include "force.h"

/**
 * @brief Compute the correction to the apparent state of the body due to the
 * gravitational light bending.
 */
void get_glb_correction(PropSimulation *propSim, const size_t &interpIdx,
                        const real &tInterpGeom,
                        std::vector<real> &xInterpApparentBary);

/**
 * @brief Get the relevant measurement (optical/radar) for a given measurement time.
 */
void get_measurement(PropSimulation *propSim, const size_t &interpIdx,
                     const real tInterpGeom,
                     const std::vector<real> &xInterpGeom,
                     const std::vector<real> &xInterpApparent);

/**
 * @brief Get the optical measurement and partials.
 */
void get_optical_measurement(PropSimulation *propSim,
                             const std::vector<real> &xInterpApparent,
                             std::vector<real> &opticalMeasurement,
                             std::vector<real> &opticalMeasurementDot,
                             std::vector<real> &opticalPartials);

/**
 * @brief Get the photocenter-barycenter correction for an optical measurement.
 */
void get_photocenter_correction(PropSimulation *propSim, const size_t &interpIdx,
                                const real &tInterpGeom,
                                const std::vector<real> &xInterpApparent,
                                std::vector<real> &photocenterCorr);

/**
 * @brief Get the radar measurement and partials.
 */
void get_radar_measurement(PropSimulation *propSim, const size_t &interpIdx,
                           const real tInterpGeom,
                           const std::vector<real> &xInterpGeom,
                           std::vector<real> &radarMeasurement,
                           std::vector<real> &radarPartials);

/**
 * @brief Get the radar delay measurement and partials.
 */
void get_delay_measurement(PropSimulation *propSim, const size_t &interpIdx,
                           const size_t &i, const real tInterpGeom,
                           const std::vector<real> &xInterpGeom,
                           const real &receiveTimeTDB, real &transmitTimeTDB,
                           std::vector<real> &xObsBaryRcv,
                           std::vector<real> &xTrgtBaryBounce,
                           std::vector<real> &xObsBaryTx,
                           real &delayMeasurement,
                           std::vector<real> &delayPartials);

/**
 * @brief Get the relativistic delay measurement correction.
 */
void get_delta_delay_relativistic(PropSimulation *propSim,
                                  const real &tForSpice,
                                  const std::vector<real> &targetState,
                                  real &deltaDelayRelativistic);

/**
 * @brief Get the Doppler measurement and partials.
 */
void get_doppler_measurement(PropSimulation *propSim, const size_t &i,
                             const real receiveTimeTDB,
                             const real transmitTimeTDB,
                             const std::vector<real> xObsBaryRcv,
                             const std::vector<real> xTrgtBaryBounce,
                             const std::vector<real> xObsBaryTx,
                             const real transmitFreq, real &dopplerMeasurement,
                             std::vector<real> &dopplerPartials);

/**
 * @brief Interpolate the integrator state for one evaluation time.
 */
void evaluate_one_interpolation(
    const PropSimulation *propSim, const real &tInterp,
    std::vector<real> &xInterp);  // defined in interpolate.cpp

#endif
