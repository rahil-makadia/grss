#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include "force.h"

/**
 * @brief Compute the sum of two real numbers with a compensated summation.
 * 
 * @param[in] num Number to add to the sum.
 * @param[inout] sum Sum of the numbers.
 * @param[inout] compCoeff Compensation coefficient.
 */
static inline void comp_sum(real num, real *sum, real *compCoeff) {
    const real y = num - *compCoeff;
    const real t = *sum + y;
    *compCoeff = (t - *sum) - y;
    *sum = t;
}

/**
 * @brief Evaluate the Gauss-Radau polynomial.
 */
void approx_xInteg_math(const std::vector<real> &xInteg0,
                        const std::vector<real> &accInteg0, const real &dt,
                        const real &h, const std::vector<std::vector<real>> &b,
                        const size_t starti, const size_t startb,
                        const size_t &iterStep, std::vector<real> &xIntegNext,
                        std::vector<real> &xIntegCompCoeffs);

/**
 * @brief Use the Gauss-Radau polynomial to approximate the integral of the
 * acceleration.
 */
void approx_xInteg(const std::vector<real> &xInteg0,
                   const std::vector<real> &accInteg0, const real &dt,
                   const real &h, const std::vector<std::vector<real>> &b,
                   const std::vector<IntegBody> &integBodies,
                   std::vector<real> &xIntegNext,
                   std::vector<real> &xIntegCompCoeffs);

/**
 * @brief Interpolate the integrator state for the time step that was just
 * completed.
 */
void interpolate_on_the_fly(PropSimulation *propSim, const real &t, const real &dt);

/**
 * @brief Interpolate the integrator state for one evaluation time.
 */
void evaluate_one_interpolation(const PropSimulation *propSim, const real &t,
                                const real &dt, const real &tInterp,
                                std::vector<real> &xInterp);

/**
 * @brief Determine whether the next interpolation index is within the window
 * of the time step that was just completed.
 */
void get_interpIdxInWindow(const PropSimulation *propSim,
                           const real &tWindowStart, const real &tNext,
                           const bool &forwardProp,
                           const bool &backwardProp,
                           bool &interpIdxInWindow);

/**
 * @brief Compute the light time and apparent position of the target body.
 */
void get_lightTime_and_xRelative(PropSimulation *propSim,
                                 const size_t interpIdx, const real &t,
                                 const real &dt, const real tInterpGeom,
                                 const std::vector<real> &xInterpGeom,
                                 std::vector<real> &lightTime,
                                 std::vector<real> &xInterpApparent);

/**
 * @brief Compute the light time to the target body.
 */
void get_lightTimeOneBody(PropSimulation *propSim, const size_t &i,
                          const real tInterpGeom, std::vector<real> xInterpGeom,
                          std::vector<real> xObserver,
                          const bool bouncePointAtLeadingEdge, const real &t,
                          const real &dt, real &lightTimeOneBody);

/**
 * @brief Compute the correction to the apparent state of the body due to the
 * gravitational light bending.
 */
void get_glb_correction(PropSimulation *propSim, const real &tInterpGeom,
                        std::vector<real> &xInterpApparentBary);

/**
 * @brief Get the relevant measurement (optical/radar) for a given measurement time.
 */
void get_measurement(PropSimulation *propSim, const size_t &interpIdx,
                     const real &t, const real &dt, const real tInterpGeom,
                     const std::vector<real> &xInterpGeom,
                     const std::vector<real> &xInterpApparent);

/**
 * @brief Get the optical measurement and partials.
 */
void get_optical_measurement(PropSimulation *propSim,
                             const std::vector<real> &xInterpApparent,
                             std::vector<real> &opticalMeasurement,
                             std::vector<real> &opticalPartials);

/**
 * @brief Get the radar measurement and partials.
 */
void get_radar_measurement(PropSimulation *propSim, const size_t &interpIdx,
                           const real &t, const real &dt,
                           const real tInterpGeom,
                           const std::vector<real> &xInterpGeom,
                           std::vector<real> &radarMeasurement,
                           std::vector<real> &radarPartials);

/**
 * @brief Get the radar delay measurement and partials.
 */
void get_delay_measurement(PropSimulation *propSim, const size_t &interpIdx,
                           const real &t, const real &dt, const size_t &i,
                           const real tInterpGeom,
                           const std::vector<real> &xInterpGeom,
                           const real &receiveTimeTDB, real &transmitTimeTDB,
                           std::vector<real> &xObsBaryRcv,
                           std::vector<real> &xTrgtBaryBounce,
                           std::vector<real> &xObsBaryTx, real &delayMeasurement,
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

#endif
