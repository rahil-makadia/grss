#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include "force.h"

inline void comp_sum(real num, real *sum, real *compCoeff) {
    const real y = num - *compCoeff;
    const real t = *sum + y;
    *compCoeff = (t - *sum) - y;
    *sum = t;
}
void approx_xInteg_math(const std::vector<real> &xInteg0,
                        const std::vector<real> &accInteg0, const real &dt,
                        const real &h, const std::vector<std::vector<real>> &b,
                        const size_t starti, const size_t startb,
                        const size_t &iterStep, std::vector<real> &xIntegNext,
                        std::vector<real> &xIntegCompCoeffs);
void approx_xInteg(const std::vector<real> &xInteg0,
                   const std::vector<real> &accInteg0, const real &dt,
                   const real &h, const std::vector<std::vector<real>> &b,
                   const std::vector<IntegBody> &integBodies,
                   std::vector<real> &xIntegNext,
                   std::vector<real> &xIntegCompCoeffs);
void interpolate_on_the_fly(propSimulation *propSim, const real &t, const real &dt);
void evaluate_one_interpolation(const propSimulation *propSim, const real &t,
                                const real &dt, const real &tInterp,
                                std::vector<real> &xInterp);
void get_interpIdxInWindow(const propSimulation *propSim,
                           const real &tWindowStart, const real &tNext,
                           const bool &forwardIntegrate,
                           const bool &backwardIntegrate,
                           bool &interpIdxInWindow);
void get_lightTime_and_xRelative(propSimulation *propSim,
                                 const size_t interpIdx, const real &t,
                                 const real &dt, const real tInterpGeom,
                                 const std::vector<real> &xInterpGeom,
                                 std::vector<real> &lightTime,
                                 std::vector<real> &xInterpApparent);
void get_lightTimeOneBody(propSimulation *propSim, const size_t &i,
                          const real tInterpGeom, std::vector<real> xInterpGeom,
                          std::vector<real> xObserver,
                          const bool bouncePointAtLeadingEdge, const real &t,
                          const real &dt, real &lightTimeOneBody);
void get_glb_correction(propSimulation *propSim, const real &tInterpGeom,
                        std::vector<real> &xInterpApparentBary);
void get_radar_measurement(propSimulation *propSim, const size_t &interpIdx,
                           const real &t, const real &dt,
                           const real tInterpGeom,
                           const std::vector<real> &xInterpGeom,
                           std::vector<real> &radarMeasurement,
                           std::vector<real> &radarPartials);
void get_measurement(propSimulation *propSim, const size_t &interpIdx,
                     const real &t, const real &dt, const real tInterpGeom,
                     const std::vector<real> &xInterpGeom,
                     const std::vector<real> &xInterpApparent);
void get_optical_measurement(propSimulation *propSim,
                             const std::vector<real> &xInterpApparent,
                             std::vector<real> &opticalMeasurement,
                             std::vector<real> &opticalPartials);
void get_delay_measurement(propSimulation *propSim, const size_t &interpIdx,
                           const real &t, const real &dt, const size_t &i,
                           const real tInterpGeom,
                           const std::vector<real> &xInterpGeom,
                           const real &receiveTimeTDB, real &transmitTimeTDB,
                           std::vector<real> &xObsBaryRcv,
                           std::vector<real> &xTrgtBaryBounce,
                           std::vector<real> &xObsBaryTx, real &delayMeasurement,
                           std::vector<real> &delayPartials);
void get_delta_delay_relativistic(propSimulation *propSim,
                                  const real &tForSpice,
                                  const std::vector<real> &targetState,
                                  real &deltaDelayRelativistic);
void get_doppler_measurement(propSimulation *propSim, const size_t &i,
                             const real receiveTimeTDB,
                             const real transmitTimeTDB,
                             const std::vector<real> xObsBaryRcv,
                             const std::vector<real> xTrgtBaryBounce,
                             const std::vector<real> xObsBaryTx,
                             const real transmitFreq, real &dopplerMeasurement,
                             std::vector<real> &dopplerPartials);

#endif
