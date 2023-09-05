#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include "force.h"
#include "gr15.h"

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
void get_radar_measurement(propSimulation *propSim, const size_t interpIdx,
                           const real &t, const real &dt,
                           const real tInterpGeom,
                           const std::vector<real> &xInterpGeom,
                           std::vector<real> &radarMeasurement);
void get_delta_delay_relativistic(propSimulation *propSim,
                                  const real &tForSpice,
                                  const std::vector<real> &targetState,
                                  real &deltaDelayRelativistic);
void get_doppler_measurement(propSimulation *propSim, const real receiveTimeTDB,
                             const real transmitTimeTDB,
                             const std::vector<real> xObsBaryRcv,
                             const std::vector<real> xTrgtBaryBounce,
                             const std::vector<real> xObsBaryTx,
                             const real transmitFreq, real &dopplerMeasurement);

#endif
