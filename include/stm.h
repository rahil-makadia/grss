#ifndef STM_H
#define STM_H

#include "simulation.h"

void bcd_and_dot(const std::vector<real> &stm, real *B, real *Bdot, real *C,
                 real *Cdot, real *D, real *Ddot);
void bcd_2dot(STMParameters &stmParams, size_t numParams,
              size_t stmStarti, std::vector<real> &accInteg);
void stm_newton(STMParameters &stmParams, const real &gm, const real &dx,
                const real &dy, const real &dz);
void stm_ppn_simple(STMParameters &stmParams, const real &gm, const real &c,
                    const real &beta, const real &gamma, const real &dx,
                    const real &dy, const real &dz, const real &dvx,
                    const real &dvy, const real &dvz);
void stm_J2(STMParameters &stmParams, const real &gm, const real &J2,
            const real &dxBody, const real &dyBody, const real &dzBody,
            const real &radius, const real &sinRA, const real &cosRA,
            const real &sinDec, const real &cosDec,
            const real &smoothing_threshold);
void stm_nongrav(STMParameters &stmParams, const real &g,
                 const NongravParamaters &ngParams, const real &dx,
                 const real &dy, const real &dz, const real &dvx,
                 const real &dvy, const real &dvz, real *rVec, real *nVec);

#endif
