#ifndef STM_H
#define STM_H

#include "simulation.h"

void stm_newton(const IntegBody &bodyi, const real &gm, const real &dx,
                const real &dy, const real &dz, const size_t &stmStarti,
                std::vector<real> &accInteg);
void stm_ppn_simple(const IntegBody &bodyi, const real &gm, const real &c,
                    const real &beta, const real &gamma, const real &dx,
                    const real &dy, const real &dz, const real &dvx,
                    const real &dvy, const real &dvz, const size_t &stmStarti,
                    std::vector<real> &accInteg);
void stm_J2(const IntegBody &bodyi, const real &gm, const real &J2,
            const real &dxBody, const real &dyBody, const real &dzBody,
            const real &radius, const real &sinRA, const real &cosRA,
            const real &sinDec, const real &cosDec,
            const real &smoothing_threshold, const size_t &stmStarti,
            std::vector<real> &accInteg);
void stm_nongrav(const IntegBody &bodyi, const real &g,
                 const NongravParamaters &ngParams, const real &dx,
                 const real &dy, const real &dz, const real &dvx,
                 const real &dvy, const real &dvz, real *rVec, real *nVec,
                 const size_t &stmStarti, std::vector<real> &accInteg);

#endif
