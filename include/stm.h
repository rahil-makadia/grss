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
void stm_nongrav(const IntegBody &bodyi, const real &g,
                 const NongravParamaters &ngParams, const real &dx,
                 const real &dy, const real &dz, const real &dvx,
                 const real &dvy, const real &dvz, real *rVec, real *nVec,
                 const size_t &stmStarti, std::vector<real> &accInteg);

#endif
