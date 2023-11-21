#include "gr15.h"

real get_initial_timestep(propSimulation *propSim){
    real dt0 = propSim->integParams.dtMin;
    if (propSim->integParams.dt0 != 0.0) {
        dt0 = fabs(propSim->integParams.dt0);
    }
    real span = fabs(propSim->integParams.tf - propSim->integParams.t0);
    if (span < dt0) {
        dt0 = span;
    }
    if (propSim->integParams.tf < propSim->integParams.t0) {
        dt0 *= -1.0;
    }
    return dt0;
}

static inline void comp_sum(real num, real *sum, real *compCoeff) {
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
                        std::vector<real> &xIntegCompCoeffs) {
    for (size_t j = 0; j < iterStep; j++) {
        real inc0, inc1;
        if (h == 1.0) {
            inc0 = 0;
            real dummy = 0;
            comp_sum(b[6][startb+j]*dt*dt/72.0, &inc0, &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[5][startb+j]*dt*dt/56.0, &inc0, &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[4][startb+j]*dt*dt/42.0, &inc0, &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[3][startb+j]*dt*dt/30.0, &inc0, &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[2][startb+j]*dt*dt/20.0, &inc0, &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[1][startb+j]*dt*dt/12.0, &inc0, &(xIntegCompCoeffs[starti+j]));
            comp_sum(b[0][startb+j]*dt*dt/6.0 , &inc0, &(xIntegCompCoeffs[starti+j]));
            comp_sum(accInteg0[startb+j]*dt*dt/2.0, &inc0, &(xIntegCompCoeffs[starti+j]));
            comp_sum(xInteg0[starti+j+iterStep]*dt, &inc0, &(xIntegCompCoeffs[starti+j]));
            inc1 = 0;
            comp_sum(b[6][startb+j]*dt/8.0, &inc1, &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(b[5][startb+j]*dt/7.0, &inc1, &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(b[4][startb+j]*dt/6.0, &inc1, &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(b[3][startb+j]*dt/5.0, &inc1, &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(b[2][startb+j]*dt/4.0, &inc1, &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(b[1][startb+j]*dt/3.0, &inc1, &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(b[0][startb+j]*dt/2.0, &inc1, &(xIntegCompCoeffs[starti+j+iterStep]));
            comp_sum(accInteg0[startb+j]*dt, &inc1, &(xIntegCompCoeffs[starti+j+iterStep]));
        } else {
            inc0 = ((((((((b[6][startb+j]*7.*h/9. + b[5][startb+j])*3.*h/4. + b[4][startb+j])*5.*h/7. + b[3][startb+j])*2.*h/3. + b[2][startb+j])*3.*h/5. + b[1][startb+j])*h/2.    + b[0][startb+j])*h/3. + accInteg0[startb+j])*dt*h/2. + xInteg0[starti+j+iterStep])*dt*h;;
            inc1 = (((((((b[6][startb+j]*7.*h/8.  + b[5][startb+j])*6.*h/7. + b[4][startb+j])*5.*h/6. + b[3][startb+j])*4.*h/5. + b[2][startb+j])*3.*h/4. + b[1][startb+j])*2.*h/3. + b[0][startb+j])*h/2. + accInteg0[startb+j])*dt*h;
        }
        xIntegNext[starti+j]          = xInteg0[starti+j] + inc0 - xIntegCompCoeffs[starti+j];
        xIntegNext[starti+j+iterStep] = xInteg0[starti+j+iterStep] + inc1 - xIntegCompCoeffs[starti+j+iterStep];
    }
}

void approx_xInteg(const std::vector<real> &xInteg0,
                   const std::vector<real> &accInteg0, const real &dt,
                   const real &h, const std::vector<std::vector<real>> &b,
                   const std::vector<IntegBody> &integBodies,
                   std::vector<real> &xIntegNext,
                   std::vector<real> &xIntegCompCoeffs) {
    size_t starti = 0;
    size_t startb = 0;
    for (size_t i = 0; i < integBodies.size(); i++) {
        // do pos/vel first, then STM
        approx_xInteg_math(xInteg0, accInteg0, dt, h, b, starti, startb, 3, xIntegNext, xIntegCompCoeffs);
        starti += 6;
        startb += 3;
        if (integBodies[i].propStm) {
            // do 6x6 STM first, then 6x1 fitted parameters
            approx_xInteg_math(xInteg0, accInteg0, dt, h, b, starti, startb, 18, xIntegNext, xIntegCompCoeffs);
            starti += 36;
            startb += 18;
            if (integBodies[i].stm.size() > 36) {
                // do fitted parameters
                const size_t numParams = (integBodies[i].stm.size()-36)/6;
                for (size_t param = 0; param < numParams; param++) {
                    approx_xInteg_math(xInteg0, accInteg0, dt, h, b, starti, startb, 3, xIntegNext, xIntegCompCoeffs);
                    starti += 6;
                    startb += 3;
                }
            }
        }
    }
}

void update_g_with_b(const std::vector<std::vector<real>> &b, const size_t &dim, real *g) {
    for (size_t i = 0; i < dim; i++) {
        g[0*dim+i] = b[6][i]*dVec[15] + b[5][i]*dVec[10] + b[4][i]*dVec[6] + b[3][i]*dVec[3] + b[2][i]*dVec[1] + b[1][i]*dVec[0] + b[0][i];
        g[1*dim+i] = b[6][i]*dVec[16] + b[5][i]*dVec[11] + b[4][i]*dVec[7] + b[3][i]*dVec[4] + b[2][i]*dVec[2] + b[1][i];
        g[2*dim+i] = b[6][i]*dVec[17] + b[5][i]*dVec[12] + b[4][i]*dVec[8] + b[3][i]*dVec[5] + b[2][i];
        g[3*dim+i] = b[6][i]*dVec[18] + b[5][i]*dVec[13] + b[4][i]*dVec[9] + b[3][i];
        g[4*dim+i] = b[6][i]*dVec[19] + b[5][i]*dVec[14] + b[4][i];
        g[5*dim+i] = b[6][i]*dVec[20] + b[5][i];
        g[6*dim+i] = b[6][i];
    }
}

void compute_g_and_b(const std::vector<std::vector<real>> &AccIntegArr,
                     const size_t &hIdx, real *g, real *bCompCoeffs,
                     std::vector<std::vector<real>> &b, const size_t &dim) {
    switch (hIdx) {
        case 0:
            throw std::runtime_error("hIdx cannot be 0");
        case 1:
            for (size_t i = 0; i < dim; i++) {
                // g[0*dim+i] = (AccIntegArr[1][i] - AccIntegArr[0][i]) * rMat[1][0];
                // b[0][i] = cMat[0][0]*g[0*dim+i] + cMat[1][0]*g[1*dim+i] + cMat[2][0]*g[2*dim+i] + cMat[3][0]*g[3*dim+i] + cMat[4][0]*g[4*dim+i] + cMat[5][0]*g[5*dim+i] + cMat[6][0]*g[6*dim+i];
                
                real dummy = 0;
                real temp = g[0*dim+i];
                real gVal = AccIntegArr[hIdx][i];
                comp_sum(-AccIntegArr[0][i], &gVal, &dummy);
                g[0*dim+i] = gVal/rVec[0];
                comp_sum(g[0*dim+i]-temp, &(b[0][i]), &(bCompCoeffs[0*dim+i]));
            }
            break;
        case 2:
            for (size_t i = 0; i < dim; i++) {
                // g[1*dim+i] = ((AccIntegArr[2][i] - AccIntegArr[0][i]) * rMat[2][0] - g[0*dim+i]) * rMat[2][1];
                // b[0][i] = cMat[0][0]*g[0*dim+i] + cMat[1][0]*g[1*dim+i] + cMat[2][0]*g[2*dim+i] + cMat[3][0]*g[3*dim+i] + cMat[4][0]*g[4*dim+i] + cMat[5][0]*g[5*dim+i] + cMat[6][0]*g[6*dim+i];
                // b[1][i] =                       + cMat[1][1]*g[1*dim+i] + cMat[2][1]*g[2*dim+i] + cMat[3][1]*g[3*dim+i] + cMat[4][1]*g[4*dim+i] + cMat[5][1]*g[5*dim+i] + cMat[6][1]*g[6*dim+i];
                
                real dummy = 0;
                real temp = g[1*dim+i];
                real gVal = AccIntegArr[hIdx][i];
                comp_sum(-AccIntegArr[0][i], &gVal, &dummy);
                g[1*dim+i] = (gVal/rVec[1] - g[0*dim+i])/rVec[2];
                temp = g[1*dim+i] - temp;
                comp_sum(temp*cVec[0], &(b[0][i]), &(bCompCoeffs[0*dim+i]));
                comp_sum(temp        , &(b[1][i]), &(bCompCoeffs[1*dim+i]));
            }
            break;
        case 3:
            for (size_t i = 0; i < dim; i++) {
                // g[2*dim+i] = (((AccIntegArr[3][i] - AccIntegArr[0][i]) * rMat[3][0] - g[0*dim+i]) * rMat[3][1] - g[1*dim+i]) * rMat[3][2];
                // b[0][i] = cMat[0][0]*g[0*dim+i] + cMat[1][0]*g[1*dim+i] + cMat[2][0]*g[2*dim+i] + cMat[3][0]*g[3*dim+i] + cMat[4][0]*g[4*dim+i] + cMat[5][0]*g[5*dim+i] + cMat[6][0]*g[6*dim+i];
                // b[1][i] =                       + cMat[1][1]*g[1*dim+i] + cMat[2][1]*g[2*dim+i] + cMat[3][1]*g[3*dim+i] + cMat[4][1]*g[4*dim+i] + cMat[5][1]*g[5*dim+i] + cMat[6][1]*g[6*dim+i];
                // b[2][i] =                                               + cMat[2][2]*g[2*dim+i] + cMat[3][2]*g[3*dim+i] + cMat[4][2]*g[4*dim+i] + cMat[5][2]*g[5*dim+i] + cMat[6][2]*g[6*dim+i];

                real dummy = 0;
                real temp = g[2*dim+i];
                real gVal = AccIntegArr[hIdx][i];
                comp_sum(-AccIntegArr[0][i], &gVal, &dummy);
                g[2*dim+i] = ((gVal/rVec[3] - g[0*dim+i])/rVec[4] - g[1*dim+i])/rVec[5];
                temp = g[2*dim+i] - temp;
                comp_sum(temp*cVec[1], &(b[0][i]), &(bCompCoeffs[0*dim+i]));
                comp_sum(temp*cVec[2], &(b[1][i]), &(bCompCoeffs[1*dim+i]));
                comp_sum(temp        , &(b[2][i]), &(bCompCoeffs[2*dim+i]));
            }
            break;
        case 4:
            for (size_t i = 0; i < dim; i++) {
                // g[3*dim+i] = ((((AccIntegArr[4][i] - AccIntegArr[0][i]) * rMat[4][0] - g[0*dim+i]) * rMat[4][1] - g[1*dim+i]) * rMat[4][2] - g[2*dim+i]) * rMat[4][3];
                // b[0][i] = cMat[0][0]*g[0*dim+i] + cMat[1][0]*g[1*dim+i] + cMat[2][0]*g[2*dim+i] + cMat[3][0]*g[3*dim+i] + cMat[4][0]*g[4*dim+i] + cMat[5][0]*g[5*dim+i] + cMat[6][0]*g[6*dim+i];
                // b[1][i] =                       + cMat[1][1]*g[1*dim+i] + cMat[2][1]*g[2*dim+i] + cMat[3][1]*g[3*dim+i] + cMat[4][1]*g[4*dim+i] + cMat[5][1]*g[5*dim+i] + cMat[6][1]*g[6*dim+i];
                // b[2][i] =                                               + cMat[2][2]*g[2*dim+i] + cMat[3][2]*g[3*dim+i] + cMat[4][2]*g[4*dim+i] + cMat[5][2]*g[5*dim+i] + cMat[6][2]*g[6*dim+i];
                // b[3][i] =                                                                       + cMat[3][3]*g[3*dim+i] + cMat[4][3]*g[4*dim+i] + cMat[5][3]*g[5*dim+i] + cMat[6][3]*g[6*dim+i];

                real dummy = 0;
                real temp = g[3*dim+i];
                real gVal = AccIntegArr[hIdx][i];
                comp_sum(-AccIntegArr[0][i], &gVal, &dummy);
                g[3*dim+i] = (((gVal/rVec[6] - g[0*dim+i])/rVec[7] - g[1*dim+i])/rVec[8] - g[2*dim+i])/rVec[9];
                temp = g[3*dim+i] - temp;
                comp_sum(temp*cVec[3], &(b[0][i]), &(bCompCoeffs[0*dim+i]));
                comp_sum(temp*cVec[4], &(b[1][i]), &(bCompCoeffs[1*dim+i]));
                comp_sum(temp*cVec[5], &(b[2][i]), &(bCompCoeffs[2*dim+i]));
                comp_sum(temp        , &(b[3][i]), &(bCompCoeffs[3*dim+i]));
            }
            break;
        case 5:
            for (size_t i = 0; i < dim; i++) {
                // g[4*dim+i] = (((((AccIntegArr[5][i] - AccIntegArr[0][i]) * rMat[5][0] - g[0*dim+i]) * rMat[5][1] - g[1*dim+i]) * rMat[5][2] - g[2*dim+i]) * rMat[5][3] - g[3*dim+i]) * rMat[5][4];
                // b[0][i] = cMat[0][0]*g[0*dim+i] + cMat[1][0]*g[1*dim+i] + cMat[2][0]*g[2*dim+i] + cMat[3][0]*g[3*dim+i] + cMat[4][0]*g[4*dim+i] + cMat[5][0]*g[5*dim+i] + cMat[6][0]*g[6*dim+i];
                // b[1][i] =                       + cMat[1][1]*g[1*dim+i] + cMat[2][1]*g[2*dim+i] + cMat[3][1]*g[3*dim+i] + cMat[4][1]*g[4*dim+i] + cMat[5][1]*g[5*dim+i] + cMat[6][1]*g[6*dim+i];
                // b[2][i] =                                               + cMat[2][2]*g[2*dim+i] + cMat[3][2]*g[3*dim+i] + cMat[4][2]*g[4*dim+i] + cMat[5][2]*g[5*dim+i] + cMat[6][2]*g[6*dim+i];
                // b[3][i] =                                                                       + cMat[3][3]*g[3*dim+i] + cMat[4][3]*g[4*dim+i] + cMat[5][3]*g[5*dim+i] + cMat[6][3]*g[6*dim+i];
                // b[4][i] =                                                                                               + cMat[4][4]*g[4*dim+i] + cMat[5][4]*g[5*dim+i] + cMat[6][4]*g[6*dim+i];

                real dummy = 0;
                real temp = g[4*dim+i];
                real gVal = AccIntegArr[hIdx][i];
                comp_sum(-AccIntegArr[0][i], &gVal, &dummy);
                g[4*dim+i] = ((((gVal/rVec[10] - g[0*dim+i])/rVec[11] - g[1*dim+i])/rVec[12] - g[2*dim+i])/rVec[13] - g[3*dim+i])/rVec[14];
                temp = g[4*dim+i] - temp;
                comp_sum(temp*cVec[6], &(b[0][i]), &(bCompCoeffs[0*dim+i]));
                comp_sum(temp*cVec[7], &(b[1][i]), &(bCompCoeffs[1*dim+i]));
                comp_sum(temp*cVec[8], &(b[2][i]), &(bCompCoeffs[2*dim+i]));
                comp_sum(temp*cVec[9], &(b[3][i]), &(bCompCoeffs[3*dim+i]));
                comp_sum(temp        , &(b[4][i]), &(bCompCoeffs[4*dim+i]));
            }
            break;
        case 6:
            for (size_t i = 0; i < dim; i++) {
                // g[5*dim+i] = ((((((AccIntegArr[6][i] - AccIntegArr[0][i]) * rMat[6][0] - g[0*dim+i]) * rMat[6][1] - g[1*dim+i]) * rMat[6][2] - g[2*dim+i]) * rMat[6][3] - g[3*dim+i]) * rMat[6][4] - g[4*dim+i]) * rMat[6][5];
                // b[0][i] = cMat[0][0]*g[0*dim+i] + cMat[1][0]*g[1*dim+i] + cMat[2][0]*g[2*dim+i] + cMat[3][0]*g[3*dim+i] + cMat[4][0]*g[4*dim+i] + cMat[5][0]*g[5*dim+i] + cMat[6][0]*g[6*dim+i];
                // b[1][i] =                       + cMat[1][1]*g[1*dim+i] + cMat[2][1]*g[2*dim+i] + cMat[3][1]*g[3*dim+i] + cMat[4][1]*g[4*dim+i] + cMat[5][1]*g[5*dim+i] + cMat[6][1]*g[6*dim+i];
                // b[2][i] =                                               + cMat[2][2]*g[2*dim+i] + cMat[3][2]*g[3*dim+i] + cMat[4][2]*g[4*dim+i] + cMat[5][2]*g[5*dim+i] + cMat[6][2]*g[6*dim+i];
                // b[3][i] =                                                                       + cMat[3][3]*g[3*dim+i] + cMat[4][3]*g[4*dim+i] + cMat[5][3]*g[5*dim+i] + cMat[6][3]*g[6*dim+i];
                // b[4][i] =                                                                                               + cMat[4][4]*g[4*dim+i] + cMat[5][4]*g[5*dim+i] + cMat[6][4]*g[6*dim+i];
                // b[5][i] =                                                                                                                       + cMat[5][5]*g[5*dim+i] + cMat[6][5]*g[6*dim+i];

                real dummy = 0;
                real temp = g[5*dim+i];
                real gVal = AccIntegArr[hIdx][i];
                comp_sum(-AccIntegArr[0][i], &gVal, &dummy);
                g[5*dim+i] = (((((gVal/rVec[15] - g[0*dim+i])/rVec[16] - g[1*dim+i])/rVec[17] - g[2*dim+i])/rVec[18] - g[3*dim+i])/rVec[19] - g[4*dim+i])/rVec[20];
                temp = g[5*dim+i] - temp;
                comp_sum(temp*cVec[10], &(b[0][i]), &(bCompCoeffs[0*dim+i]));
                comp_sum(temp*cVec[11], &(b[1][i]), &(bCompCoeffs[1*dim+i]));
                comp_sum(temp*cVec[12], &(b[2][i]), &(bCompCoeffs[2*dim+i]));
                comp_sum(temp*cVec[13], &(b[3][i]), &(bCompCoeffs[3*dim+i]));
                comp_sum(temp*cVec[14], &(b[4][i]), &(bCompCoeffs[4*dim+i]));
                comp_sum(temp         , &(b[5][i]), &(bCompCoeffs[5*dim+i]));
            }
            break;
        case 7:
            for (size_t i = 0; i < dim; i++) {
                // g[6*dim+i] = (((((((AccIntegArr[7][i] - AccIntegArr[0][i]) * rMat[7][0] - g[0*dim+i]) * rMat[7][1] - g[1*dim+i]) * rMat[7][2] - g[2*dim+i]) * rMat[7][3] - g[3*dim+i]) * rMat[7][4] - g[4*dim+i]) * rMat[7][5] - g[5*dim+i]) * rMat[7][6];
                // b[0][i] = cMat[0][0]*g[0*dim+i] + cMat[1][0]*g[1*dim+i] + cMat[2][0]*g[2*dim+i] + cMat[3][0]*g[3*dim+i] + cMat[4][0]*g[4*dim+i] + cMat[5][0]*g[5*dim+i] + cMat[6][0]*g[6*dim+i];
                // b[1][i] =                       + cMat[1][1]*g[1*dim+i] + cMat[2][1]*g[2*dim+i] + cMat[3][1]*g[3*dim+i] + cMat[4][1]*g[4*dim+i] + cMat[5][1]*g[5*dim+i] + cMat[6][1]*g[6*dim+i];
                // b[2][i] =                                               + cMat[2][2]*g[2*dim+i] + cMat[3][2]*g[3*dim+i] + cMat[4][2]*g[4*dim+i] + cMat[5][2]*g[5*dim+i] + cMat[6][2]*g[6*dim+i];
                // b[3][i] =                                                                       + cMat[3][3]*g[3*dim+i] + cMat[4][3]*g[4*dim+i] + cMat[5][3]*g[5*dim+i] + cMat[6][3]*g[6*dim+i];
                // b[4][i] =                                                                                               + cMat[4][4]*g[4*dim+i] + cMat[5][4]*g[5*dim+i] + cMat[6][4]*g[6*dim+i];
                // b[5][i] =                                                                                                                       + cMat[5][5]*g[5*dim+i] + cMat[6][5]*g[6*dim+i];
                // b[6][i] =                                                                                                                                               + cMat[6][6]*g[6*dim+i];

                real dummy = 0;
                real temp = g[6*dim+i];
                real gVal = AccIntegArr[hIdx][i];
                comp_sum(-AccIntegArr[0][i], &gVal, &dummy);
                g[6*dim+i] = ((((((gVal/rVec[21] - g[0*dim+i])/rVec[22] - g[1*dim+i])/rVec[23] - g[2*dim+i])/rVec[24] - g[3*dim+i])/rVec[25] - g[4*dim+i])/rVec[26] - g[5*dim+i])/rVec[27];
                temp = g[6*dim+i] - temp;
                comp_sum(temp*cVec[15], &(b[0][i]), &(bCompCoeffs[0*dim+i]));
                comp_sum(temp*cVec[16], &(b[1][i]), &(bCompCoeffs[1*dim+i]));
                comp_sum(temp*cVec[17], &(b[2][i]), &(bCompCoeffs[2*dim+i]));
                comp_sum(temp*cVec[18], &(b[3][i]), &(bCompCoeffs[3*dim+i]));
                comp_sum(temp*cVec[19], &(b[4][i]), &(bCompCoeffs[4*dim+i]));
                comp_sum(temp*cVec[20], &(b[5][i]), &(bCompCoeffs[5*dim+i]));
                comp_sum(temp         , &(b[6][i]), &(bCompCoeffs[6*dim+i]));
            }
            break;
        default:
            throw std::runtime_error("hIdx is out of range");
    }
}

void refine_b(std::vector<std::vector<real>> &b,
              real *e, const real &dtRatio,
              const size_t &dim) {
    std::vector<std::vector<real>> bDiff(7, std::vector<real>(dim, 0.0L));
    for (size_t i = 0; i < dim; i++){
        bDiff[0][i] = b[0][i] - e[0*dim+i];
        bDiff[1][i] = b[1][i] - e[1*dim+i];
        bDiff[2][i] = b[2][i] - e[2*dim+i];
        bDiff[3][i] = b[3][i] - e[3*dim+i];
        bDiff[4][i] = b[4][i] - e[4*dim+i];
        bDiff[5][i] = b[5][i] - e[5*dim+i];
        bDiff[6][i] = b[6][i] - e[6*dim+i];
    }

    real q = dtRatio;
    real q2 = q * q;
    real q3 = q2 * q;
    real q4 = q2 * q2;
    real q5 = q2 * q3;
    real q6 = q3 * q3;
    real q7 = q2 * q5;

    for (size_t i = 0; i < dim; i++) {
        e[0*dim+i] = q  * (b[6][i] * 7.0  + b[5][i] * 6.0  + b[4][i] * 5.0  + b[3][i] * 4.0 + b[2][i] * 3.0 + b[1][i] * 2.0 + b[0][i]);
        e[1*dim+i] = q2 * (b[6][i] * 21.0 + b[5][i] * 15.0 + b[4][i] * 10.0 + b[3][i] * 6.0 + b[2][i] * 3.0 + b[1][i]);
        e[2*dim+i] = q3 * (b[6][i] * 35.0 + b[5][i] * 20.0 + b[4][i] * 10.0 + b[3][i] * 4.0 + b[2][i]);
        e[3*dim+i] = q4 * (b[6][i] * 35.0 + b[5][i] * 15.0 + b[4][i] * 5.0  + b[3][i]);
        e[4*dim+i] = q5 * (b[6][i] * 21.0 + b[5][i] * 6.0  + b[4][i]);
        e[5*dim+i] = q6 * (b[6][i] * 7.0  + b[5][i]);
        e[6*dim+i] = q7 * (b[6][i]);
    }

    for (size_t i = 0; i < dim; i++) {
        b[0][i] = e[0*dim+i] + bDiff[0][i];
        b[1][i] = e[1*dim+i] + bDiff[1][i];
        b[2][i] = e[2*dim+i] + bDiff[2][i];
        b[3][i] = e[3*dim+i] + bDiff[3][i];
        b[4][i] = e[4*dim+i] + bDiff[4][i];
        b[5][i] = e[5*dim+i] + bDiff[5][i];
        b[6][i] = e[6*dim+i] + bDiff[6][i];
    }
}

void check_and_apply_events(propSimulation *propSim, const real &t,
                            real &tNextEvent, size_t &nextEventIdx,
                            std::vector<real> &xInteg) {
    while (nextEventIdx < propSim->events.size() && t == tNextEvent) {
        // apply events for the state just reached by the integrator
        real propDir;
        if (propSim->integParams.t0 < propSim->integParams.tf) {
            propDir = 1.0L;
        } else {
            propDir = -1.0L;
        }
        propSim->events[nextEventIdx].apply(t, xInteg, propDir);
        // update next event index and time
        nextEventIdx += 1;
        if (nextEventIdx < propSim->events.size()) {
            tNextEvent = propSim->events[nextEventIdx].t;
        } else {
            tNextEvent = propSim->integParams.tf;
        }
    }
}

static real root7(real num){
    real fac = 1.0;
    if (num<1e-7){
        fac = 0.1;
        num *= 1e7;
    }
    if (num>1e2){
        fac = 10.0;
        num *= 1e-7;
    }
    real root = 1.0;
    for (size_t i = 0; i < 20; i++) {
        root += (num/(root*root*root*root*root*root)-root)/7.0;
    }
    return fac*root;
}

void gr15(propSimulation *propSim) {
    if (!std::isfinite(propSim->t)) {
        throw std::runtime_error("t is not finite");
    }
    for (size_t i = 0; i < propSim->xInteg.size(); i++) {
        if (!std::isfinite(propSim->xInteg[i])) {
            throw std::runtime_error("xInteg is not finite");
        }
    }
    real t = propSim->t;
    std::vector<real> xInteg0 = propSim->xInteg;
    size_t nh = 8;
    size_t dim = propSim->integParams.n2Derivs;
    real dt = get_initial_timestep(propSim);
    propSim->integParams.timestepCounter = 0;
    std::vector<real> accInteg0 = get_state_der(t, xInteg0, propSim);
    std::vector<real> accIntegNext = std::vector<real>(accInteg0.size(), 0.0);
    std::vector<real> xInteg(propSim->xInteg.size(), 0.0);
    std::vector<real> xIntegCompCoeffs(propSim->xInteg.size(), 0.0);
    std::vector<std::vector<real> > bOld(7, std::vector<real>(dim, 0.0));
    std::vector<std::vector<real> > b(7, std::vector<real>(dim, 0.0));
    real *bCompCoeffs = new real [7 * dim];
    memset(bCompCoeffs, 0.0, 7 * dim * sizeof(real));
    real *g = new real[7 * dim];
    memset(g, 0.0, 7 * dim * sizeof(real));
    real *e = new real[7 * dim];
    memset(e, 0.0, 7 * dim * sizeof(real));
    std::vector<std::vector<real> > accIntegArr(nh,
                                                std::vector<real>(dim, 0.0));
    real *b6Tilde = new real[dim];
    memset(b6Tilde, 0.0, dim * sizeof(real));
    real b6TildeMax, accIntegArr7Max;
    real b6Max, accIntegNextMax;
    real relError, dtReq;
    real tNextEvent = propSim->integParams.tf;
    static size_t nextEventIdx = 0;
    if (t == propSim->integParams.t0) {
        nextEventIdx = 0;
    }
    if (propSim->events.size() != 0) {
        tNextEvent = propSim->events[0].t;
    }
    check_and_apply_events(propSim, t, tNextEvent, nextEventIdx, xInteg0);
    if ((propSim->integParams.tf > propSim->integParams.t0 &&
         t + dt > tNextEvent) ||
        (propSim->integParams.tf < propSim->integParams.t0 &&
         t + dt < tNextEvent)) {
        dt = tNextEvent - t;
    }
    size_t PCmaxIter = 12;
    int keepStepping = 1;
    int oneStepDone = 0;
    if (propSim->integParams.t0 == propSim->integParams.tf) {
        keepStepping = 0;
    }
    while (keepStepping) {
        t = propSim->t;
        xInteg0 = propSim->xInteg;
        oneStepDone = 0;
        update_g_with_b(b, dim, g);
        while (!oneStepDone) {
            real PCerr = 1.0/propSim->integParams.tolPC;
            real PCerrPrev = PCerr + 1.0;
            accIntegArr[0] = accInteg0;
            size_t PCIter = 0;
            while (true) {
                if (PCerr < propSim->integParams.tolPC) {
                    break;
                }
                if (PCerr > PCerrPrev && PCIter > 2) {
                    break;
                }
                if (PCIter > PCmaxIter) {
                    break;
                }
                PCerrPrev = PCerr;
                PCIter++;
                for (size_t hIdx = 1; hIdx < nh; hIdx++) {
                    approx_xInteg(xInteg0, accInteg0, dt, hVec[hIdx], b,
                                  propSim->integBodies, xInteg, xIntegCompCoeffs);
                    accIntegArr[hIdx] = get_state_der(
                        t + hVec[hIdx] * dt, xInteg, propSim);
                    compute_g_and_b(accIntegArr, hIdx, g, bCompCoeffs, b, dim);
                }
                for (size_t i = 0; i < dim; i++) {
                    b6Tilde[i] = b[6][i] - bOld[6][i];
                }
                vabs_max(b6Tilde, dim, b6TildeMax);
                vabs_max(accIntegArr[7], accIntegArr7Max);
                PCerr = b6TildeMax / accIntegArr7Max;
                bOld = b;
            }
            // printf("g0[0] = %0.6e, g0[1] = %0.6e, g0[2] = %0.6e\n",g[0], g[1], g[2]);
            // printf("g1[0] = %0.6e, g1[1] = %0.6e, g1[2] = %0.6e\n",g[3], g[4], g[5]);
            // printf("g2[0] = %0.6e, g2[1] = %0.6e, g2[2] = %0.6e\n",g[6], g[7], g[8]);
            // printf("g3[0] = %0.6e, g3[1] = %0.6e, g3[2] = %0.6e\n",g[9], g[10], g[11]);
            // printf("g4[0] = %0.6e, g4[1] = %0.6e, g4[2] = %0.6e\n",g[12], g[13], g[14]);
            // printf("g5[0] = %0.6e, g5[1] = %0.6e, g5[2] = %0.6e\n",g[15], g[16], g[17]);
            // printf("g6[0] = %0.6e, g6[1] = %0.6e, g6[2] = %0.6e\n",g[18], g[19], g[20]);
            // printf("b0[0] = %0.6e, b0[1] = %0.6e, b0[2] = %0.6e\n",b[0][0], b[0][1], b[0][2]);
            // printf("b1[0] = %0.6e, b1[1] = %0.6e, b1[2] = %0.6e\n",b[1][0], b[1][1], b[1][2]);
            // printf("b2[0] = %0.6e, b2[1] = %0.6e, b2[2] = %0.6e\n",b[2][0], b[2][1], b[2][2]);
            // printf("b3[0] = %0.6e, b3[1] = %0.6e, b3[2] = %0.6e\n",b[3][0], b[3][1], b[3][2]);
            // printf("b4[0] = %0.6e, b4[1] = %0.6e, b4[2] = %0.6e\n",b[4][0], b[4][1], b[4][2]);
            // printf("b5[0] = %0.6e, b5[1] = %0.6e, b5[2] = %0.6e\n",b[5][0], b[5][1], b[5][2]);
            // printf("b6[0] = %0.6e, b6[1] = %0.6e, b6[2] = %0.6e\n",b[6][0], b[6][1], b[6][2]);
            // printf("ks: %d. osd: %d. t: %0.10e, dt: %0.8e", keepStepping, oneStepDone,t, dt);
            // printf(". iterations: %zu",PCIter);
            // printf(". predictor_corrector_error: %e\n",PCerr);
            vabs_max(b[6], b6Max);
            vabs_max(accIntegArr[7], accIntegNextMax);
            relError = b6Max / accIntegNextMax;
            if (propSim->integParams.adaptiveTimestep) {
                if (std::isnormal(relError)) {
                    dtReq = root7(propSim->integParams.tolInteg/relError)*dt;
                } else {
                    dtReq = dt/propSim->integParams.dtChangeFactor;
                }
            } else {
                dtReq = dt;
            }
            // printf(". dtReq: %e. relError: %e", dtReq,relError);
            if (fabs(dtReq) < propSim->integParams.dtMin) {
                dtReq = copysign(propSim->integParams.dtMin, dtReq);
            }
            if (fabs(dtReq/dt) < propSim->integParams.dtChangeFactor) {
                if (propSim->integParams.timestepCounter > 1) {
                    refine_b(b, e, dtReq / dt, dim);
                }
                dt = dtReq;
                std::cout << "Restarting next while loop at time t = " << t << std::endl;
                continue;
            }
            // accept step
            if (dtReq / dt > 1.0 / propSim->integParams.dtChangeFactor) {
                dtReq = dt / propSim->integParams.dtChangeFactor;
            }
            propSim->interpParams.bStack.push_back(bOld);
            propSim->interpParams.accIntegStack.push_back(accInteg0);
            if (propSim->tEval.size() != propSim->xIntegEval.size()) {
                interpolate_on_the_fly(propSim, t, dt);
            }
            approx_xInteg(xInteg0, accInteg0, dt, 1.0, b,
                        propSim->integBodies, xInteg, xIntegCompCoeffs);
            accInteg0 = get_state_der(t + dt, xInteg, propSim);
            t += dt;
            propSim->t = t;
            check_and_apply_events(propSim, t, tNextEvent, nextEventIdx,
                                    xInteg);
            propSim->xInteg = xInteg;
            propSim->interpParams.tStack.push_back(t);
            propSim->interpParams.xIntegStack.push_back(xInteg);
            propSim->integParams.timestepCounter++;
            if (propSim->integParams.timestepCounter > 1) {
                refine_b(b, e, dtReq / dt, dim);
            }
            check_ca_or_impact(propSim, t-dt, xInteg0, t, xInteg, keepStepping);
            if ((propSim->integParams.tf > propSim->integParams.t0 &&
                    t >= propSim->integParams.tf) ||
                (propSim->integParams.tf < propSim->integParams.t0 &&
                    t <= propSim->integParams.tf)) {
                keepStepping = 0;
            }
            dt = dtReq;
            if ((propSim->integParams.tf > propSim->integParams.t0 &&
                 t + dt > tNextEvent) ||
                (propSim->integParams.tf < propSim->integParams.t0 &&
                 t + dt < tNextEvent)) {
                dt = tNextEvent - t;
            }
            oneStepDone = 1;
        }
    }
    delete[] bCompCoeffs;
    delete[] g;
    delete[] e;
    delete[] b6Tilde;
}
