#ifndef FORCE_H
#define FORCE_H
#include "utilities.h"

struct ForceParameters{
    std::vector<real> masses;
    std::vector<real> radii;
    std::vector<int> spiceIdList;
    std::vector<NongravParams> ngParamsList;
    std::vector<bool> isPPNList;
    std::vector<bool> isJ2List;
    std::vector<real> J2List;
    std::vector<real> obliquityList;
    std::vector<bool> isNongravList;
};

std::vector<real> get_state_der(const real &t, const std::vector<real> &xInteg, const ForceParameters &forceParams, const IntegrationParameters &integParams, const Constants &consts);
void force_newton(const std::vector<real> &posAll, std::vector<real> &xDotInteg, const ForceParameters &forceParams, const IntegrationParameters &integParams, const Constants &consts);
void force_ppn(const std::vector<real> &posAll, const std::vector<real> &velAll, std::vector<real> &xDotInteg, const ForceParameters &forceParams, const IntegrationParameters &integParams, const Constants &consts);
void force_J2(const std::vector<real> &posAll, std::vector<real> &xDotInteg, const ForceParameters &forceParams, const IntegrationParameters &integParams, const Constants &consts);
void force_nongrav(const std::vector<real> &posAll, const std::vector<real> &velAll, std::vector<real> &xDotInteg, const ForceParameters &forceParams, const IntegrationParameters &integParams, const Constants &consts);
#endif