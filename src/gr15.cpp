#include "gr15.h"

void gr15(real t, std::vector<real> xInteg, ForceParameters forceparams, IntegrationParameters integParams, Constants consts){
    std::vector<real> xDotInteg(6*integParams.nInteg, 0.0);
    get_state_der(t, xInteg, xDotInteg, forceparams, integParams, consts);
}
