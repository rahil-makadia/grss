#ifndef GRSS_H
#define GRSS_H

#include "gr15.h"

void sort_and_clean_up_tEval(Simulation &sim){
    bool forwardProp = sim.integParams.t0 < sim.integParams.tf;
    bool backwardProp = sim.integParams.t0 > sim.integParams.tf;
    if (forwardProp){
        std::sort(sim.tEval.begin(), sim.tEval.end()); // sort sim.tEval into ascending order
        int removeCounter = 0;
        while (sim.tEval[0] < sim.integParams.t0 - sim.tEvalMargin){
            // remove any tEval values that are too small
            sim.tEval.erase(sim.tEval.begin());
            removeCounter++;
        }
        while (sim.tEval.back() > sim.integParams.tf + sim.tEvalMargin){
            // remove any tEval values that are too large
            sim.tEval.pop_back();
            removeCounter++;
        }
        if (removeCounter > 0){
            std::cout << "WARNING: " << removeCounter << " tEval value(s) were removed because they were outside the interpolation range, i.e., integration range with a margin of "<< sim.tEvalMargin << " day(s)." << std::endl;
        }
    }
    else if (backwardProp){
        std::sort(sim.tEval.begin(), sim.tEval.end(), std::greater<real>()); // sort sim.tEval into descending order
        int removeCounter = 0;
        while (sim.tEval[0] > sim.integParams.t0 + sim.tEvalMargin){
            // remove any tEval values that are too large
            sim.tEval.erase(sim.tEval.begin());
            removeCounter++;
        }
        while (sim.tEval.back() < sim.integParams.tf - sim.tEvalMargin){
            // remove any tEval values that are too small
            sim.tEval.pop_back();
            removeCounter++;
        }
        if (removeCounter > 0){
            std::cout << "WARNING: " << removeCounter << " tEval value(s) were removed because they were outside the interpolation range, i.e., integration range with a margin of "<< sim.tEvalMargin << " day(s)." << std::endl;
        }
    }
}

void Simulation::integrate(){
    // check that xInteg has a valid size. If not, raise error
    if (this->integParams.nInteg < 1){
        throw std::runtime_error("\n\ngrss.h: ERROR: There are no integration bodies in the simulation. Need at least one body with a full state vector to integrate.\n");
    }
    sort_and_clean_up_tEval(*this);
    bool backwardProp = this->integParams.t0 > this->integParams.tf;
    if (backwardProp){
        std::reverse(this->events.begin(),this->events.end());
    }
    furnsh_c(this->DEkernelPath.c_str());
    gr15(this->t, this->xInteg, *this);
    unload_c(this->DEkernelPath.c_str());
}
#endif
