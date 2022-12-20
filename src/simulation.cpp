#include "utilities.h"
#include "simulation.h"

void Simulation::set_sim_params(const size_t nInteg, const size_t nSpice){
    this->nbParams.nInteg = nInteg;
    this->nbParams.nSpice = nSpice;
    this->nbParams.nTotal = nInteg+nSpice;
    this->nbParams.dim = 6*nbParams.nTotal;
};

void Simulation::set_sim_constants(real du2m, real tu2sec, real G, real clight){
    this->consts.du2m = du2m;
    this->consts.tu2sec = tu2sec;
    this->consts.G = G;
    this->consts.clight = clight;
    this->consts.j2000Jd = 2451545.0;
    this->consts.JdMinusMjd = 2400000.5;
};

void Simulation::set_integration_parameters(real t0, real tf, real dt0, real dtMax, real dtChangeFactor, bool adaptiveTimestep, real tolPC, real tolInteg){
    this->integPrms.t0 = t0;
    this->integPrms.tf = tf;
    this->integPrms.dt0 = dt0;
    this->integPrms.dtMax = dtMax;
    this->integPrms.dtChangeFactor = dtChangeFactor;
    this->integPrms.adaptiveTimestep = adaptiveTimestep;
    this->integPrms.tolPC = tolPC;
    this->integPrms.tolInteg = tolInteg;
};

std::vector<size_t> Simulation::get_sim_params(){
    std::cout << "The number of integrated bodies is: " << nbParams.nInteg << std::endl;
    std::cout << "The number of bodies whose states are taken from SPICE kernels is: " << nbParams.nSpice << std::endl;
    std::cout << "The total number of bodies in this simulation is: " << nbParams.nTotal << std::endl;
    
    std::vector<size_t> params = {nbParams.nInteg, nbParams.nSpice, nbParams.nTotal};
    return params;
};

std::vector<real> Simulation::get_sim_constants(){
    std::cout << "The conversion from distance units to meters is: " << consts.du2m << std::endl;
    std::cout << "The conversion from time units to seconds is: " << consts.tu2sec << std::endl;
    std::cout << "The universal gravitational constant in nondimensional units is: " << consts.G << std::endl;
    std::cout << "The speed of light in nondimensional units is: " << consts.clight << std::endl;
    std::cout << "The Julian date of J2000 is: " << consts.j2000Jd << std::endl;
    std::cout << "The difference between the Julian date and the Modified Julian date is: " << consts.JdMinusMjd << std::endl;

    std::vector<real> constants = {consts.du2m, consts.tu2sec, consts.G, consts.clight, consts.j2000Jd, consts.JdMinusMjd};
    return constants;
};
