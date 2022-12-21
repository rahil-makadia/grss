#include "utilities.h"
#include "simulation.h"

SpiceBody::SpiceBody(std::string name, int spiceID, real t0, real mass, real radius, Constants consts){
    this->name = std::to_string(spiceID) + " " + name;
    this->spiceID = spiceID;
    this->t0 = t0;
    this->mass = mass;
    this->radius = radius;
    this->isNongrav = false;
    if (this->isSpice){
        double state[6];
        double lt;
        get_spice_state_lt(this->spiceID, this->t0, consts, state, lt);
        this->pos = {state[0], state[1], state[2]};
        this->vel = {state[3], state[4], state[5]};
    }
}

IntegBody::IntegBody(std::string name, real t0, real mass, real radius, std::vector<real> cometaryState, std::vector< std::vector<real> > covariance, NongravParams ngParams, Constants consts){
    this->name = name;
    this->t0 = t0;
    this->mass = mass;
    this->radius = radius;
    this->pos = pos;
    this->vel = vel;
    this->covariance = covariance;
    if (ngParams.a1!=0.0L || ngParams.a2!=0.0L || ngParams.a3!=0.0L){
        this->ngParams.a1 = ngParams.a1;
        this->ngParams.a2 = ngParams.a2;
        this->ngParams.a3 = ngParams.a3;
        this->ngParams.alpha = ngParams.alpha;
        this->ngParams.k = ngParams.k;
        this->ngParams.m = ngParams.m;
        this->ngParams.n = ngParams.n;
        this->ngParams.r0_au = ngParams.r0_au;
        this->isNongrav = true;
    }
}


Simulation::Simulation(){
    this->nbParams.nInteg = 0;
    this->nbParams.nSpice = 0;
    this->nbParams.nTotal = 0;
    this->nbParams.dim = 0;
    this->consts.du2m = 0.0L;
    this->consts.tu2sec = 0.0L;
    this->consts.G = 0.0L;
    this->consts.clight = 0.0L;
    this->consts.j2000Jd = 0.0L;
    this->consts.JdMinusMjd = 0.0L;
    this->integParams.t0 = 0.0L;
    this->integParams.tf = 0.0L;
    this->integParams.dt0 = 0.0L;
    this->integParams.dtMax = 0.0L;
    this->integParams.dtChangeFactor = 0.0L;
    this->integParams.adaptiveTimestep = false;
    this->integParams.tolPC = 0.0L;
    this->integParams.tolInteg = 0.0L;
};

void Simulation::add_spice_body(std::string name, int spiceID, real t0, real mass, real radius, Constants consts){
    SpiceBody body(name, spiceID, t0, mass, radius, consts);
    this->spiceBodies.push_back(body);
    this->nbParams.nSpice++;
    this->nbParams.nTotal++;
    this->nbParams.dim = 3*nbParams.nTotal;
};

void Simulation::add_spice_body(SpiceBody body){
    this->spiceBodies.push_back(body);
    this->nbParams.nSpice++;
    this->nbParams.nTotal++;
    this->nbParams.dim = 3*nbParams.nTotal;
};

void Simulation::add_integ_body(std::string name, real t0, real mass, real radius, std::vector<real> cometaryState, std::vector< std::vector<real> > covariance, NongravParams ngParams, Constants consts){
    IntegBody body(name, t0, mass, radius, cometaryState, covariance, ngParams, consts);
    this->integBodies.push_back(body);
    this->nbParams.nInteg++;
    this->nbParams.nTotal++;
    this->nbParams.dim = 3*nbParams.nTotal;
    this->nbParams.nOde = 6*nbParams.nInteg;
};

void Simulation::add_integ_body(IntegBody body){
    this->integBodies.push_back(body);
    this->nbParams.nInteg++;
    this->nbParams.nTotal++;
    this->nbParams.dim = 3*nbParams.nTotal;
    this->nbParams.nOde = 6*nbParams.nInteg;
};

void Simulation::remove_body(std::string name){
    for (size_t i=0; i<this->spiceBodies.size(); i++){
        if (this->spiceBodies[i].name == name){
            this->spiceBodies.erase(this->spiceBodies.begin()+i);
            this->nbParams.nSpice--;
            this->nbParams.nTotal--;
            this->nbParams.dim = 3*nbParams.nTotal;
            return;
        }
    }
    for (size_t i=0; i<this->integBodies.size(); i++){
        if (this->integBodies[i].name == name){
            this->integBodies.erase(this->integBodies.begin()+i);
            this->nbParams.nInteg--;
            this->nbParams.nTotal--;
            this->nbParams.dim = 3*nbParams.nTotal;
            this->nbParams.nOde = 6*nbParams.nInteg;
            return;
        }
    }
    std::cout << "Error: Body " << name << " not found." << std::endl;
};

// void Simulation::set_sim_params(const size_t nInteg, const size_t nSpice){
//     this->nbParams.nInteg = nInteg;
//     this->nbParams.nSpice = nSpice;
//     this->nbParams.nTotal = nInteg+nSpice;
//     this->nbParams.dim = 3*nbParams.nTotal;
//     this->nbParams.nOde = 6*nbParams.nInteg;
// };

void Simulation::set_sim_constants(real du2m, real tu2sec, real G, real clight){
    this->consts.du2m = du2m;
    this->consts.tu2sec = tu2sec;
    this->consts.G = G;
    this->consts.clight = clight;
    this->consts.j2000Jd = 2451545.0;
    this->consts.JdMinusMjd = 2400000.5;
};

void Simulation::set_integration_parameters(real t0, real tf, real dt0, real dtMax, real dtChangeFactor, bool adaptiveTimestep, real tolPC, real tolInteg){
    this->integParams.t0 = t0;
    this->integParams.tf = tf;
    this->integParams.dt0 = dt0;
    this->integParams.dtMax = dtMax;
    this->integParams.dtChangeFactor = dtChangeFactor;
    this->integParams.adaptiveTimestep = adaptiveTimestep;
    this->integParams.tolPC = tolPC;
    this->integParams.tolInteg = tolInteg;
};

std::vector<size_t> Simulation::get_sim_params(){
    // std::cout << "The number of integrated bodies is: " << nbParams.nInteg << std::endl;
    // std::cout << "The number of bodies whose states are taken from SPICE kernels is: " << nbParams.nSpice << std::endl;
    // std::cout << "The total number of bodies in this simulation is: " << nbParams.nTotal << std::endl;
    
    std::vector<size_t> params = {nbParams.nInteg, nbParams.nSpice, nbParams.nTotal};
    return params;
};

std::vector<real> Simulation::get_sim_constants(){
    // std::cout << "The conversion from distance units to meters is: " << consts.du2m << std::endl;
    // std::cout << "The conversion from time units to seconds is: " << consts.tu2sec << std::endl;
    // std::cout << "The universal gravitational constant in nondimensional units is: " << consts.G << std::endl;
    // std::cout << "The speed of light in nondimensional units is: " << consts.clight << std::endl;
    // std::cout << "The Julian date of J2000 is: " << consts.j2000Jd << std::endl;
    // std::cout << "The difference between the Julian date and the Modified Julian date is: " << consts.JdMinusMjd << std::endl;

    std::vector<real> constants = {consts.du2m, consts.tu2sec, consts.G, consts.clight, consts.j2000Jd, consts.JdMinusMjd};
    return constants;
};

std::vector<real> Simulation::get_integration_parameters(){
    // std::cout << "The initial time of the simulation is: " << integParams.t0 << std::endl;
    // std::cout << "The final time of the simulation is: " << integParams.tf << std::endl;
    // std::cout << "The initial timestep of the simulation is: " << integParams.dt0 << std::endl;
    // std::cout << "The maximum timestep of the simulation is: " << integParams.dtMax << std::endl;
    // std::cout << "The factor by which the timestep is changed is: " << integParams.dtChangeFactor << std::endl;
    // std::cout << "The adaptive timestep flag is: " << integParams.adaptiveTimestep << std::endl;
    // std::cout << "The tolerance for the position and velocity error is: " << integParams.tolPC << std::endl;
    // std::cout << "The tolerance for the integration error is: " << integParams.tolInteg << std::endl;

    std::vector<real> integration_parameters = {integParams.t0, integParams.tf, integParams.dt0, integParams.dtMax, integParams.dtChangeFactor, (real) integParams.adaptiveTimestep, integParams.tolPC, integParams.tolInteg};
    return integration_parameters;
};
