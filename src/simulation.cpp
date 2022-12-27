#include "simulation.h"

void Body::set_J2(real J2, real obliquityToEcliptic){
    this->J2 = J2;
    if (this->J2 != 0.0L){
        this->isJ2 = true;
    }
    this->obliquityToEcliptic = obliquityToEcliptic;
}

SpiceBody::SpiceBody(std::string name, int spiceId, real t0, real mass, real radius, Constants consts){
    this->name = std::to_string(spiceId) + " " + name;
    this->spiceId = spiceId;
    this->t0 = t0;
    this->mass = mass;
    this->radius = radius/consts.du2m;
    this->isNongrav = false;
    if (this->isSpice){
        double state[6];
        double lt;
        get_spice_state_lt(this->spiceId, this->t0, consts, state, lt);
        this->pos = {state[0], state[1], state[2]};
        this->vel = {state[3], state[4], state[5]};
    }
}

IntegBody::IntegBody(std::string name, real t0, real mass, real radius, std::vector<real> cometaryState, std::vector< std::vector<real> > covariance, NongravParams ngParams, Constants consts){
    this->name = name;
    this->t0 = t0;
    this->mass = mass;
    this->radius = radius/consts.du2m;
    std::vector<real> cartesianStateEclip(6);
    std::vector<real> cartesianPos(3);
    std::vector<real> cartesianVel(3);

    cometary_to_cartesian(t0, cometaryState, cartesianStateEclip, consts.G);
    // rotate to eme2000
    std::vector< std::vector<real> > eclipToEquatorial(3, std::vector<real>(3));
    rot_mat_x(-EARTH_OBLIQUITY, eclipToEquatorial); // rotate by negative earth obliquity for ecliptic to equatorial
    mat_vec_mul(eclipToEquatorial, {cartesianStateEclip[0],cartesianStateEclip[1],cartesianStateEclip[2]}, cartesianPos);
    mat_vec_mul(eclipToEquatorial, {cartesianStateEclip[3],cartesianStateEclip[4],cartesianStateEclip[5]}, cartesianVel);
    // shift heliocentric to barycentric
    double sunState[6];
    double lt;
    get_spice_state_lt(10, t0, consts, sunState, lt);
    for (size_t i=0; i<3; i++){
        cartesianPos[i] += sunState[i];
        cartesianVel[i] += sunState[i+3];
    }
    
    this->pos = cartesianPos;
    this->vel = cartesianVel;
    this->covariance = covariance;
    this->isNongrav = false;
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

IntegBody::IntegBody(std::string name, real t0, real mass, real radius, std::vector<real> pos, std::vector<real> vel, std::vector< std::vector<real> > covariance, NongravParams ngParams, Constants consts){
    this->name = name;
    this->t0 = t0;
    this->mass = mass;
    this->radius = radius/consts.du2m;
    this->pos = pos;
    this->vel = vel;
    this->covariance = covariance;
    this->isNongrav = false;
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


Simulation::Simulation(std::string name){
    this->name = name;
    this->integParams.nInteg = 0;
    this->integParams.nSpice = 0;
    this->integParams.nTotal = 0;
    // this->integParams.dim = 0;
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
    this->integParams.timestepCounter = 0;
    this->integParams.tolPC = 0.0L;
    this->integParams.tolInteg = 0.0L;
};

void Simulation::add_spice_body(std::string name, int spiceId, real t0, real mass, real radius, Constants consts){
    // check if body already exists. if so, throw error
    for (size_t i=0; i<this->spiceBodies.size(); i++){
        if (this->spiceBodies[i].name == name){
            throw std::invalid_argument("SPICE Body with name " + name + " already exists in simulation " + this->name);
        }
    }
    SpiceBody body(name, spiceId, t0, mass, radius, consts);
    this->spiceBodies.push_back(body);
    this->integParams.nSpice++;
    this->integParams.nTotal++;
    // this->integParams.dim = 3*integParams.nTotal;
};

void Simulation::add_spice_body(SpiceBody body){
    // check if body already exists. if so, throw error
    for (size_t i=0; i<this->spiceBodies.size(); i++){
        if (this->spiceBodies[i].name == body.name){
            throw std::invalid_argument("SPICE Body with name " + body.name + " already exists in simulation " + this->name);
        }
    }
    this->spiceBodies.push_back(body);
    this->integParams.nSpice++;
    this->integParams.nTotal++;
    // this->integParams.dim = 3*integParams.nTotal;
};

void Simulation::add_integ_body(std::string name, real t0, real mass, real radius, std::vector<real> cometaryState, std::vector< std::vector<real> > covariance, NongravParams ngParams, Constants consts){
    // check if body already exists. if so, throw error
    for (size_t i=0; i<this->integBodies.size(); i++){
        if (this->integBodies[i].name == name){
            throw std::invalid_argument("Integration body with name " + name + " already exists in simulation " + this->name);
        }
    }
    IntegBody body(name, t0, mass, radius, cometaryState, covariance, ngParams, consts);
    this->integBodies.push_back(body);
    this->integParams.nInteg++;
    this->integParams.nTotal++;
    // this->integParams.dim = 3*integParams.nTotal;
};

void Simulation::add_integ_body(std::string name, real t0, real mass, real radius, std::vector<real> pos, std::vector<real> vel, std::vector< std::vector<real> > covariance, NongravParams ngParams, Constants consts){
    // check if body already exists. if so, throw error
    for (size_t i=0; i<this->integBodies.size(); i++){
        if (this->integBodies[i].name == name){
            throw std::invalid_argument("Integration body with name " + name + " already exists in simulation " + this->name);
        }
    }
    IntegBody body(name, t0, mass, radius, pos, vel, covariance, ngParams, consts);
    this->integBodies.push_back(body);
    this->integParams.nInteg++;
    this->integParams.nTotal++;
    // this->integParams.dim = 3*integParams.nTotal;
};

void Simulation::add_integ_body(IntegBody body){
    // check if body already exists. if so, throw error
    for (size_t i=0; i<this->integBodies.size(); i++){
        if (this->integBodies[i].name == body.name){
            throw std::invalid_argument("Integration body with name " + body.name + " already exists in simulation " + this->name);
        }
    }
    this->integBodies.push_back(body);
    this->integParams.nInteg++;
    this->integParams.nTotal++;
    // this->integParams.dim = 3*integParams.nTotal;
};

void Simulation::remove_body(std::string name){
    for (size_t i=0; i<this->spiceBodies.size(); i++){
        if (this->spiceBodies[i].name == name){
            this->spiceBodies.erase(this->spiceBodies.begin()+i);
            this->integParams.nSpice--;
            this->integParams.nTotal--;
            // this->integParams.dim = 3*integParams.nTotal;
            return;
        }
    }
    for (size_t i=0; i<this->integBodies.size(); i++){
        if (this->integBodies[i].name == name){
            this->integBodies.erase(this->integBodies.begin()+i);
            this->integParams.nInteg--;
            this->integParams.nTotal--;
            // this->integParams.dim = 3*integParams.nTotal;
            return;
        }
    }
    std::cout << "Error: Body " << name << " not found." << std::endl;
};

void Simulation::set_sim_constants(real du2m, real tu2sec, real G, real clight){
    this->consts.du2m = du2m;
    this->consts.tu2sec = tu2sec;
    this->consts.G = G;
    this->consts.clight = clight;
    this->consts.j2000Jd = 2451545.0;
    this->consts.JdMinusMjd = 2400000.5;
};

void Simulation::set_integration_parameters(real t0, real tf, bool adaptiveTimestep, real dt0, real dtMax, real dtMin, real dtChangeFactor, real tolInteg, real tolPC){
    this->integParams.t0 = t0;
    this->integParams.tf = tf;
    this->integParams.dt0 = dt0;
    this->integParams.dtMax = dtMax;
    this->integParams.dtMin = dtMin;
    this->integParams.dtChangeFactor = dtChangeFactor;
    this->integParams.adaptiveTimestep = adaptiveTimestep;
    this->integParams.tolPC = tolPC;
    this->integParams.tolInteg = tolInteg;
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
    // std::cout << "The number of integrated bodies is: " << integParams.nInteg << std::endl;
    // std::cout << "The number of bodies whose states are taken from SPICE kernels is: " << integParams.nSpice << std::endl;
    // std::cout << "The total number of bodies in this simulation is: " << integParams.nTotal << std::endl;
    // std::cout << "The initial time of the simulation is: " << integParams.t0 << std::endl;
    // std::cout << "The final time of the simulation is: " << integParams.tf << std::endl;
    // std::cout << "The initial timestep of the simulation is: " << integParams.dt0 << std::endl;
    // std::cout << "The maximum timestep of the simulation is: " << integParams.dtMax << std::endl;
    // std::cout << "The factor by which the timestep is changed is: " << integParams.dtChangeFactor << std::endl;
    // std::cout << "The adaptive timestep flag is: " << integParams.adaptiveTimestep << std::endl;
    // std::cout << "The tolerance for the position and velocity error is: " << integParams.tolPC << std::endl;
    // std::cout << "The tolerance for the integration error is: " << integParams.tolInteg << std::endl;

    std::vector<real> integration_parameters = {(real) integParams.nInteg, (real) integParams.nSpice, (real) integParams.nTotal, integParams.t0, integParams.tf, integParams.dt0, integParams.dtMax, integParams.dtMin, integParams.dtChangeFactor, (real) integParams.adaptiveTimestep, integParams.tolPC, integParams.tolInteg};
    return integration_parameters;
};

void Simulation::preprocess(){
    this->t = this->integParams.t0;
    for (size_t i = 0; i < this->integParams.nInteg; i++){
        for (size_t j = 0; j < 3; j++){
            this->xInteg.push_back(integBodies[i].pos[j]);
        }
        for (size_t j = 0; j < 3; j++){
            this->xInteg.push_back(integBodies[i].vel[j]);
        }
        this->forceParams.masses.push_back(integBodies[i].mass);
        this->forceParams.radii.push_back(integBodies[i].radius);
        this->forceParams.spiceIdList.push_back(-99999);
        this->forceParams.isPPNList.push_back(integBodies[i].isPPN);
        this->forceParams.isJ2List.push_back(integBodies[i].isJ2);
        this->forceParams.J2List.push_back(integBodies[i].J2);
        this->forceParams.obliquityList.push_back(integBodies[i].obliquityToEcliptic);
        this->forceParams.isNongravList.push_back(integBodies[i].isNongrav);
        this->forceParams.ngParamsList.push_back(integBodies[i].ngParams);
    }
    for (size_t i = 0; i < this->integParams.nSpice; i++){
        this->forceParams.masses.push_back(spiceBodies[i].mass);
        this->forceParams.radii.push_back(spiceBodies[i].radius);
        this->forceParams.spiceIdList.push_back(spiceBodies[i].spiceId);
        this->forceParams.isPPNList.push_back(spiceBodies[i].isPPN);
        this->forceParams.isJ2List.push_back(spiceBodies[i].isJ2);
        this->forceParams.J2List.push_back(spiceBodies[i].J2);
        this->forceParams.obliquityList.push_back(spiceBodies[i].obliquityToEcliptic);
    }
}
