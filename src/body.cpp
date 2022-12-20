#include "utilities.h"
#include "body.h"

SpiceBody::SpiceBody(std::string name, int spiceID, real t0, real mass, real radius){
    this->name = std::to_string(spiceID) + " " + name;
    this->spiceID = spiceID;
    this->t0 = t0;
    this->mass = mass;
    this->radius = radius;
    this->isNongrav = false;
    if (this->isSpice){
        double state[6];
        double lt;
        get_spice_state_lt(this->spiceID, this->t0, state, lt);
        this->pos = {state[0], state[1], state[2]};
        this->vel = {state[3], state[4], state[5]};
    }
}

IntegBody::IntegBody(std::string name, real t0, real mass, real radius, std::vector<real> pos, std::vector<real> vel, std::vector< std::vector<real> > covariance, real a1, real a2, real a3, real alpha, real k, real m, real n, real r0_au){
    this->name = name;
    this->t0 = t0;
    this->mass = mass;
    this->radius = radius;
    this->pos = pos;
    this->vel = vel;
    this->covariance = covariance;
    if (a1!=0.0L || a2!=0.0L || a3!=0.0L){
        this->ngParams.a1 = a1;
        this->ngParams.a2 = a2;
        this->ngParams.a3 = a3;
        this->ngParams.alpha = alpha;
        this->ngParams.k = k;
        this->ngParams.m = m;
        this->ngParams.n = n;
        this->ngParams.r0_au = r0_au;
        this->isNongrav = true;
    }
}
