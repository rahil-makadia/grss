#ifndef SIMULATION_H
#define SIMULATION_H

class Body{
    private:

    public:
        real t0;
        real mass;
        real radius;
        std::string name;
        std::vector<real> pos;
        std::vector<real> vel;
        bool isNongrav=false;
};

class SpiceBody: public Body{
    private:

    public:
        int spiceID;
        bool isSpice=true;
        // constructor
        SpiceBody(std::string name, int spiceID, real t0, real mass, real radius, Constants consts);

};

class IntegBody: public Body{
    private:

    public:
        bool isInteg=true;
        std::vector< std::vector<real> > covariance;
        NongravParams ngParams;
        // constructor
        IntegBody(std::string name, real t0, real mass, real radius, std::vector<real> cometaryState, std::vector< std::vector<real> > covariance, NongravParams ngParams, Constants consts);

};

class Simulation
{
    private:

    public:
        // constructor
        Simulation();

        // constants
        Constants consts;
        // setup n-body variables
        NBodyParameters nbParams;
        // integration parameters
        IntegrationParameters integParams;

        // bodies
        std::vector<SpiceBody> spiceBodies;
        std::vector<IntegBody> integBodies;

        // add and remove bodies
        void add_spice_body(std::string name, int spiceID, real t0, real mass, real radius, Constants consts);
        void add_spice_body(SpiceBody body);
        void add_integ_body(std::string name, real t0, real mass, real radius, std::vector<real> cometaryState, std::vector< std::vector<real> > covariance, NongravParams ngParams, Constants consts);
        void add_integ_body(IntegBody body);
        void remove_body(std::string name);

        // setters
        // void set_sim_params(const size_t nInteg, const size_t nSpice);
        void set_sim_constants(real du2m=149597870700.0L, real tu2sec=86400.0L, real G=6.6743e-11L/(149597870700.0L*149597870700.0L*149597870700.0L)*86400.0L*86400.0L, real clight=299792458.0L/149597870700.0L*86400.0L);
        void set_integration_parameters(real t0, real tf, real dt0=0.0L, real dtMax=6.0L, real dtChangeFactor=0.25L, bool adaptiveTimestep=true, real tolPC=1.0e-16L, real tolInteg=1.0e-6L);

        // getters
        std::vector<size_t> get_sim_params();
        std::vector<real> get_sim_constants();
        std::vector<real> get_integration_parameters();
};

#endif