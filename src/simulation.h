#ifndef SIMULATION_H
#define SIMULATION_H

typedef class Simulation
{
    private:

    public:
        // constants
        constants consts;
        // setup n-body variables
        nBodyParameters nbParams;
        // integration parameters
        integrationParameters integPrms;

        // setters
        void set_sim_params(const size_t nInteg, const size_t nSpice);
        void set_sim_constants(real du2m, real tu2sec, real G, real clight);
        void set_integration_parameters(real t0, real tf, real dt0=0.0L, real dtMax=6.0L, real dtChangeFactor=0.25L, bool adaptiveTimestep=true, real tolPC=1.0e-16L, real tolInteg=1.0e-6L);

        // getters
        std::vector<size_t> get_sim_params();
        std::vector<real> get_sim_constants();
} simulation;

#endif