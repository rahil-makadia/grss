#include "parallel.h"

/**
 * @param[in] refSim Reference simulation to use as a template for the parallel propagation.
 * @param[in] isCometary Flag to indicate whether the bodies are cometary states or Cartesian states.
 * @param[in] allBodies Array of information for bodies to propagate in parallel.
 * @return std::vector<PropSimulation> Array of propagated simulations.
 */
std::vector<PropSimulation> propSim_parallel_omp(
    const PropSimulation refSim, const bool isCometary,
    const std::vector<std::vector<real> > &allBodies) {
    size_t numBodies = allBodies.size();
    std::vector<PropSimulation> allSims(numBodies, refSim);

    // make sure refSim is in parallel mode
    if (!refSim.parallelMode) {
        throw std::runtime_error(
            "ERROR: The reference simulation must be in parallel mode to "
            "propagate in parallel.");
    }

    // parallel for loop to first create an integBody for each entry in the
    // allBodies vector, then integrate each integBody using the reference
    // simulation
    omp_set_num_threads(omp_get_max_threads());
    #pragma omp parallel shared(allBodies, refSim)
    {
        #pragma omp for schedule(static)
        for (size_t i = 0; i < numBodies; i++) {
            std::vector<real> data = allBodies[i];
            std::string name = refSim.name+" clone "+std::to_string(i);
            NongravParameters ngParams;
            ngParams.a1 = data[9];
            ngParams.a2 = data[10];
            ngParams.a3 = data[11];
            ngParams.alpha = data[12];
            ngParams.k = data[13];
            ngParams.m = data[14];
            ngParams.n = data[15];
            ngParams.r0_au = data[16];
            PropSimulation sim(name, refSim);
            if (isCometary) {
                std::vector<real> com = {data[3], data[4], data[5],
                                         data[6], data[7], data[8]};
                IntegBody body(name, data[0], data[1], data[2], com, ngParams);
                sim.add_integ_body(body);
            } else {
                std::vector<real> pos = {data[3], data[4], data[5]};
                std::vector<real> vel = {data[6], data[7], data[8]};
                IntegBody body(name, data[0], data[1], data[2], pos, vel,
                               ngParams);
                sim.add_integ_body(body);
            }
            sim.integrate();
            allSims[i] = sim;
        }
    }
    // in serial, compute the body-fixed impact coordinates (lat,lon,alt) for each instance
    furnsh_c(refSim.DEkernelPath.c_str());
    for (size_t i = 0; i < allSims.size(); i++) {
        for (size_t j = 0; j < allSims[i].impactParams.size(); j++) {
            allSims[i].impactParams[j].get_impact_parameters(&(allSims[i]));
        }
    }
    unload_c(refSim.DEkernelPath.c_str());
    return allSims;
}
