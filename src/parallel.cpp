#include "parallel.h"

std::vector<propSimulation> propSim_parallel_omp(
    const propSimulation refSim, const bool isCometary,
    const std::vector<std::vector<real> > &allBodies) {
    throw std::runtime_error("parallel.cpp: propSim_parallel_omp not fully implemented");
    size_t numBodies = allBodies.size();
    std::vector<propSimulation> allSims;

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
            NongravParamaters ngParams;
            ngParams.a1 = data[9];
            ngParams.a2 = data[10];
            ngParams.a3 = data[11];
            ngParams.alpha = data[12];
            ngParams.k = data[13];
            ngParams.m = data[14];
            ngParams.n = data[15];
            ngParams.r0_au = data[16];
            if (isCometary) {
                std::vector<real> com = {data[3], data[4], data[5], data[6], data[7], data[8]};
                IntegBody body(name+" body "+std::to_string(i), data[0], data[1],
                                data[2], com, ngParams);
                propSimulation sim(name, refSim);
                sim.add_integ_body(body);
                sim.integrate();
                #pragma omp critical
                {
                    allSims.push_back(sim);
                }
            } else {
                std::vector<real> pos = {data[3], data[4], data[5]};
                std::vector<real> vel = {data[6], data[7], data[8]};
                IntegBody body(name+" body "+std::to_string(i), data[0], data[1],
                                data[2], pos, vel, ngParams);
                propSimulation sim(name, refSim);
                sim.add_integ_body(body);
                sim.integrate();
                #pragma omp critical
                {
                    allSims.push_back(sim);
                }
            }
        }
    }
    return allSims;
}
