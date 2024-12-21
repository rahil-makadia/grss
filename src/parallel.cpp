/**
 * @file    parallel.cpp
 * @brief   Source file for parallel propagation of bodies using OpenMP.
 * @author  Rahil Makadia <makadia2@illinois.edu>
 *
 * @section     LICENSE
 * Copyright (C) 2022-2025 Rahil Makadia
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, see <https://www.gnu.org/licenses>.
 */

#include "parallel.h"

/**
 * @param[in] refSim Reference simulation to use as a template for the parallel propagation.
 * @param[in] isCometary Flag to indicate whether the bodies are cometary states or Cartesian states.
 * @param[in] allBodies Array of information for bodies to propagate in parallel.
 * @param[in] maxThreads Maximum number of threads to use for parallel propagation.
 */
void propSim_parallel_omp(const PropSimulation refSim, const bool isCometary,
                          const std::vector<std::vector<real> > &allBodies,
                          const int &maxThreads) {
    const size_t numBodies = allBodies.size();
    // save directory is ref_sim.name with spaces replaced by underscores
    std::string saveDir = refSim.name;
    std::replace(saveDir.begin(), saveDir.end(), ' ', '_');

    // parallel for loop to first create an integBody for each entry in the
    // allBodies vector, then integrate each integBody using the reference
    // simulation
    int numThreads = omp_get_max_threads();
    numThreads = numThreads > maxThreads ? maxThreads : numThreads;
    if (numThreads < maxThreads) {
        std::cout << "Using " << numThreads << " threads instead of the requested "
                  << maxThreads << " threads." << std::endl;
    }
    #pragma omp parallel for schedule(static) shared(allBodies, refSim, saveDir) num_threads(numThreads)
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
        // save file name as refSim.name + "/" + i but i has leading zeros for the number of digits in numBodies
        std::string num = std::to_string(i);
        std::string zeros = std::string(std::to_string(numBodies-1).size()-num.size(), '0');
        std::string filename = "./logdir_"+saveDir+"/"+zeros+num+".log";
        sim.save(filename, true);
    }
}
