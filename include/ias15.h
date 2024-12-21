/**
 * @file    ias15.h
 * @brief   Header file for the IAS15 integrator.
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

#ifndef IAS15_H
#define IAS15_H

#include "approach.h"

/**
 * @brief Vector of nodes for the Gauss-Radau quadrature.
 */
const real hVec[8] = {
    0.0,
    0.0562625605369221464656521910318,
    0.180240691736892364987579942780,
    0.352624717113169637373907769648,
    0.547153626330555383001448554766,
    0.734210177215410531523210605558,
    0.885320946839095768090359771030,
    0.977520613561287501891174488626
};

/**
 * @brief Vector of coefficients for computing the g-matrix of interpolation coefficients.
 */
const real rVec[28] = {
    0.0562625605369221464656522,
    0.1802406917368923649875799, 0.1239781311999702185219278,
    0.3526247171131696373739078, 0.2963621565762474909082556, 0.1723840253762772723863278,
    0.5471536263305553830014486, 0.4908910657936332365357964, 0.3669129345936630180138686, 0.1945289092173857456275408,
    0.7342101772154105315232106, 0.6779476166784883850575584, 0.5539694854785181665356307, 0.3815854601022408941493028, 0.1870565508848551485217621,
    0.8853209468390957680903598, 0.8290583863021736216247076, 0.7050802551022034031027798, 0.5326962297259261307164520, 0.3381673205085403850889112, 0.1511107696236852365671492,
    0.9775206135612875018911745, 0.9212580530243653554255223, 0.7972799218243951369035945, 0.6248958964481178645172667, 0.4303669872307321188897259, 0.2433104363458769703679639, 0.0921996667221917338008147
};

/**
 * @brief Vector of coefficients for computing the b-matrix of interpolation coefficients.
 */
const real cVec[21] = {
    -0.0562625605369221464656522,
    0.01014080283006362998648180399549641417413495311078, -0.2365032522738145114532321,
    -0.0035758977292516175949344589284567187362040464593728, 0.09353769525946206589574845561035371499343547051116, -0.5891279693869841488271399,
    0.0019565654099472210769005672379668610648179838140913, -0.054755386889068686440808430671055022602028382584495, 0.41588120008230686168862193041156933067050816537030, -1.1362815957175395318285885,
    -0.0014365302363708915424459554194153247134438571962198, 0.042158527721268707707297347813203202980228135395858, -0.36009959650205681228976647408968845289781580280782, 1.2501507118406910258505441186857527694077565516084, -1.8704917729329500633517991,
    0.001271790309026867749294311762296422088948466147501, -0.038760357915906770369904626849901899108502158354383, 0.36096224345284598322533983078129066420907893718190, -1.4668842084004269643701553461378480148761655599754, 2.9061362593084293014237914371173946705384212479246, -2.7558127197720458314421589
};

/**
 * @brief Vector of coefficients for updating the g-matrix of interpolation coefficients.
 */
const real dVec[21] = {
    0.0562625605369221464656522,
    0.0031654757181708292499905, 0.2365032522738145114532321,
    0.0001780977692217433881125, 0.0457929855060279188954539, 0.5891279693869841488271399,
    0.0000100202365223291272096, 0.0084318571535257015445000, 0.2535340690545692665214616, 1.1362815957175395318285885,
    0.0000005637641639318207610, 0.0015297840025004658189490, 0.0978342365324440053653648, 0.8752546646840910912297246, 1.8704917729329500633517991,
    0.0000000317188154017613665, 0.0002762930909826476593130, 0.0360285539837364596003871, 0.5767330002770787313544596, 2.2485887607691597933926895, 2.7558127197720458314421588
};

/**
 * @brief Compute the initial timestep for the PropSimulation.
 */
real get_initial_timestep(PropSimulation *propSim);

/**
 * @brief Update the g-matrix of interpolation coefficients using the b-matrix.
 */
void update_g_with_b(const std::vector<real> &b, const size_t &dim, std::vector<real> &g);

/**
 * @brief Compute the interpolation coefficients for the integration.
 */
void compute_g_and_b(const std::vector<std::vector<real>> &AccIntegArr,
                     const size_t &hIdx, std::vector<real> &g, std::vector<real> &bCompCoeffs,
                     std::vector<real> &b, const size_t &dim,
                     real &PCerr);

/**
 * @brief Refine the b-matrix of interpolation coefficients for the next timestep.
 */
void refine_b(std::vector<real> &b, std::vector<real> &e, const real &dtRatio,
              const size_t &dim);

/**
 * @brief Check if any impulsive events need to be applied after the current timestep.
 */
void check_and_apply_impulsive_events(PropSimulation *propSim, const real &t,
                                      std::vector<real> &xInteg);

/**
 * @brief Check if any continuous events need to be applied after the current timestep.
 */
void check_continuous_events(PropSimulation *propSim, const real &t);

/**
 * @brief Top level function to check for events after the current timestep.
 */
void check_events(PropSimulation *propSim, const real &t, std::vector<real> &xInteg);

/**
 * @brief Check whether timestep is too large that it would skip an event.
 */
void event_timestep_check(PropSimulation *propSim, real &dt);

/**
 * @brief 15th-order Gauss-Radau integrator for the PropSimulation.
 */
void ias15(PropSimulation *propSim);

#endif
