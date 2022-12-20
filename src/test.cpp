/* COMPILE USING
g++ utilities.cpp simulation.cpp body.cpp test.cpp -o test.out -std=c++11 -fdiagnostics-color=always -I /opt/homebrew/include -lm /opt/homebrew/Cellar/cspice/67/lib/cspice.a -g3 -Wall
*/
#include "utilities.h"
#include "body.h"
#include "simulation.h"


int main() {
    std::cout.precision(22);
    std::cout << "The size on memory of a real variable is: " << sizeof(real) << std::endl;
    
    // simulation testSim0;
    // const size_t ni = 1;
    // const size_t ns = 27;
    // testSim0.set_sim_params(ni, ns);
    // testSim0.set_sim_constants();

    // std::vector<size_t> params = testSim0.get_sim_params();
    // std::vector<real> consts = testSim0.get_sim_constants();

    // real jd = testSim0.consts.j2000Jd;
    // real mjd=0;
    // real et=-10;
    // std::cout<< mjd << std::endl;
    // std::cout<< et << std::endl;
    // jd_to_mjd(jd, mjd);
    // jd_to_et(jd, et);
    // std::cout<< mjd << std::endl;
    // std::cout<< et << std::endl;

    // real et = 0;
    // real mjd=-12091;
    // real jd=2190938;
    // std::cout<< mjd << std::endl;
    // std::cout<< jd << std::endl;
    // et_to_mjd(et, mjd);
    // et_to_jd(et, jd);
    // std::cout<< mjd << std::endl;
    // std::cout<< jd << std::endl;

    // real mjd = 51544.5;
    // real jd = mjd_to_jd(mjd);
    // real et = mjd_to_et(mjd);
    // std::cout<< jd << std::endl;
    // std::cout<< et << std::endl;
    // std::cout<< mjd_to_jd(mjd) << std::endl;
    // std::cout<< mjd_to_et(mjd) << std::endl;
    // jd = 2938476;
    // et = 238974;
    // std::cout<< jd << std::endl;
    // std::cout<< et << std::endl;
    // mjd_to_jd(mjd, jd);
    // mjd_to_et(mjd, et);
    // std::cout<< jd << std::endl;
    // std::cout<< et << std::endl;

    furnsh_c ( "../spice/planets_big16_de441_1950_2350.bsp" );
    real t0Sim = 0.0L;
    spiceBody Sun("Sun", 10, t0Sim, 1.989e30L, 6.957e5L);
    std::cout << "Sun name: " << Sun.name << std::endl;
    std::cout << "Sun spiceID: " << Sun.spiceID << std::endl;
    std::cout << "Sun t0: " << Sun.t0 << std::endl;
    std::cout << "Sun mass: " << Sun.mass << std::endl;
    std::cout << "Sun radius: " << Sun.radius << std::endl;
    std::cout << "Sun pos: " << Sun.pos[0] << " " << Sun.pos[1] << " " << Sun.pos[2] << std::endl;
    std::cout << "Sun vel: " << Sun.vel[0] << " " << Sun.vel[1] << " " << Sun.vel[2] << std::endl;

    spiceBody Earth("Earth", 399, t0Sim, 5.972e24L, 6.378137e3L);
    std::cout << "Earth name: " << Earth.name << std::endl;
    std::cout << "Earth spiceID: " << Earth.spiceID << std::endl;
    std::cout << "Earth t0: " << Earth.t0 << std::endl;
    std::cout << "Earth mass: " << Earth.mass << std::endl;
    std::cout << "Earth radius: " << Earth.radius << std::endl;
    std::cout << "Earth pos: " << Earth.pos[0] << " " << Earth.pos[1] << " " << Earth.pos[2] << std::endl;
    std::cout << "Earth vel: " << Earth.vel[0] << " " << Earth.vel[1] << " " << Earth.vel[2] << std::endl;

    std::vector<real> pos = {1.0L, 2.0L, 3.0L};
    std::vector<real> vel = {4.0L, 5.0L, 6.0L};
    std::vector<std::vector<real>> covariance;
    covariance = {  {1.0L, 2.0L, 3.0L, 4.0L, 5.0L, 6.0L},
                    {2.0L, 3.0L, 4.0L, 5.0L, 6.0L, 7.0L},
                    {3.0L, 4.0L, 5.0L, 6.0L, 7.0L, 8.0L},
                    {4.0L, 5.0L, 6.0L, 7.0L, 8.0L, 9.0L},
                    {5.0L, 6.0L, 7.0L, 8.0L, 9.0L, 10.0L},
                    {6.0L, 7.0L, 8.0L, 9.0L, 10.0L, 11.0L}};
    integBody Didymos("65803 Didymos", t0Sim, 5.5e11L, 0.780L, pos, vel, covariance, 0.0L, -1.8e-14L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 1.0L);
    std::cout << "Didymos name: " << Didymos.name << std::endl;
    std::cout << "Didymos t0: " << Didymos.t0 << std::endl;
    std::cout << "Didymos mass: " << Didymos.mass << std::endl;
    std::cout << "Didymos radius: " << Didymos.radius << std::endl;
    std::cout << "Didymos pos: " << Didymos.pos[0] << " " << Didymos.pos[1] << " " << Didymos.pos[2] << std::endl;
    std::cout << "Didymos vel: " << Didymos.vel[0] << " " << Didymos.vel[1] << " " << Didymos.vel[2] << std::endl;
    std::cout << "Didymos covariance row 1: " << Didymos.covariance[0][0] << " " << Didymos.covariance[0][1] << " " << Didymos.covariance[0][2] << " " << Didymos.covariance[0][3] << " " << Didymos.covariance[0][4] << " " << Didymos.covariance[0][5] << std::endl;
    std::cout << "Didymos covariance row 2: " << Didymos.covariance[1][0] << " " << Didymos.covariance[1][1] << " " << Didymos.covariance[1][2] << " " << Didymos.covariance[1][3] << " " << Didymos.covariance[1][4] << " " << Didymos.covariance[1][5] << std::endl;
    std::cout << "Didymos covariance row 3: " << Didymos.covariance[2][0] << " " << Didymos.covariance[2][1] << " " << Didymos.covariance[2][2] << " " << Didymos.covariance[2][3] << " " << Didymos.covariance[2][4] << " " << Didymos.covariance[2][5] << std::endl;
    std::cout << "Didymos covariance row 4: " << Didymos.covariance[3][0] << " " << Didymos.covariance[3][1] << " " << Didymos.covariance[3][2] << " " << Didymos.covariance[3][3] << " " << Didymos.covariance[3][4] << " " << Didymos.covariance[3][5] << std::endl;
    std::cout << "Didymos covariance row 5: " << Didymos.covariance[4][0] << " " << Didymos.covariance[4][1] << " " << Didymos.covariance[4][2] << " " << Didymos.covariance[4][3] << " " << Didymos.covariance[4][4] << " " << Didymos.covariance[4][5] << std::endl;
    std::cout << "Didymos covariance row 6: " << Didymos.covariance[5][0] << " " << Didymos.covariance[5][1] << " " << Didymos.covariance[5][2] << " " << Didymos.covariance[5][3] << " " << Didymos.covariance[5][4] << " " << Didymos.covariance[5][5] << std::endl;
    std::cout << "Didymos a1: " << Didymos.ngParams.a1 << std::endl;
    std::cout << "Didymos a2: " << Didymos.ngParams.a2 << std::endl;
    std::cout << "Didymos a3: " << Didymos.ngParams.a3 << std::endl;

    return 0;
};
