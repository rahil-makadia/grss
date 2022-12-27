/* COMPILE AND RUN USING
clear && g++ utilities.cpp force.cpp simulation.cpp gr15.cpp test.cpp -o test.out -std=c++11 -O3 -fdiagnostics-color=always -I /opt/homebrew/include -lm /opt/homebrew/Cellar/cspice/67/lib/cspice.a -g3 -Wall && ./test.out
clear && /opt/homebrew/Cellar/gcc/12.2.0/bin/g++-12 utilities.cpp force.cpp simulation.cpp gr15.cpp test.cpp -o test.out -std=c++11 -O3 -fdiagnostics-color=always -I /opt/homebrew/include -lm /opt/homebrew/Cellar/cspice/67/lib/cspice.a -g3 -Wall && ./test.out
*/
#include "simulation.h"


int main() {
    std::cout.precision(22);
    // std::cout << "The size on memory of a real variable is: " << sizeof(real) << std::endl;
    
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

    // furnsh_c ( "../spice/planets_big16_de441_1950_2350.bsp" );
    // real t0SimMjd = 51544.5;
    // SpiceBody Sun("Sun", 10, t0SimMjd, 1.989e30L, 6.957e5L);
    // std::cout << "Sun name: " << Sun.name << std::endl;
    // std::cout << "Sun spiceID: " << Sun.spiceID << std::endl;
    // std::cout << "Sun t0: " << Sun.t0 << std::endl;
    // std::cout << "Sun mass: " << Sun.mass << std::endl;
    // std::cout << "Sun radius: " << Sun.radius << std::endl;
    // std::cout << "Sun pos: " << Sun.pos[0] << " " << Sun.pos[1] << " " << Sun.pos[2] << std::endl;
    // std::cout << "Sun vel: " << Sun.vel[0] << " " << Sun.vel[1] << " " << Sun.vel[2] << std::endl;

    // SpiceBody Earth("Earth", 399, t0SimMjd, 5.972e24L, 6.378137e3L);
    // std::cout << "Earth name: " << Earth.name << std::endl;
    // std::cout << "Earth spiceID: " << Earth.spiceID << std::endl;
    // std::cout << "Earth t0: " << Earth.t0 << std::endl;
    // std::cout << "Earth mass: " << Earth.mass << std::endl;
    // std::cout << "Earth radius: " << Earth.radius << std::endl;
    // std::cout << "Earth pos: " << Earth.pos[0] << " " << Earth.pos[1] << " " << Earth.pos[2] << std::endl;
    // std::cout << "Earth vel: " << Earth.vel[0] << " " << Earth.vel[1] << " " << Earth.vel[2] << std::endl;

    // std::vector<real> pos = {1.0L, 2.0L, 3.0L};
    // std::vector<real> vel = {4.0L, 5.0L, 6.0L};
    // std::vector<std::vector<real>> covariance;
    // covariance = {  {1.0L, 2.0L, 3.0L, 4.0L, 5.0L, 6.0L},
    //                 {2.0L, 3.0L, 4.0L, 5.0L, 6.0L, 7.0L},
    //                 {3.0L, 4.0L, 5.0L, 6.0L, 7.0L, 8.0L},
    //                 {4.0L, 5.0L, 6.0L, 7.0L, 8.0L, 9.0L},
    //                 {5.0L, 6.0L, 7.0L, 8.0L, 9.0L, 10.0L},
    //                 {6.0L, 7.0L, 8.0L, 9.0L, 10.0L, 11.0L}};
    // IntegBody Didymos("65803 Didymos", t0SimMjd, 5.5e11L, 0.780L, pos, vel, covariance, 0.0L, -1.8e-14L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 1.0L);
    // std::cout << "Didymos name: " << Didymos.name << std::endl;
    // std::cout << "Didymos t0: " << Didymos.t0 << std::endl;
    // std::cout << "Didymos mass: " << Didymos.mass << std::endl;
    // std::cout << "Didymos radius: " << Didymos.radius << std::endl;
    // std::cout << "Didymos pos: " << Didymos.pos[0] << " " << Didymos.pos[1] << " " << Didymos.pos[2] << std::endl;
    // std::cout << "Didymos vel: " << Didymos.vel[0] << " " << Didymos.vel[1] << " " << Didymos.vel[2] << std::endl;
    // std::cout << "Didymos covariance row 1: " << Didymos.covariance[0][0] << " " << Didymos.covariance[0][1] << " " << Didymos.covariance[0][2] << " " << Didymos.covariance[0][3] << " " << Didymos.covariance[0][4] << " " << Didymos.covariance[0][5] << std::endl;
    // std::cout << "Didymos covariance row 2: " << Didymos.covariance[1][0] << " " << Didymos.covariance[1][1] << " " << Didymos.covariance[1][2] << " " << Didymos.covariance[1][3] << " " << Didymos.covariance[1][4] << " " << Didymos.covariance[1][5] << std::endl;
    // std::cout << "Didymos covariance row 3: " << Didymos.covariance[2][0] << " " << Didymos.covariance[2][1] << " " << Didymos.covariance[2][2] << " " << Didymos.covariance[2][3] << " " << Didymos.covariance[2][4] << " " << Didymos.covariance[2][5] << std::endl;
    // std::cout << "Didymos covariance row 4: " << Didymos.covariance[3][0] << " " << Didymos.covariance[3][1] << " " << Didymos.covariance[3][2] << " " << Didymos.covariance[3][3] << " " << Didymos.covariance[3][4] << " " << Didymos.covariance[3][5] << std::endl;
    // std::cout << "Didymos covariance row 5: " << Didymos.covariance[4][0] << " " << Didymos.covariance[4][1] << " " << Didymos.covariance[4][2] << " " << Didymos.covariance[4][3] << " " << Didymos.covariance[4][4] << " " << Didymos.covariance[4][5] << std::endl;
    // std::cout << "Didymos covariance row 6: " << Didymos.covariance[5][0] << " " << Didymos.covariance[5][1] << " " << Didymos.covariance[5][2] << " " << Didymos.covariance[5][3] << " " << Didymos.covariance[5][4] << " " << Didymos.covariance[5][5] << std::endl;
    // std::cout << "Didymos a1: " << Didymos.ngParams.a1 << std::endl;
    // std::cout << "Didymos a2: " << Didymos.ngParams.a2 << std::endl;
    // std::cout << "Didymos a3: " << Didymos.ngParams.a3 << std::endl;

    Simulation simTest1("simTest1");
    furnsh_c ( "../spice/planets_big16_de441_1950_2350.bsp" );
    real t0SimJd = 2457380.0L;
    real t0SimMjd;
    jd_to_mjd(t0SimJd, t0SimMjd);
    real num_years = 120.0L;
    real num_days = num_years * 365.25L;
    real tfSimMjd = t0SimMjd + num_days;
    simTest1.set_sim_constants();
    real G = simTest1.consts.G;
    simTest1.set_integration_parameters(t0SimMjd, tfSimMjd, true, 1.0L, 6.0L, 1e-7L, 0.25L);

    // add planets and planetary bodies from DE441 header (https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_tech-comments.txt)
    SpiceBody Sun("Sun", 10, simTest1.integParams.t0, 2.9591220828411956e-04L/G, 6.957e8L, simTest1.consts);
    SpiceBody MercuryBarycenter("Mercury Barycenter", 1, simTest1.integParams.t0, 4.9125001948893182e-11L/G, 0.0L, simTest1.consts);
    SpiceBody VenusBarycenter("Venus Barycenter", 2, simTest1.integParams.t0, 7.2434523326441187e-10L/G, 0.0L, simTest1.consts);
    SpiceBody Earth("Earth", 399, simTest1.integParams.t0, 8.8876924467071022e-10L/G, 6378136.3L, simTest1.consts);
    SpiceBody Moon("Moon", 301, simTest1.integParams.t0, 1.0931894624024351e-11L/G, 0.0L, simTest1.consts);
    SpiceBody MarsBarycenter("Mars Barycenter", 4, simTest1.integParams.t0, 9.5495488297258119e-11L/G, 0.0L, simTest1.consts);
    SpiceBody JupiterBarycenter("Jupiter Barycenter", 5, simTest1.integParams.t0, 2.8253458252257917e-07L/G, 0.0L, simTest1.consts);
    SpiceBody SaturnBarycenter("Saturn Barycenter", 6, simTest1.integParams.t0, 8.4597059933762903e-08L/G, 0.0L, simTest1.consts);
    SpiceBody UranusBarycenter("Uranus Barycenter", 7, simTest1.integParams.t0, 1.2920265649682399e-08L/G, 0.0L, simTest1.consts);
    SpiceBody NeptuneBarycenter("Neptune Barycenter", 8, simTest1.integParams.t0, 1.5243573478851939e-08L/G, 0.0L, simTest1.consts);
    SpiceBody PlutoBarycenter("Pluto Barycenter", 9, simTest1.integParams.t0, 2.1750964648933581e-12L/G, 0.0L, simTest1.consts);
    Sun.isPPN = true;
    Sun.set_J2(2.1106088532726840e-7, 7.25*DEG2RAD);
    Earth.isPPN = true;
    Earth.set_J2(1.0826266835532573e-3, EARTH_OBLIQUITY);
    JupiterBarycenter.isPPN = true;
    simTest1.add_spice_body(Sun);
    simTest1.add_spice_body(MercuryBarycenter);
    simTest1.add_spice_body(VenusBarycenter);
    simTest1.add_spice_body(Earth);
    simTest1.add_spice_body(Moon);
    simTest1.add_spice_body(MarsBarycenter);
    simTest1.add_spice_body(JupiterBarycenter);
    simTest1.add_spice_body(SaturnBarycenter);
    simTest1.add_spice_body(UranusBarycenter);
    simTest1.add_spice_body(NeptuneBarycenter);
    simTest1.add_spice_body(PlutoBarycenter);

    // add asteroids
    // DE441 big16 asteroids from JPL SSD IOM 392R-21-005 (ftp://ssd.jpl.nasa.gov/pub/eph/small_bodies/asteroids_de441/SB441_IOM392R-21-005_perturbers.pdf)
    SpiceBody Ceres("Ceres", 2000001, simTest1.integParams.t0, 1.3964518123081070e-13L/G, 0.0L, simTest1.consts);
    SpiceBody Vesta("Vesta", 2000004, simTest1.integParams.t0, 3.8548000225257904e-14L/G, 0.0L, simTest1.consts);
    SpiceBody Pallas("Pallas", 2000002, simTest1.integParams.t0, 3.0471146330043200e-14L/G, 0.0L, simTest1.consts);
    SpiceBody Hygiea("Hygiea", 2000010, simTest1.integParams.t0, 1.2542530761640810e-14L/G, 0.0L, simTest1.consts);
    SpiceBody Davida("Davida", 2000511, simTest1.integParams.t0, 8.6836253492286545e-15L/G, 0.0L, simTest1.consts);
    SpiceBody Interamnia("Interamnia", 2000704, simTest1.integParams.t0, 6.3110343420878887e-15L/G, 0.0L, simTest1.consts);
    SpiceBody Europa("Europa", 2000052, simTest1.integParams.t0, 5.9824315264869841e-15L/G, 0.0L, simTest1.consts);
    SpiceBody Sylvia("Sylvia", 2000087, simTest1.integParams.t0, 4.8345606546105521e-15L/G, 0.0L, simTest1.consts);
    SpiceBody Eunomia("Eunomia", 2000015, simTest1.integParams.t0, 4.5107799051436795e-15L/G, 0.0L, simTest1.consts);
    SpiceBody Juno("Juno", 2000003, simTest1.integParams.t0, 4.2823439677995011e-15L/G, 0.0L, simTest1.consts);
    SpiceBody Psyche("Psyche", 2000016, simTest1.integParams.t0, 3.5445002842488978e-15L/G, 0.0L, simTest1.consts);
    SpiceBody Camilla("Camilla", 2000107, simTest1.integParams.t0, 3.2191392075878588e-15L/G, 0.0L, simTest1.consts);
    SpiceBody Thisbe("Thisbe", 2000088, simTest1.integParams.t0, 2.6529436610356353e-15L/G, 0.0L, simTest1.consts);
    SpiceBody Iris("Iris", 2000007, simTest1.integParams.t0, 2.5416014973471498e-15L/G, 0.0L, simTest1.consts);
    SpiceBody Euphrosyne("Euphrosyne", 2000031, simTest1.integParams.t0, 2.4067012218937576e-15L/G, 0.0L, simTest1.consts);
    SpiceBody Cybele("Cybele", 2000065, simTest1.integParams.t0, 2.0917175955133682e-15L/G, 0.0L, simTest1.consts);
    simTest1.add_spice_body(Ceres);
    simTest1.add_spice_body(Vesta);
    simTest1.add_spice_body(Pallas);
    simTest1.add_spice_body(Hygiea);
    simTest1.add_spice_body(Davida);
    simTest1.add_spice_body(Interamnia);
    simTest1.add_spice_body(Europa);
    simTest1.add_spice_body(Sylvia);
    simTest1.add_spice_body(Eunomia);
    simTest1.add_spice_body(Juno);
    simTest1.add_spice_body(Psyche);
    simTest1.add_spice_body(Camilla);
    simTest1.add_spice_body(Thisbe);
    simTest1.add_spice_body(Iris);
    simTest1.add_spice_body(Euphrosyne);
    simTest1.add_spice_body(Cybele);

    // add integration bodies
    // Didymos solution 181 from JPL Horizons
    real didyEcc = 3.838828022218969e-01;
    real didyPeriDist = 1.013062336281239e+00;
    real didyIncDeg = 3.407768167104253e+00;
    real didyLongAscNodeDeg = 7.322791476483350e+01;
    real didyArgPeriDeg = 3.192333230140641e+02;
    real didyInc, didyLongAscNode, didyArgPeri;
    deg_to_rad(didyIncDeg, didyInc);
    deg_to_rad(didyLongAscNodeDeg, didyLongAscNode);
    deg_to_rad(didyArgPeriDeg, didyArgPeri);
    real didyT0 = 2.457563407789108e+06-simTest1.consts.JdMinusMjd;

    std::vector<real> didyComState = {didyEcc, didyPeriDist, didyT0, didyInc, didyLongAscNode, didyArgPeri};
    std::vector<std::vector<real>> didyCov =   {{7.583475051964887e-18, -1.328908097977630e-17, -6.972304924721117e-15,  2.283173935512507e-14, -2.686764859997025e-14, -2.616429398334697e-15, -1.581004123520010e-24},
                                                {-1.328908097977630e-17,  2.412510709691607e-17,  1.390386748241651e-14, -3.950950552041731e-14,  4.650581628454678e-14,  4.556738803984023e-15,  7.214046260003573e-24},
                                                {-6.972304924721117e-15,  1.390386748241651e-14,  3.619570291905178e-11, -2.916044151815316e-11,  3.517556751183672e-11,  3.493038361472854e-12, -1.692677626695315e-20},
                                                {2.283173935512507e-14, -3.950950552041731e-14, -2.916044151815316e-11,  9.350629567728970e-11, -1.059971509732549e-10, -1.173201125655876e-11,  8.257884060755273e-21},
                                                {-2.686764859997025e-14,  4.650581628454678e-14,  3.517556751183672e-11, -1.059971509732549e-10,  1.209283974305069e-10,  1.324551477944010e-11, -1.034612590981367e-20},
                                                {-2.616429398334697e-15,  4.556738803984023e-15,  3.493038361472854e-12, -1.173201125655876e-11,  1.324551477944010e-11,  1.762153561997360e-12, -9.405832454921979e-22},
                                                {-1.581004123520010e-24,  7.214046260003573e-24, -1.692677626695315e-20,  8.257884060755273e-21, -1.034612590981367e-20, -9.405832454921979e-22,  5.222314155442253e-29}};
    NongravParams didyNongrav;
    didyNongrav.a1 = 0.0L;
    didyNongrav.a2 = -1.885839515169381e-14;
    didyNongrav.a3 = 0.0L;
    didyNongrav.alpha = 1.0L;
    didyNongrav.k = 0.0L;
    didyNongrav.m = 2.0L;
    didyNongrav.n = 0.0L;
    didyNongrav.r0_au = 1.0L;
    IntegBody Didymos("Didymos", simTest1.integParams.t0, 0.0L, 0.0L, didyComState, didyCov, didyNongrav, simTest1.consts);
    simTest1.add_integ_body(Didymos);

    std::cout<< "Size of simulation: "<< sizeof(simTest1) << " bytes" << std::endl;
    // // check orbital elements and their conversion
    // std::cout<< "Input Cometary state: " << std::endl;
    // for (int i = 0; i < 6; i++){
    //     std::cout<< didyComState[i] << std::endl;
    // }
    // std::vector<real> didyComState0 = didyComState;
    // std::vector<real> didyKepState;
    // std::vector<real> didyCartState;
    // cometary_to_keplerian(simTest1.integParams.t0, didyComState0, didyKepState, simTest1.consts.G);
    // std::cout<< "Keplerian state: " << std::endl;
    // for (int i = 0; i < 6; i++){
    //     std::cout<< didyKepState[i] << std::endl;
    // }
    // keplerian_to_cometary(simTest1.integParams.t0, didyKepState, didyComState0, simTest1.consts.G);
    // std::cout<< "Cometary state: " << std::endl;
    // for (int i = 0; i < 6; i++){
    //     std::cout<< didyComState0[i] << std::endl;
    // }
    // keplerian_to_cartesian(simTest1.integParams.t0, didyKepState, didyCartState, simTest1.consts.G);
    // std::cout<< "Cartesian state: " << std::endl;
    // for (int i = 0; i < 6; i++){
    //     std::cout<< didyCartState[i] << std::endl;
    // }
    // cartesian_to_keplerian(simTest1.integParams.t0, didyCartState, didyKepState, simTest1.consts.G);
    // std::cout<< "Keplerian state: " << std::endl;
    // for (int i = 0; i < 6; i++){
    //     std::cout<< didyKepState[i] << std::endl;
    // }
    // cometary_to_cartesian(simTest1.integParams.t0, didyComState0, didyCartState, simTest1.consts.G);
    // std::cout<< "Cartesian state: " << std::endl;
    // for (int i = 0; i < 6; i++){
    //     std::cout<< didyCartState[i] << std::endl;
    // }
    // cartesian_to_cometary(simTest1.integParams.t0, didyCartState, didyComState0, simTest1.consts.G);
    // std::cout<< "Cometary state: " << std::endl;
    // for (int i = 0; i < 6; i++){
    //     std::cout<< didyComState0[i] << std::endl;
    // }
    // std::cout<< "Input Cometary state: " << std::endl;
    // for (int i = 0; i < 6; i++){
    //     std::cout<< didyComState[i] << std::endl;
    // }

    // // print the covariance matrix
    // std::cout<< "Input covariance matrix: " << std::endl;
    // for (int i = 0; i < 7; i++){
    //     for (int j = 0; j < 7; j++){
    //         std::cout<< didyCov[i][j] << " ";
    //     }
    //     std::cout<< std::endl;
    // }

    real VhEcc = 1.443297491498406E-01;
    real VhPeriDist = 9.732482184639847E-01;
    real VhIncDeg = 1.391152461868084E+01;
    real VhLongAscNodeDeg = 1.393694568920910E+02;
    real VhArgPeriDeg = 2.069227657017103E+02;
    real VhInc, VhLongAscNode, VhArgPeri;
    deg_to_rad(VhIncDeg, VhInc);
    deg_to_rad(VhLongAscNodeDeg, VhLongAscNode);
    deg_to_rad(VhArgPeriDeg, VhArgPeri);
    real VhT0 = 2457375.443007633090-simTest1.consts.JdMinusMjd;

    std::vector<real> VhComState = {VhEcc, VhPeriDist, VhT0, VhInc, VhLongAscNode, VhArgPeri};
    // fill 7x7 covariance matrix with zeros
    std::vector<std::vector<real>> VhCov(7, std::vector<real>(7, 0.0L));
    NongravParams VhNongrav;
    VhNongrav.a1 = 0.0L;
    VhNongrav.a2 = 1.883069273845E-15;
    VhNongrav.a3 = 0.0L;
    VhNongrav.alpha = 1.0L;
    VhNongrav.k = 0.0L;
    VhNongrav.m = 2.0L;
    VhNongrav.n = 0.0L;
    VhNongrav.r0_au = 1.0L;
    IntegBody Vh("1991 VH", simTest1.integParams.t0, 0.0L, 0.0L, VhComState, VhCov, VhNongrav, simTest1.consts);
    simTest1.add_integ_body(Vh);

    simTest1.preprocess();
    // real conv = 1;
    real conv = simTest1.consts.du2m/1000.0L;
    std::cout << "t: " << simTest1.t << std::endl;
    std::cout << "xInteg: " << "np.array([";
    std::cout << simTest1.xInteg[0]*conv << ", ";
    std::cout << simTest1.xInteg[1]*conv << ", ";
    std::cout << simTest1.xInteg[2]*conv << ", ";
    std::cout << simTest1.xInteg[3]*conv << ", ";
    std::cout << simTest1.xInteg[4]*conv << ", ";
    std::cout << simTest1.xInteg[5]*conv << " ";
    std::cout << "])" << std::endl;
    std::cout << "timestepCounter: " << simTest1.integParams.timestepCounter << std::endl;

    simTest1.integrate();
    
    std::cout << "t: " << simTest1.t << std::endl;
    std::cout << "xInteg: " << "np.array([";
    std::cout << simTest1.xInteg[0]*conv << ", ";
    std::cout << simTest1.xInteg[1]*conv << ", ";
    std::cout << simTest1.xInteg[2]*conv << ", ";
    std::cout << simTest1.xInteg[3]*conv << ", ";
    std::cout << simTest1.xInteg[4]*conv << ", ";
    std::cout << simTest1.xInteg[5]*conv << " ";
    std::cout << "])" << std::endl;
    std::cout << "timestepCounter: " << simTest1.integParams.timestepCounter << std::endl;

    simTest1.set_integration_parameters(tfSimMjd, t0SimMjd, true, 1.0L, 6.0L, 1e-7L, 0.25L);

    simTest1.integrate();

    std::cout << "t: " << simTest1.t << std::endl;
    std::cout << "xInteg: " << "np.array([";
    std::cout << simTest1.xInteg[0]*conv << ", ";
    std::cout << simTest1.xInteg[1]*conv << ", ";
    std::cout << simTest1.xInteg[2]*conv << ", ";
    std::cout << simTest1.xInteg[3]*conv << ", ";
    std::cout << simTest1.xInteg[4]*conv << ", ";
    std::cout << simTest1.xInteg[5]*conv << " ";
    std::cout << "])" << std::endl;
    std::cout << "timestepCounter: " << simTest1.integParams.timestepCounter << std::endl;

    std::cout<< "Size of simulation: "<< sizeof(simTest1) << " bytes" << std::endl;
    return 0;
};
