#include "grss.h"
#include <sys/time.h>
#ifdef PROFILE_YES
#include <gperftools/profiler.h>
#endif
/*
export CPUPROFILE_FREQUENCY=10000 && g++ -DPROFILE_YES -std=c++11 -stdlib=libc++ -g3 apophis.cpp ../../../src/approach.cpp ../../../src/force.cpp ../../../src/gr15.cpp ../../../src/interpolate.cpp ../../../src/simulation.cpp ../../../src/spk.cpp ../../../src/utilities.cpp ../../../extern/cspice/lib/cspice.a -o apophis.out -I ../../../include/ -I ../../../extern/cspice/include/ -lprofiler && ./apophis.out && pprof --web ./apophis.out ./apophis.prof
*/

int main() {
    std::cout.precision(5);
#ifdef PROFILE_YES
    ProfilerStart("apophis.prof");
#endif
    timeval t1, t2;
    gettimeofday(&t1, NULL);
    real tDiff;

    int DEkernel = 441;
    std::string DEkernelPath = "../../../grss/kernels/planets_big16_de" +
        std::to_string(DEkernel) + "_1950_2350.tm";
    if (DEkernel == 0) {
        DEkernelPath = "../../../grss/kernels/planets_big16_de441_1950_2350.tm";
    }
    real t0SimMjd = 2.4621385359989386E+06L - 2400000.5L;
    real tfSimMjd = 2.4625030372426095E+06L - 2400000.5L;
    propSimulation simTest("simTest", t0SimMjd, DEkernel, DEkernelPath);

    std::vector<real> tEval = {};
    bool tEvalUTC = false;
    bool evalApparentState = false;
    bool convergedLightTime = false;
    simTest.set_integration_parameters(
        tfSimMjd, tEval, tEvalUTC, evalApparentState, convergedLightTime,
        std::vector<std::vector<real>>{}, true, 1e-3, 25, 5e-3);
    std::vector<real> pos = {-5.58232604283634858966e-01,
                             8.55200571132647247019e-01,
                             3.03631949052953764578e-01};
    std::vector<real> vel = {-1.38187429130724199339e-02,
                             -6.00401017433447106719e-03,
                             -2.57842571752728202256e-03};
    NongravParamaters ngPrms;
    ngPrms.a1 = 4.999999873689E-13L;
    ngPrms.a2 = -2.901085508711E-14;
    ngPrms.a3 = 0.0L;
    ngPrms.alpha = 1.0L;
    ngPrms.k = 0.0L;
    ngPrms.m = 2.0L;
    ngPrms.n = 0.0L;
    ngPrms.r0_au = 1.0L;
    IntegBody Apophis("(99942) Apophis", simTest.integParams.t0, 0.0L, 0.0L,
                      pos, vel, ngPrms);
    simTest.add_integ_body(Apophis);
    simTest.integrate();

    std::cout
        << "/////////////////////// Apophis comparison ///////////////////////"
        << std::endl
        << std::endl;
    std::vector<real> jplFromAssist = {
        1.75755421201475581228e-02, 1.21966728504766619423e+00,
        4.78396516543700689450e-01, -1.35392833847556674776e-02,
        5.35536634997958596940e-04, -1.50699043360527989245e-05};
    std::vector<real> rm = simTest.xInteg;
    std::cout << "difference between JPL code and GRSS: " << std::endl;
    std::cout << "Position (m): [";
    for (size_t i = 0; i < 3; i++) {
        std::cout << (jplFromAssist[i] - rm[i]) * simTest.consts.du2m << ", ";
    }
    std::cout << "]" << std::endl;

    std::cout << "Velocity (m/s): [";
    for (size_t i = 3; i < 6; i++) {
        std::cout << (jplFromAssist[i] - rm[i]) * simTest.consts.duptu2mps
                  << ", ";
    }
    std::cout << "]" << std::endl;

    real distDiff = sqrt(pow(jplFromAssist[0] - rm[0], 2) +
                         pow(jplFromAssist[1] - rm[1], 2) +
                         pow(jplFromAssist[2] - rm[2], 2)) *
        simTest.consts.du2m;
    std::cout << "Distance (m): " << distDiff << std::endl;
    // make sure the difference is less than 10km
    assert(distDiff < 1.0e4L);

    std::cout
        << std::endl
        << "/////////////////////// Apophis comparison ///////////////////////"
        << std::endl;

#ifdef PROFILE_YES
    ProfilerStop();
#endif
    gettimeofday(&t2, NULL);
    tDiff = (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.0L;
    std::cout << "elapsed time: " << tDiff << " sec" << std::endl;

    return 0;
}
