#include "grss.h"
#include <sys/time.h>
#ifdef PROFILE_YES
#include <gperftools/profiler.h>
#endif
/*
export CPUPROFILE_FREQUENCY=10000 && g++ -DPROFILE_YES -std=c++11 -stdlib=libc++ -g3 didymos.cpp ../../../src/approach.cpp ../../../src/force.cpp ../../../src/gr15.cpp ../../../src/interpolate.cpp ../../../src/simulation.cpp ../../../src/spk.cpp ../../../src/utilities.cpp ../../../extern/cspice/lib/cspice.a -o didymos.out -I ../../../include/ -I ../../../extern/cspice/include/ -lprofiler && ./didymos.out && pprof --web ./didymos.out ./didymos.prof
*/

int main() {
    std::cout.precision(5);
    #ifdef PROFILE_YES
        ProfilerStart("didymos.prof");
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
    real t0SimMjd = 59956.5L;
    real numDays = 3000.0L;
    real tfSimMjd = t0SimMjd + numDays;
    propSimulation simTestForward("simTestForward", t0SimMjd, DEkernel, DEkernelPath);
    propSimulation simTestBackward("simTestBackward", tfSimMjd, DEkernel, DEkernelPath);

    std::vector<real> tEval = {};
    bool tEvalUTC = false;
    bool evalApparentState = false;
    bool convergedLightTime = false;
    simTestForward.set_integration_parameters(tfSimMjd, tEval, tEvalUTC,
                                       evalApparentState, convergedLightTime);
    simTestBackward.set_integration_parameters(t0SimMjd, tEval, tEvalUTC,
                                        evalApparentState, convergedLightTime);

    // handle forward sim
    // Didymos solution 204 from JPL Horizons
    std::vector<real> posf = {-4.739057923859238E-01, 1.095212783390737E+00,
                             5.270678006987460E-01};
    std::vector<real> velf = {-1.653652752454120E-02, -8.564528317342794E-04,
                             6.471815241318629E-04};
    NongravParamaters ngPrms;
    ngPrms.a1 = 0.0L;
    ngPrms.a2 = -1.042309691002E-14;
    ngPrms.a3 = 0.0L;
    ngPrms.alpha = 1.0L;
    ngPrms.k = 0.0L;
    ngPrms.m = 2.0L;
    ngPrms.n = 0.0L;
    ngPrms.r0_au = 1.0L;
    IntegBody Didymos204f("(65803) Didymos204", simTestForward.integParams.t0, 0.0L,
                         390.0L, posf, velf, ngPrms);
    simTestForward.add_integ_body(Didymos204f);
    simTestForward.integrate();

    // handle backward sim
    // Didymos state from forward sim
    std::vector<real> posb = {simTestForward.xInteg[0], simTestForward.xInteg[1],
            simTestForward.xInteg[2]};
    std::vector<real> velb = {simTestForward.xInteg[3], simTestForward.xInteg[4],
            simTestForward.xInteg[5]};
    IntegBody Didymos204b("(65803) Didymos204", simTestBackward.integParams.t0, 0.0L,
                         0.0L, posb, velb, ngPrms);
    simTestBackward.add_integ_body(Didymos204b);
    simTestBackward.integrate();

    std::cout
        << "/////////////////////// Didymos comparison ///////////////////////"
        << std::endl
        << std::endl;
    std::vector<real> sun = {-2.137253647445074E+08, 1.820691741416802E+08,
                             9.084856358728494E+07,  -8.872945056199384E+00,
                             8.222904852571648E-01,  5.360178568968080E-01};
    std::vector<real> ast = {1.1623330140015338E+11, 9.1581845691978887E+10,
                             3.4340792838671610E+10, -2.1324483422419686E+04,
                             2.4444706234341073E+04, 1.2445751473790578E+04};
    std::vector<real> sc = {sun[0] + ast[0], sun[1] + ast[1], sun[2] + ast[2],
                            sun[3] + ast[3], sun[4] + ast[4], sun[5] + ast[5]};
    std::vector<real> rm = simTestForward.xInteg;
    rm[0] *= simTestForward.consts.du2m;
    rm[1] *= simTestForward.consts.du2m;
    rm[2] *= simTestForward.consts.du2m;
    rm[3] *= simTestForward.consts.du2m / 86400;
    rm[4] *= simTestForward.consts.du2m / 86400;
    rm[5] *= simTestForward.consts.du2m / 86400;
    std::cout << "difference between JPL code and GRSS: " << std::endl;
    std::cout << "Position (m): [";
    for (size_t i = 0; i < 3; i++) {
        std::cout << sc[i] - rm[i] << ", ";
    }
    std::cout << "]" << std::endl;

    std::cout << "Velocity (m/s): [";
    for (size_t i = 3; i < 6; i++) {
        std::cout << sc[i] - rm[i] << ", ";
    }
    std::cout << "]" << std::endl;

    real distDiff = sqrt(pow(sc[0] - rm[0], 2) + pow(sc[1] - rm[1], 2) +
                         pow(sc[2] - rm[2], 2));
    std::cout << "Distance (m): " << distDiff << std::endl << std::endl;
    // make sure the difference is less than 100m
    assert(distDiff < 100.0L);

    std::cout << "difference from roundtrip integration: " << std::endl;
    std::cout << "Position (m): [";
    for (size_t i = 0; i < 3; i++) {
        std::cout << (simTestBackward.xInteg[i] - posf[i])*simTestForward.consts.du2m
                  << ", ";
    }
    std::cout << "]" << std::endl;

    std::cout << "Velocity (m/s): [";
    for (size_t i = 3; i < 6; i++) {
        std::cout << (simTestBackward.xInteg[i] - velf[i-3])*simTestForward.consts.du2m / 86400
                  << ", ";
    }
    std::cout << "]" << std::endl;

    real roundTripDiff = sqrt(pow(simTestBackward.xInteg[0] - posf[0], 2) +
                              pow(simTestBackward.xInteg[1] - posf[1], 2) +
                              pow(simTestBackward.xInteg[2] - posf[2], 2))*simTestForward.consts.du2m;
    std::cout << "Distance (m): " << roundTripDiff << std::endl;
    // make sure the difference is less than 1% of distDiff
    assert(roundTripDiff < 0.25L*distDiff);

    std::cout
        << std::endl
        << "/////////////////////// Didymos comparison ///////////////////////"
        << std::endl;

    #ifdef PROFILE_YES
    ProfilerStop();
    #endif
    gettimeofday(&t2, NULL);
    tDiff = (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.0L;
    std::cout << "elapsed time: " << tDiff << " sec" << std::endl;

    return 0;
}
