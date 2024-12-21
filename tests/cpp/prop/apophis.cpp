/**
 * Apophis propagation unit test inspired by the ASSIST library
 * See https://github.com/matthewholman/assist/tree/main/unit_tests/apophis
 */

#include "grss.h"
#include <sys/time.h>
#include <assert.h>
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

    int DEkernel = 440;
    std::string DEkernelPath = "../../../grss/kernels/";
    const real t0SimMjd = 2462130.5L - 2400000.5L;
    const real tfSimMjd = 2462495.5L - 2400000.5L;
    PropSimulation simTest("simTest", t0SimMjd, DEkernel, DEkernelPath);

    std::vector<real> tEval = {};
    bool tEvalUTC = false;
    bool evalApparentState = false;
    bool convergedLightTime = false;
    simTest.set_integration_parameters(
        tfSimMjd, tEval, tEvalUTC, evalApparentState, convergedLightTime);
    const std::vector<real> pos = {
        -0.4430312741660833, 0.8965062452933791, 0.3218771965080118
    };
    const std::vector<real> vel = {
        -0.0148212842112486, -0.0042500259738234, -0.0019519585329488
    };
    NongravParameters ngPrms;
    ngPrms.a1 = 5.E-13L;
    ngPrms.a2 = -2.901085583204654E-14L;
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
    const std::vector<real> df = {
            0.1194303822347949, 1.2110683438150556, 0.4767184664960469,
           -0.0134714830959050, 0.0017461826766498, 0.0004605752588853
    };
    std::vector<real> rm = simTest.xInteg;
    std::cout << "difference between JPL code and GRSS: " << std::endl;
    std::cout << "Position (m): [";
    for (size_t i = 0; i < 3; i++) {
        std::cout << (df[i] - rm[i]) * simTest.consts.du2m << ", ";
    }
    std::cout << "]" << std::endl;

    std::cout << "Velocity (m/s): [";
    for (size_t i = 3; i < 6; i++) {
        std::cout << (df[i] - rm[i]) * simTest.consts.duptu2mps
                  << ", ";
    }
    std::cout << "]" << std::endl;

    real distDiff = sqrt(pow(df[0] - rm[0], 2) +
                         pow(df[1] - rm[1], 2) +
                         pow(df[2] - rm[2], 2)) *
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
