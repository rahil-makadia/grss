#include "grss.h"
#include <sys/time.h>
#include <assert.h>
// #include <gperftools/profiler.h>

int main() {
    std::cout.precision(5);
    // ProfilerStart("test.prof");
    timeval t1, t2;
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
    propSimulation simTest("simTest", t0SimMjd, DEkernel, DEkernelPath);

    std::vector<real> tEval = {};
    bool tEvalUTC = false;
    bool evalApparentState = false;
    bool convergedLightTime = false;
    simTest.set_integration_parameters(tfSimMjd, tEval, tEvalUTC,
                                       evalApparentState, convergedLightTime);

    // add integration bodies
    // Didymos solution 204 from JPL Horizons
    std::vector<real> pos = {-4.739057923859238E-01, 1.095212783390737E+00,
                             5.270678006987460E-01};
    std::vector<real> vel = {-1.653652752454120E-02, -8.564528317342794E-04,
                             6.471815241318629E-04};
    std::vector<std::vector<real> > cov(6, std::vector<real>(6, 0.0L));
    NongravParamaters ngPrms;
    ngPrms.a1 = 0.0L;
    ngPrms.a2 = -1.042309691002E-14;
    ngPrms.a3 = 0.0L;
    ngPrms.alpha = 1.0L;
    ngPrms.k = 0.0L;
    ngPrms.m = 2.0L;
    ngPrms.n = 0.0L;
    ngPrms.r0_au = 1.0L;
    IntegBody Didymos204("(65803) Didymos204", simTest.integParams.t0, 0.0L,
                         0.0L, pos, vel, cov, ngPrms, simTest.consts);
    simTest.add_integ_body(Didymos204);

    // propagate
    gettimeofday(&t1, NULL);

    simTest.preprocess();
    simTest.integrate();

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
    std::vector<real> rm = simTest.xInteg;
    rm[0] *= simTest.consts.du2m;
    rm[1] *= simTest.consts.du2m;
    rm[2] *= simTest.consts.du2m;
    rm[3] *= simTest.consts.du2m / 86400;
    rm[4] *= simTest.consts.du2m / 86400;
    rm[5] *= simTest.consts.du2m / 86400;
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
    std::cout << "Distance (m): " << distDiff << std::endl;
    // make sure the difference is less than 50m
    assert(distDiff < 50.0L);

    std::cout
        << std::endl
        << "/////////////////////// Didymos comparison ///////////////////////"
        << std::endl;

    // ProfilerStop();
    gettimeofday(&t2, NULL);
    tDiff = (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.0L;
    std::cout << "elapsed time: " << tDiff << " sec" << std::endl;

    return 0;
}
