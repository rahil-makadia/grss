#include "grss.h"
#include <sys/time.h>
#include <assert.h>
#include <random>
#include "SpiceUsr.h"

int main(){
    timeval t1, t2;
    gettimeofday(&t1, NULL);
    real tDiff;

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    // limits of distribution are dictated by the time span of sb431-n16s.bsp
    std::uniform_real_distribution<> dis(-94576.0, 234192.0);

    std::cout
        << "/////////////////////// SPK map accuracy test ///////////////////////"
        << std::endl
        << std::endl;
    std::cout.precision(8);
    SpiceDouble mjd = dis(gen);
    std::cout << "MJD: " << mjd << std::endl;
    int kernel = 431;
    std::string kernel_sb, kernel_mb;
    std::vector<int> spiceIds;
    kernel_sb = "../../../grss/kernels/sb431-n16s.bsp";
    kernel_mb = "../../../grss/kernels/de430.bsp";
    spiceIds = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 199, 299, 301, 399,
                // de430 asteroids
                2000451, 2000001, 2000065, 2000511, 2000015, 2000031, 2000052, 2000010,
                2000704, 2000048, 2000003, 2000002, 2000016, 2000087, 2000088, 2000004};
    furnsh_c(kernel_sb.c_str());
    furnsh_c(kernel_mb.c_str());
    SpkInfo* mbInfo431 = spk_init(kernel_mb);
    SpkInfo* sbInfo431 = spk_init(kernel_sb);
    SpkEphemeris eph431;
    eph431.mb = mbInfo431;
    eph431.sb = sbInfo431;
    double pos_error, vel_error;
    pos_error = vel_error = 0;
    for (size_t i = 0; i < spiceIds.size(); i++){
        int spiceId = spiceIds[i];
        // Regular SPICE calls
        SpiceDouble state[6];
        SpiceDouble lt;
        SpiceDouble et = mjd_to_et(mjd);
        spkez_c(spiceId, et, "J2000", "NONE", 0, state, &lt);
        state[0] /= 149597870.7;
        state[1] /= 149597870.7;
        state[2] /= 149597870.7;
        state[3] /= 149597870.7/86400;
        state[4] /= 149597870.7/86400;
        state[5] /= 149597870.7/86400;
        double mapState[9];
        get_spk_state(spiceId, mjd, eph431, mapState);
        // std::cout << "id: " << spiceId << ", ";
        // printf("spkstate: %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f\n", state[0], state[1], state[2], state[3], state[4], state[5]);
        // printf("mapstate: %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f\n", mapState[0], mapState[1], mapState[2], mapState[3], mapState[4], mapState[5]);
        // printf("factors: %0.5e %0.5e %0.5e %0.5e %0.5e %0.5e\n", state[0]-mapState[0], state[1]-mapState[1], state[2]-mapState[2], state[3]-mapState[3], state[4]-mapState[4], state[5]-mapState[5]);
        for (size_t j = 0; j < 3; j++){
            pos_error += fabs((state[j]-mapState[j])/state[j]);
            vel_error += fabs((state[j+3]-mapState[j+3])/state[j+3]);
        }
    }
    kclear_c();
    spk_free(mbInfo431);
    spk_free(sbInfo431);
    std::cout << "DE" << kernel << " errors: " << std::endl;
    std::cout << "Position Cumulative Relative Error: " << pos_error*100 << "%" << std::endl;
    std::cout << "Velocity Cumulative Relative Error: " << vel_error*100 << "%" << std::endl;
    assert(pos_error < 1e-5); // 0.001%
    assert(vel_error < 1e-5); // 0.001%

    kernel = 441;
    kernel_sb = "../../../grss/kernels/sb441-n16s.bsp";
    kernel_mb = "../../../grss/kernels/de440.bsp";
    spiceIds = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 199, 299, 301, 399,
                // de440 asteroids
                2000107, 2000001, 2000065, 2000511, 2000015, 2000031, 2000052, 2000010,
                2000704, 2000007, 2000003, 2000002, 2000016, 2000087, 2000088, 2000004};
    furnsh_c(kernel_sb.c_str());
    furnsh_c(kernel_mb.c_str());
    SpkInfo* mbInfo441 = spk_init(kernel_mb);
    SpkInfo* sbInfo441 = spk_init(kernel_sb);
    SpkEphemeris eph441;
    eph441.mb = mbInfo441;
    eph441.sb = sbInfo441;
    pos_error = vel_error = 0;
    for (size_t i = 0; i < spiceIds.size(); i++){
        int spiceId = spiceIds[i];
        // Regular SPICE calls
        SpiceDouble state[6];
        SpiceDouble lt;
        SpiceDouble et = mjd_to_et(mjd);
        spkez_c(spiceId, et, "J2000", "NONE", 0, state, &lt);
        state[0] /= 149597870.7;
        state[1] /= 149597870.7;
        state[2] /= 149597870.7;
        state[3] /= 149597870.7/86400;
        state[4] /= 149597870.7/86400;
        state[5] /= 149597870.7/86400;
        double mapState[9];
        get_spk_state(spiceId, mjd, eph441, mapState);
        // std::cout << "id: " << spiceId << ", ";
        // printf("spkstate: %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f\n", state[0], state[1], state[2], state[3], state[4], state[5]);
        // printf("mapstate: %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f\n", mapState[0], mapState[1], mapState[2], mapState[3], mapState[4], mapState[5]);
        // printf("factors: %0.5e %0.5e %0.5e %0.5e %0.5e %0.5e\n", state[0]-mapState[0], state[1]-mapState[1], state[2]-mapState[2], state[3]-mapState[3], state[4]-mapState[4], state[5]-mapState[5]);
        for (size_t j = 0; j < 3; j++){
            pos_error += fabs((state[j]-mapState[j])/state[j]);
            vel_error += fabs((state[j+3]-mapState[j+3])/state[j+3]);
        }
    }
    kclear_c();
    spk_free(mbInfo441);
    spk_free(sbInfo441);
    std::cout << "DE" << kernel << " errors: " << std::endl;
    std::cout << "Position Cumulative Relative Error: " << pos_error*100 << "%" << std::endl;
    std::cout << "Velocity Cumulative Relative Error: " << vel_error*100 << "%" << std::endl;
    assert(pos_error < 1e-5); // 0.001%
    assert(vel_error < 1e-5); // 0.001%

    std::cout
        << std::endl
        << "/////////////////////// SPK map accuracy test ///////////////////////"
        << std::endl;
    gettimeofday(&t2, NULL);
    tDiff = (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.0L;
    std::cout << "elapsed time: " << tDiff << " sec" << std::endl;
    return 0;
}
