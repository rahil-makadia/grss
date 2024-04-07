#include "grss.h"
#include <sys/time.h>
#include <assert.h>
#include <random>

int main(){
    timeval t1, t2;
    gettimeofday(&t1, NULL);
    real tDiff;

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    // limits of distribution are dictated by the time span of 
    // earth_720101_230601.bpc and earth_200101_990825_predict.bpc
    std::uniform_real_distribution<> dis(41317.0005, 87940.0008);

    std::cout
        << "/////////////////////// PCK map accuracy test ///////////////////////"
        << std::endl
        << std::endl;
    std::cout.precision(8);
    SpiceDouble mjd = dis(gen);
    std::string from = "ITRF93";
    std::string to = "J2000";
    std::string pck_text = "../../../grss/kernels/pck00011.tpc";
    std::string pck_hist = "../../../grss/kernels/earth_720101_230601.bpc";
    std::string pck_latest = "../../../grss/kernels/earth_latest_high_prec.bpc";
    std::string pck_predict = "../../../grss/kernels/earth_200101_990825_predict.bpc";

    //////////// GRSS ////////////
    PckEphemeris pckEphem;
    pckEphem.histPckPath = pck_hist;
    pckEphem.latestPckPath = pck_latest;
    pckEphem.predictPckPath = pck_predict;
    pckEphem.histPck = pck_init(pck_hist);
    pckEphem.latestPck = pck_init(pck_latest);
    pckEphem.predictPck = pck_init(pck_predict);
    std::vector<std::vector<real>> mapMat(6, std::vector<real>(6, 0.0));
    get_pck_rotMat(from, to, mjd, pckEphem, mapMat);
    pck_free(pckEphem.histPck);
    pck_free(pckEphem.latestPck);
    pck_free(pckEphem.predictPck);

    //////////// SPICE ////////////
    // the kernel furninshing order matters, don't mess with it
    furnsh_c(pck_text.c_str());
    furnsh_c(pck_predict.c_str());
    furnsh_c(pck_latest.c_str());
    furnsh_c(pck_hist.c_str());
    SpiceDouble et = mjd_to_et(mjd);
    std::cout << "MJD: " << mjd << std::endl;
    std::cout << "From: " << from << std::endl;
    std::cout << "To: " << to << std::endl;
    SpiceDouble spiceMat[6][6];
    sxform_c(from.c_str(), to.c_str(), et, spiceMat);
    kclear_c();

    // compute relative error
    real error = 0;
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){
            if (spiceMat[i][j] != 0){
                error += fabs((spiceMat[i][j] - mapMat[i][j])/spiceMat[i][j]);
            }
        }
    }
    std::cout << "Cumulative Relative Error: " << error*100 << "%" << std::endl;
    assert(error < 1e-5); // 0.001%
    // // print both matrices
    // std::cout << "Spice Matrix:" << std::endl;
    // for (int i = 0; i < 6; i++){
    //     for (int j = 0; j < 6; j++){
    //         std::cout << std::setw(20) << std::scientific << spiceMat[i][j];
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "Map Matrix:" << std::endl;
    // for (int i = 0; i < 6; i++){
    //     for (int j = 0; j < 6; j++){
    //         std::cout << std::setw(20) << std::scientific << mapMat[i][j];
    //     }
    //     std::cout << std::endl;
    // }
    std::cout
        << std::endl
        << "/////////////////////// PCK map accuracy test ///////////////////////"
        << std::endl;
    gettimeofday(&t2, NULL);
    tDiff = (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.0L;
    std::cout << "elapsed time: " << tDiff << " sec" << std::endl;
    return 0;
}
