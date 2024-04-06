#include "grss.h"
#include <sys/time.h>
#include <assert.h>

int main(){
    timeval t1, t2;
    gettimeofday(&t1, NULL);
    real tDiff;

    std::cout
        << "/////////////////////// PCK map accuracy test ///////////////////////"
        << std::endl
        << std::endl;
    std::cout.precision(8);
    SpiceDouble mjd = 51544.5;
    std::string from = "ITRF93";
    std::string to = "ECLIPJ2000";
    std::string pck_kernel = "../../../grss/kernels/earth_720101_230601.bpc";

    furnsh_c(pck_kernel.c_str());
    SpiceDouble et = mjd_to_et(mjd);
    std::cout << "MJD: " << mjd << std::endl;
    std::cout << "From: " << from << std::endl;
    std::cout << "To: " << to << std::endl;
    SpiceDouble spiceMat[6][6];
    sxform_c(from.c_str(), to.c_str(), et, spiceMat);
    unload_c(pck_kernel.c_str());

    PckInfo* bpc = pck_init(pck_kernel);
    std::vector<std::vector<real>> mapMat(6, std::vector<real>(6, 0.0));
    get_pck_rotMat2(from, to, mjd, bpc, mapMat);
    pck_free(bpc);

    // compute relative error
    real error = 0;
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){
            if (spiceMat[i][j] != 0){
                error += fabs((spiceMat[i][j] - mapMat[i][j])/spiceMat[i][j]);
            }
        }
    }
    // // print both matrices
    // std::cout << "Spice Matrix:" << std::endl;
    // for (int i = 0; i < 6; i++){
    //     for (int j = 0; j < 6; j++){
    //         std::cout << std::setw(20) << spiceMat[i][j];
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "Map Matrix:" << std::endl;
    // for (int i = 0; i < 6; i++){
    //     for (int j = 0; j < 6; j++){
    //         std::cout << std::setw(20) << mapMat[i][j];
    //     }
    //     std::cout << std::endl;
    // }
    std::cout << "Cumulative Relative Error: " << error*100 << "%" << std::endl;
    assert(error < 1e-5); // 0.001%
    std::cout
        << std::endl
        << "/////////////////////// PCK map accuracy test ///////////////////////"
        << std::endl;
    gettimeofday(&t2, NULL);
    tDiff = (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.0L;
    std::cout << "elapsed time: " << tDiff << " sec" << std::endl;
    return 0;
}
