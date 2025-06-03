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
    // limits of distribution are dictated by the time span of 
    // earth_620120_240827.bpc and earth_200101_990827_predict.bpc
    std::uniform_real_distribution<> dis(37684.0004767L, 87942.0008007L);

    std::cout
        << "/////////////////////// PCK map accuracy test ///////////////////////"
        << std::endl
        << std::endl;
    std::cout.precision(8);
    SpiceDouble mjd = dis(gen);
    std::string from = "MOON_ME_DE440_ME421";
    std::string to = "J2000";
    if (from.substr(0,5) == "MOON_" || to.substr(0,5) == "MOON_"){
        std::uniform_real_distribution<> dis(-112816.0L, 288976.0L);
        mjd = dis(gen);
    }
    std::string pck_text = "../../../grss/kernels/pck00011.tpc";
    std::string pck_hist = "../../../grss/kernels/earth_historic.bpc";
    std::string pck_latest = "../../../grss/kernels/earth_latest.bpc";
    std::string pck_predict = "../../../grss/kernels/earth_predict.bpc";
    std::string pck_moon = "../../../grss/kernels/moon_pa_de440.bpc";

    //////////// GRSS ////////////
    PckEphemeris pckEphem;
    pckEphem.histPckPath = pck_hist;
    pckEphem.latestPckPath = pck_latest;
    pckEphem.predictPckPath = pck_predict;
    pckEphem.moonPckPath = pck_moon;
    pckEphem.histPck = pck_init(pck_hist);
    pckEphem.latestPck = pck_init(pck_latest);
    pckEphem.predictPck = pck_init(pck_predict);
    pckEphem.moonPck = pck_init(pck_moon);
    std::vector<std::vector<real>> mapMat(6, std::vector<real>(6, 0.0));
    get_pck_rotMat(from, to, mjd, pckEphem, mapMat);
    pck_free(pckEphem.histPck);
    pck_free(pckEphem.latestPck);
    pck_free(pckEphem.predictPck);
    pck_free(pckEphem.moonPck);

    //////////// SPICE ////////////
    // following variable defined from https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/a_old_versions/moon_de440_220930.tf
    SpiceInt center = 301;
    SpiceInt frclass = 2;
    SpiceInt frclsid = 31008;
    SpiceInt frcode = 31008;
    pipool_c("FRAME_MOON_PA_DE440", 1, &frcode);
    pcpool_c("FRAME_31008_NAME", 1, 14, "MOON_PA_DE440");
    pipool_c("FRAME_31008_CLASS", 1, &frclass);
    pipool_c("FRAME_31008_CLASS_ID", 1, &frclsid);
    pipool_c("FRAME_31008_CENTER", 1, &center);
    center = 301;
    frclass = 4;
    frclsid = 31009;
    frcode = 31009;
    pipool_c("FRAME_MOON_ME_DE440_ME421", 1, &frcode);
    pcpool_c("FRAME_31009_NAME", 1, 20, "MOON_ME_DE440_ME421");
    pipool_c("FRAME_31009_CLASS", 1, &frclass);
    pipool_c("FRAME_31009_CLASS_ID", 1, &frclsid);
    pipool_c("FRAME_31009_CENTER", 1, &center);
    pcpool_c("TKFRAME_31009_SPEC", 1, 6, "ANGLES");
    pcpool_c("TKFRAME_31009_RELATIVE", 1, 14, "MOON_PA_DE440");
    pcpool_c("TKFRAME_31009_UNITS", 1, 11, "ARCSECONDS");
    SpiceDouble angles[3] = {67.8526, 78.6944, 0.2785};
    SpiceInt axes[3] = {3, 2, 1};
    pipool_c("TKFRAME_31009_AXES", 3, axes);
    pdpool_c("TKFRAME_31009_ANGLES", 3, angles);

    // the kernel furninshing order matters, don't mess with it
    furnsh_c(pck_text.c_str());
    furnsh_c(pck_moon.c_str());
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
    // convert SPICE frame so it rotates au, au/day instead of km, km/s
    for (int i = 3; i < 6; i++){
        for (int j = 0; j < 3; j++){
            spiceMat[i][j] *= 86400.0;
        }
    }

    // // print both matrices
    // std::cout.precision(16);
    // std::cout << "Spice Matrix:" << std::endl;
    // for (int i = 0; i < 6; i++){
    //     for (int j = 0; j < 6; j++){
    //         std::cout << std::setw(25) << std::scientific << spiceMat[i][j];
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "Map Matrix:" << std::endl;
    // for (int i = 0; i < 6; i++){
    //     for (int j = 0; j < 6; j++){
    //         std::cout << std::setw(25) << std::scientific << mapMat[i][j];
    //     }
    //     std::cout << std::endl;
    // }

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
    std::cout
        << std::endl
        << "/////////////////////// PCK map accuracy test ///////////////////////"
        << std::endl;
    gettimeofday(&t2, NULL);
    tDiff = (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.0L;
    std::cout << "elapsed time: " << tDiff << " sec" << std::endl;
    return 0;
}
