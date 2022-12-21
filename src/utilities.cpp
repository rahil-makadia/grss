#include "utilities.h"

void get_spice_state_lt(int spiceID, real t0_mjd, Constants consts, double state[6], double &lt){
    real t0_et = mjd_to_et(t0_mjd);
    SpiceInt center = 0; // 0 = solar system barycenter
    ConstSpiceChar *frame="J2000"; // Earth mean equator and equinox of J2000, states output will be ICRF-EME2000 frame
    ConstSpiceChar *abcorr="NONE"; // No aberration correction
    spkez_c(spiceID, t0_et, frame, abcorr, center, state, &lt);
    for (int i=0; i<6; i++){
        state[i] /= consts.du2m;
    }
    for (int i=3; i<6; i++){
        state[i] *= consts.tu2sec;
    }
};

void jd_to_et(const real jd, real &et){
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    et = (jd-j2000)*day2sec;
};

real jd_to_et(const real jd){
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    return (jd-j2000)*day2sec;
};

void jd_to_mjd(const real jd, real &mjd){
    real offset = 2400000.5;
    mjd = jd-offset;
};

real jd_to_mjd(const real jd){
    real offset = 2400000.5;
    return jd-offset;
};

void et_to_jd(const real et, real &jd){
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    jd = (et/day2sec)+j2000;
};

real et_to_jd(const real et){
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    return (et/day2sec)+j2000;
};

void et_to_mjd(const real et, real &mjd){
    real offset = 2400000.5;
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    mjd = (et/day2sec)-offset+j2000;
};

real et_to_mjd(const real et){
    real offset = 2400000.5;
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    return (et/day2sec)-offset+j2000;
};

void mjd_to_jd(const real mjd, real &jd){
    real offset = 2400000.5;
    jd = mjd+offset;
};

real mjd_to_jd(const real mjd){
    real offset = 2400000.5;
    real jd = mjd+offset;
    return jd;
};

void mjd_to_et(const real mjd, real &et){
    real offset = 2400000.5;
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    et = (mjd+offset-j2000)*day2sec;
};

real mjd_to_et(const real mjd){
    real offset = 2400000.5;
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    real et = (mjd+offset-j2000)*day2sec;
    return et;
};
