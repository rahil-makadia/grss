#include "timeconvert.h"

void jd_to_et(const real jd, real &et) {
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    et = (jd - j2000) * day2sec;
}

real jd_to_et(const real jd) {
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    return (jd - j2000) * day2sec;
}

void jd_to_mjd(const real jd, real &mjd) {
    real offset = 2400000.5;
    mjd = jd - offset;
}

real jd_to_mjd(const real jd) {
    real offset = 2400000.5;
    return jd - offset;
}

void et_to_jd(const real et, real &jd) {
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    jd = (et / day2sec) + j2000;
}

real et_to_jd(const real et) {
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    return (et / day2sec) + j2000;
}

void et_to_mjd(const real et, real &mjd) {
    real offset = 2400000.5;
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    mjd = (et / day2sec) - offset + j2000;
}

real et_to_mjd(const real et) {
    real offset = 2400000.5;
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    return (et / day2sec) - offset + j2000;
}

void mjd_to_jd(const real mjd, real &jd) {
    real offset = 2400000.5;
    jd = mjd + offset;
}

real mjd_to_jd(const real mjd) {
    real offset = 2400000.5;
    real jd = mjd + offset;
    return jd;
}

void mjd_to_et(const real mjd, real &et) {
    real offset = 2400000.5;
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    et = (mjd + offset - j2000) * day2sec;
}

real mjd_to_et(const real mjd) {
    real offset = 2400000.5;
    real j2000 = 2451545.0;
    real day2sec = 86400.0;
    real et = (mjd + offset - j2000) * day2sec;
    return et;
}

void delta_at(const real mjdUtc, real &deltaAt) {
    const real leapDatesMjdUtc[] = {
        36934.0, 37300.0, 37512.0, 37665.0, 38334.0, 38395.0, 38486.0,
        38639.0, 38761.0, 38820.0, 38942.0, 39004.0, 39126.0, 39887.0,
        41317.0, 41499.0, 41683.0, 42048.0, 42413.0, 42778.0, 43144.0,
        43509.0, 43874.0, 44239.0, 44786.0, 45151.0, 45516.0, 46247.0,
        47161.0, 47892.0, 48257.0, 48804.0, 49169.0, 49534.0, 50083.0,
        50630.0, 51179.0, 53736.0, 54832.0, 56109.0, 57204.0, 57754.0};
    const real leapSecs[] = {
        1.4178180, 1.4228180, 1.3728180, 1.8458580, 1.9458580, 3.2401300,
        3.3401300, 3.4401300, 3.5401300, 3.6401300, 3.7401300, 3.8401300,
        4.3131700, 4.2131700, 10,        11,        12,        13,
        14,        15,        16,        17,        18,        19,
        20,        21,        22,        23,        24,        25,
        26,        27,        28,        29,        30,        31,
        32,        33,        34,        35,        36,        37};
    const real leapSecDriftMjdUtc[] = {
        37300.0, 37300.0, 37300.0, 37665.0, 37665.0, 38761.0, 38761.0,
        38761.0, 38761.0, 38761.0, 38761.0, 38761.0, 39126.0, 39126.0};
    const real leapSecDrift[] = {
        0.0012960, 0.0012960, 0.0012960, 0.0011232, 0.0011232,
        0.0012960, 0.0012960, 0.0012960, 0.0012960, 0.0012960,
        0.0012960, 0.0012960, 0.0025920, 0.0025920};
    deltaAt = 0;
    if (mjdUtc < leapDatesMjdUtc[0]) {
        return;
    }
    size_t len = sizeof(leapDatesMjdUtc) / sizeof(leapDatesMjdUtc[0]);
    size_t i;
    for (i = len - 1; i >= 0; i--) {
        if (mjdUtc >= leapDatesMjdUtc[i]) {
            break;
        }
    }
    deltaAt = leapSecs[i];
    if (mjdUtc >= leapDatesMjdUtc[0] && mjdUtc < 41317.0) {
        deltaAt +=
            leapSecDrift[i] * (mjdUtc - leapSecDriftMjdUtc[i]);
    }
}

real delta_et_utc(const real mjdUtc) {
    // start of SPICE leapseconds kernels information (currently naif0012.tls)
    const real ttMinusTai = 32.184;
    const real k = 1.657e-3;
    const real eb = 1.671e-2;
    const real m0 = 6.239996;
    const real mDot = 1.99096871e-7;
    // end of SPICE leapseconds kernel information

    real taiMinusUtc;
    delta_at(mjdUtc, taiMinusUtc);
    real tdbSecPastJ2000Approx =
        (mjdUtc - 51544.5) * 86400 + ttMinusTai + taiMinusUtc;
    // equivalent of spice d_nint routine
    if (tdbSecPastJ2000Approx >= 0) {
        tdbSecPastJ2000Approx = floor(tdbSecPastJ2000Approx + 0.5);
    } else {
        tdbSecPastJ2000Approx = -floor(0.5 - tdbSecPastJ2000Approx);
    }
    // Earth-Moon barycenter eccentric anomaly contribution
    const real M = m0 + mDot * tdbSecPastJ2000Approx;
    const real E = M + eb * sin(M);
    const real tdbMinusTtApprox = k * sin(E);
    return tdbMinusTtApprox + ttMinusTai + taiMinusUtc;
}

real delta_et_tdb(const real mjdTdb) {
    // start of SPICE leapseconds kernels information (currently naif0012.tls)
    const real ttMinusTai = 32.184;
    const real k = 1.657e-3;
    const real eb = 1.671e-2;
    const real m0 = 6.239996;
    const real mDot = 1.99096871e-7;
    // end of SPICE leapseconds kernel information

    real tdbSecPastJ2000 = (mjdTdb - 51544.5) * 86400;
    // equivalent of spice d_nint routine
    if (tdbSecPastJ2000 >= 0) {
        tdbSecPastJ2000 = floor(tdbSecPastJ2000 + 0.5);
    } else {
        tdbSecPastJ2000 = -floor(0.5 - tdbSecPastJ2000);
    }
    // Earth-Moon barycenter eccentric anomaly contribution
    const real M = m0 + mDot * tdbSecPastJ2000;
    const real E = M + eb * sin(M);
    const real tdbMinusTtApprox = k * sin(E);
    const real tai_approx = mjdTdb - (tdbMinusTtApprox + ttMinusTai) / 86400;
    real utcApprox = tai_approx;
    real taiMinusUtc;
    delta_at(utcApprox, taiMinusUtc);
    utcApprox = mjdTdb - (tdbMinusTtApprox + ttMinusTai + taiMinusUtc) / 86400;
    delta_at(utcApprox, taiMinusUtc);
    return tdbMinusTtApprox + ttMinusTai + taiMinusUtc;
}