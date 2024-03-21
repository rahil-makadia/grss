#ifndef TIME_H
#define TIME_H

#include "utilities.h"

void jd_to_et(const real jd, real &et);
real jd_to_et(const real jd);
void jd_to_mjd(const real jd, real &mjd);
real jd_to_mjd(const real jd);
void et_to_jd(const real et, real &jd);
real et_to_jd(const real et);
void et_to_mjd(const real et, real &mjd);
real et_to_mjd(const real et);
void mjd_to_jd(const real mjd, real &jd);
real mjd_to_jd(const real mjd);
void mjd_to_et(const real mjd, real &et);
real mjd_to_et(const real mjd);

real delta_at_utc(const real mjdUtc);
real delta_at_tai(const real mjdTai);
real delta_et_utc(const real mjdUtc);
real delta_et_tdb(const real mjdTdb);
#endif
