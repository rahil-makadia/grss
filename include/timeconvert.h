#ifndef TIME_H
#define TIME_H

#include "utilities.h"

/**
 * @brief Convert Julian Date to TDB ephemeris time.
 */
void jd_to_et(const real jd, real &et);

/**
 * @brief Convert Julian Date to TDB ephemeris time.
 */
real jd_to_et(const real jd);

/**
 * @brief Convert Julian Date to Modified Julian Date.
 */
void jd_to_mjd(const real jd, real &mjd);

/**
 * @brief Convert Julian Date to Modified Julian Date.
 */
real jd_to_mjd(const real jd);

/**
 * @brief Convert TDB ephemeris time to Julian Date.
 */
void et_to_jd(const real et, real &jd);

/**
 * @brief Convert TDB ephemeris time to Julian Date.
 */
real et_to_jd(const real et);

/**
 * @brief Convert TDB ephemeris time to Modified Julian Date.
 */
void et_to_mjd(const real et, real &mjd);

/**
 * @brief Convert TDB ephemeris time to Modified Julian Date.
 */
real et_to_mjd(const real et);

/**
 * @brief Convert Modified Julian Date to Julian Date.
 */
void mjd_to_jd(const real mjd, real &jd);

/**
 * @brief Convert Modified Julian Date to Julian Date.
 */
real mjd_to_jd(const real mjd);

/**
 * @brief Convert Modified Julian Date to TDB ephemeris time.
 */
void mjd_to_et(const real mjd, real &et);

/**
 * @brief Convert Modified Julian Date to TDB ephemeris time.
 */
real mjd_to_et(const real mjd);

/**
 * @brief Compute TAI-UTC difference in seconds at a given UTC Modified Julian Date.
 */
real delta_at_utc(const real mjdUtc);

/**
 * @brief Compute TAI-UTC difference in seconds at a given TAI Modified Julian Date.
 */
real delta_at_tai(const real mjdTai);

/**
 * @brief Compute TDB-UTC difference in seconds at a given UTC Modified Julian Date.
 */
real delta_et_utc(const real mjdUtc);

/**
 * @brief Compute TDB-UTC difference in seconds at a given TDB Modified Julian Date.
 */
real delta_et_tdb(const real mjdTdb);

#endif
