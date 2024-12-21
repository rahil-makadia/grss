/**
 * @file    timeconvert.h
 * @brief   Header file for time conversion functions.
 * @author  Rahil Makadia <makadia2@illinois.edu>
 *
 * @section     LICENSE
 * Copyright (C) 2022-2025 Rahil Makadia
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, see <https://www.gnu.org/licenses>.
 */

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
