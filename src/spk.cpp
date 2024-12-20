/**
 * @file    spk.cpp
 * @brief   Source file for the SPK ephemeris data.
 * @author  Rahil Makadia <makadia2@illinois.edu>
 * @details This file implements the routines to read in the binary SPK files
 * provided by NAIF and to interpolate the ephemeris data to a requested time.
 * The routines here are based on the small-body ephemeris reading routines
 * implemented by Holman et al. (2023) in the ASSIST library. However, those
 * rountines are extended here to also handle the planetary ephemeris data with
 * the same code rather than separatley using different ephemeris file types.
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

// code here is heavily modified but taken from
// https://github.com/matthewholman/assist/tree/main/src/spk
#include "spk.h"

/**
 * @brief Calculate the Modified Julian Date from the Ephemeris Time.
 * 
 * @param[in] et Ephemeris Time (TDB seconds since J2000).
 * @return double Modified Julian Date (TDB).
 */
static double inline _mjd(double et) { return 51544.5 + et / 86400.0; }

/**
 * @param[in] bsp SpkInfo structure.
 */
void spk_free(SpkInfo* bsp) {
    if (bsp == nullptr) {
        return;
    }

    if (bsp->targets) {
        for (int m = 0; m < bsp->num; m++) {
            free(bsp->targets[m].one);
            free(bsp->targets[m].two);
        }
        free(bsp->targets);
    }
    munmap(bsp->map, bsp->len);
    memset((void *)bsp, 0, sizeof(SpkInfo));
    free(bsp);
    bsp = nullptr;
}

/**
 * @param[in] path Path to the SPK file.
 * @return SpkInfo* Pointer to the SpkInfo structure.
 */
SpkInfo* spk_init(const std::string &path) {
    // For file format information, see
    // https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/daf.html#The%20File%20Record
    // https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html#SPK%20File%20Structure

    // Format for one summary
    struct summary {
        double beg;  // begin epoch, seconds since J2000.0
        double end;  // ending epoch
        int code;     // target code
        int cen;     // centre code (10 = sun)
        int ref;     // reference frame (1 = J2000.0)
        int ver;     // type of ephemeris (2 = chebyshev)
        int one;     // initial array address
        int two;     // final array address
    };

    // File is split into records. We read one record at a time.
    union {
        char buf[SPK_RECORD_LEN];
        struct {
            double next;    // The record number of the next summary record in the file. Zero if this is the final summary record.
            double prev;    // The record number of the previous summary record in the file. Zero if this is the initial summary record.
            double nsum;    // Number of summaries in this record
            summary s[25];  // Summaries (25 is the maximum)
        } summaries;        // Summary record
        struct {
            char locidw[8]; // An identification word
            int nd;         // The number of double precision components in each array summary.
            int ni;         // The number of integer components in each array summary.
            char locifn[60];// The internal name or description of the array file.
            int fward;      // The record number of the initial summary record in the file.
            int bward;      // The record number of the final summary record in the file.
        } file;             // File record
    } record;

    // Try opening file.
    int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        throw std::runtime_error("spk_init: Error opening "+path+".");
    }

    // Read the file record.
    read(fd, &record, SPK_RECORD_LEN);
    // Check if the file is a valid double Precision Array File
    std::string full_file_type = "DAF/SPK";
    if (strncmp(record.file.locidw, full_file_type.c_str(), 7) != 0) {
        close(fd);
        throw std::runtime_error(
            "spk_init: Error parsing "+full_file_type+". Incorrect "
            "header.");
    }

    // Check that the size of a summary record is equal to the size of our
    // struct.
    int nc = 8 * (record.file.nd + (record.file.ni + 1) / 2);
    if (nc != sizeof(summary)) {
        close(fd);
        throw std::runtime_error(
            "spk_init: Error parsing "+full_file_type+". Wrong size of "
            "summary record.");
    }

    // Seek until the first summary record using the file record's fward pointer.
    // Record numbers start from 1 not 0 so we subtract 1 to get to the correct record.
    lseek(fd, (record.file.fward - 1) * SPK_RECORD_LEN, SEEK_SET);
    read(fd, record.buf, SPK_RECORD_LEN);

    // We are at the first summary block, validate
    if ((int64_t)record.buf[8] != 0) {
        close(fd);
        throw std::runtime_error(
            "spk_init: Error parsing "+full_file_type+". Cannot find "
            "summary block.");
    }
    // okay, here we go
    SpkInfo *bsp = (SpkInfo *)calloc(1, sizeof(SpkInfo));
    // Loop over records
    while (true) {
        // Loop over summaries
        for (int b = 0; b < (int)record.summaries.nsum; b++) {
            // get current summary
            summary *sum = &record.summaries.s[b];
            // Index in our arrays for current target
            int m = bsp->num - 1;
            // New target?
            if (bsp->num == 0 || sum->code != bsp->targets[m].code) {
                if (bsp->num == bsp->allocatedNum) {
                    bsp->allocatedNum += SPK_CACHE_ITEM_SIZE;  // increase space in batches of SPK_CACHE_ITEM_SIZE
                    bsp->targets = (SpkTarget *)realloc(
                        bsp->targets, bsp->allocatedNum * sizeof(SpkTarget));
                }
                m++;
                bsp->targets[m].code = sum->code;
                bsp->targets[m].cen = sum->cen;
                bsp->targets[m].beg = _mjd(sum->beg);
                bsp->targets[m].res = _mjd(sum->end) - bsp->targets[m].beg;
                bsp->targets[m].one = (int *)calloc(32768, sizeof(int));
                bsp->targets[m].two = (int *)calloc(32768, sizeof(int));
                bsp->targets[m].ind = 0;
                bsp->num++;
            }
            // add index for target
            bsp->targets[m].one[bsp->targets[m].ind] = sum->one;
            bsp->targets[m].two[bsp->targets[m].ind] = sum->two;
            bsp->targets[m].end = _mjd(sum->end);
            bsp->targets[m].ind++;
        }

        // Location of next record
        long n = (long)record.summaries.next - 1;
        if (n < 0) {
            // this is already the last record.
            break;
        } else {
            // Find and read next record
            lseek(fd, n * SPK_RECORD_LEN, SEEK_SET);
            read(fd, record.buf, SPK_RECORD_LEN);
        }
    }

    // Create a map of SPICE ID to index in the targets array
    std::unordered_map<int, int> spiceIdToIdx;
    for (int m = 0; m < bsp->num; m++) {
        spiceIdToIdx[bsp->targets[m].code] = m;
    }
    bsp->spiceIdToIdx = spiceIdToIdx;

    // Get file size
    struct stat sb;
    if (fstat(fd, &sb) < 0) {
        throw std::runtime_error("spk_init: Error calculating size for "+full_file_type+".");
    }
    bsp->len = sb.st_size;

    // Memory map
    bsp->map = mmap(NULL, bsp->len, PROT_READ, MAP_SHARED, fd, 0);
    if (bsp->map == NULL) {
        // this will leak memory
        throw std::runtime_error("spk_init: Error creating memory map.");
    }
    #if defined(MADV_RANDOM)
        if (madvise(bsp->map, bsp->len, MADV_RANDOM) < 0) {
            // this will leak memory
            throw std::runtime_error("spk_init: Error while calling madvise().");
        }
    #endif
    close(fd);
    return bsp;
}

/**
 * @param[in] bsp SpkInfo structure.
 * @param[in] epoch Epoch to compute the state at (MJD TDB).
 * @param[in] spiceId SPICE ID of the body.
 * @param[out] state State+acceleration of the body at the requested epoch [AU, AU/day, AU/day^2].
 */
void spk_calc(SpkInfo *bsp, double epoch, int spiceId, double *state) {
    if (spiceId == 199) spiceId = 1;
    if (spiceId == 299) spiceId = 2;
    SpkTarget *target = &(bsp->targets[bsp->spiceIdToIdx.at(spiceId)]);
    if (epoch < target->beg || epoch > target->end) {
        throw std::runtime_error(
            "The requested time is outside the coverage provided by "
            "the ephemeris file.");
    }
    for (size_t i = 0; i < 9; i++) {
        state[i] = 0.0;
    }
    if (target->cen == 3) {
        spk_calc(bsp, epoch, target->cen, state);
    }
    // find location of 'directory' describing the data records
    const int n = (int)((epoch - target->beg) / target->res);
    double *val = (double *)bsp->map + target->two[n] - 1;
    // record size and number of coefficients per coordinate
    const int R = (int)val[-1];
    const int P = (R - 2) / 3;  // must be < 32 !!
    // pick out the precise record
    const int b = (int)((epoch - _mjd(val[-3])) / (val[-2] / 86400.0));
    val = (double *)bsp->map + (target->one[n] - 1) + b * R;
    // scale to interpolation units
    const double z = (epoch - _mjd(val[0])) / (val[1] / 86400.0);
    // set up Chebyshev polynomials
    double T[32];
    double S[32];
    double U[32];
    T[0] = 1.0;
    T[1] = z;
    S[0] = 0.0;
    S[1] = 1.0;
    U[0] = 0.0;
    U[1] = 0.0;
    U[2] = 4.0;
    for (int p = 2; p < P; p++) {
        T[p] = 2.0 * z * T[p - 1] - T[p - 2];
        S[p] = 2.0 * z * S[p - 1] + 2.0 * T[p - 1] - S[p - 2];
    }
    for (int p = 3; p < P; p++) {
        U[p] = 2.0 * z * U[p - 1] + 4.0 * S[p - 1] - U[p - 2];
    }
    const double c = 1.0 / val[1]; // derivative scaling factor from SPICE/spke02.f and chbint.f
    for (size_t i = 0; i < 3; i++) {
        const int b = 2 + i * P;
        for (int p = 0; p < P; p++) {
            const double v = val[b + p];
            state[i] += v * T[p] / 149597870.7;
            state[i + 3] += v * S[p] * c / 149597870.7 * 86400.0;
            state[i + 6] += v * U[p] * c * c / 149597870.7 * 86400.0 * 86400.0;
        }
    }
}

/**
 * @param[in] spiceId SPICE ID of the body.
 * @param[in] t0_mjd Epoch to compute the state at (MJD TDB).
 * @param[in] ephem Ephemeris data from the PropSimulation.
 * @param[out] state State+acceleration of the body at the requested epoch [AU, AU/day, AU/day^2].
 * @param[in] writeCache If true, the state will be written to the cache.
 */
void get_spk_state(const int &spiceId, const double &t0_mjd, SpkEphemeris &ephem,
                   double *state, const bool &writeCache) {
    bool smallBody = spiceId > 1000000;
    SpkInfo *infoToUse = smallBody ? ephem.sb : ephem.mb;
    // find what cache index corresponds to the requested SPICE ID
    int cacheIdx = infoToUse->spiceIdToIdx.at(spiceId) + (smallBody ? ephem.mb->num : 0);
    // check if t0_mjd is in the ephem cache
    bool t0SomewhereInCache = false;
    for (size_t i = 0; i < SPK_CACHE_SIZE; i++) {
        if (ephem.cache[i].t == t0_mjd) {
            t0SomewhereInCache = true;
            if (ephem.cache[i].items[cacheIdx].t == t0_mjd &&
                ephem.cache[i].items[cacheIdx].spiceId == spiceId) {
                memcpy(state, ephem.cache[i].items[cacheIdx].state, 9 * sizeof(double));
                return;
            }
        }
    }
    // if not, calculate it from the SPK memory map,
    spk_calc(infoToUse, t0_mjd, spiceId, state);
    if (smallBody) {
        double sunState[9];
        get_spk_state(10, t0_mjd, ephem, sunState);
        for (size_t i = 0; i < 9; i++) {
            state[i] += sunState[i];
        }
    }
    if (writeCache) {
        // and add it to the cache
        if (!t0SomewhereInCache) {
            ephem.nextIdxToWrite = (ephem.nextIdxToWrite + 1) % SPK_CACHE_SIZE;
        }
        const size_t idx = ephem.nextIdxToWrite;
        ephem.cache[idx].t = t0_mjd;
        ephem.cache[idx].items[cacheIdx].t = t0_mjd;
        ephem.cache[idx].items[cacheIdx].spiceId = spiceId;
        memcpy(ephem.cache[idx].items[cacheIdx].state, state, 9 * sizeof(double));
    }
}
