#ifndef PCK_H
#define PCK_H

#include "timeconvert.h"

/**
 * @brief Length of a record in an PCK file.
 */
#define RECORD_LEN 1024

/**
 * @brief Structure to hold the data for a single frame in an PCK file.
 *
 * @param code SPICE ID of the frame.
 * @param ref Inertial reference of the frame (usually ECLIPJ2000).
 * @param beg Starting epoch.
 * @param end End epoch.
 * @param res Epoch span.
 * @param one First record index.
 * @param two Final record index.
 * @param ind Span of the records.
 */
struct PckTarget {
    int code;
    int ref;
    double beg;
    double end;
    double res;
    int *one;
    int *two;
    int ind;
};

/**
 * @brief Structure to hold the data for a single PCK file.
 *
 * @param targets Array of PckTarget.
 * @param num Number of targets.
 * @param allocatedNum Number of allocated targets.
 * @param map Memory map of the PCK file.
 * @param len Length of the memory map.
 */
struct PckInfo {
    PckTarget* targets;
    int num;
    int allocatedNum;
    void *map;
    size_t len;
};

/**
 * @brief Free the memory allocated for an PckInfo structure.
 */
void pck_free(PckInfo *bpc);

/**
 * @brief Initialise a PCK file.
 */
PckInfo* pck_init(const std::string &path);

/**
 * @brief Structure to hold all the data for the PCK files in a PropSimulation.
 *
 * @param histPckPath Path to the historical PCK file.
 * @param latestPckPath Path to the latest PCK file.
 * @param predictPckPath Path to the predicted PCK file.
 * @param histPck PckInfo structure for the historical PCK file.
 * @param latestPck PckInfo structure for the latest PCK file.
 * @param predictPck PckInfo structure for the predicted PCK file.
 */
struct PckEphemeris {
    std::string histPckPath;
    std::string latestPckPath;
    std::string predictPckPath;
    PckInfo* histPck = nullptr;
    PckInfo* latestPck = nullptr;
    PckInfo* predictPck = nullptr;
};

/**
 * @brief Compute angle,angleDot for a given frame
 * at a given time using an PckInfo structure.
 */
void pck_calc(PckInfo *bpc, real epoch, int spiceId, real *rotMat,
              real *rotMatDot);

/**
 * @brief Compute angle,angleDot for a given frame
 * at a given time using the IAU pole polynomials.
 */
void iau_to_euler(const real t0_mjd, std::string iauFrame, real *euler);

/**
 * @brief Convert a 313 Euler angle to a rotation matrix and its derivative.
 */
void euler313_to_rotMat(const real euler[6], real *rotMat, real *rotMatDot);

/**
 * @brief Get the rotation matrix from one frame to another.
 */
void get_pck_rotMat(const std::string &from, const std::string &to,
                    const real &t0_mjd, PckEphemeris &ephem,
                    std::vector<std::vector<real>> &xformMat);

#endif
