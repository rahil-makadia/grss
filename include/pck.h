#ifndef PCK_H
#define PCK_H

#include "utilities.h"

/**
 * @brief Structure to hold a single time,angle,angleDot record from an PCK file.
 *
 * @param spiceId SPICE ID of the frame.
 * @param t Time of the frame angles.
 * @param alpha Right ascension of the pole.
 * @param delta Declination of the pole.
 * @param w Prime meridian angle.
 * @param alphaDot Right ascension rate.
 * @param deltaDot Declination rate.
 * @param wDot Prime meridian angle rate.
 */
struct PckCacheItem {
    int spiceId = -99999;
    double t;
    double alpha;
    double delta;
    double w;
    double alphaDot;
    double deltaDot;
    double wDot;
};

/**
 * @brief Size of each item in the cache.
 * It is set to 32, which is a bit of a buffer on the usual number of
 * 1 frame in a binary PCK (e.g. ITRF93 for Earth).
 */
#define PCK_CACHE_ITEM_SIZE 4

/**
 * @brief Structure to hold a cache of PCK data.
 *
 * @param t Time of the cache.
 * @param items Array of CacheItems.
 */
struct PckCache {
    double t;
    PckCacheItem items[PCK_CACHE_ITEM_SIZE];
};

/**
 * @brief Number of items in the cache.
 * This is the number of PCK queries that are remembered.
 */
#define PCK_CACHE_SIZE 5

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
void pck_calc(PckInfo *bpc, real epoch, int spiceId, bool inertialJ2000,
              real *rotMat, real *rotMatDot);

void euler313_to_rotmat(const real euler[6], real *rotMat, real *rotMatDot);

/**
 * @brief Using cspice to get the rotation matrix from one frame to another.
 */
void get_pck_rotMat(const std::string &from, const std::string &to,
                    const real &et, std::vector<std::vector<real>> &rotMat);

/**
 * @brief Get the rotation matrix from one frame to another.
 */
void get_pck_rotMat2(const std::string &from, const std::string &to,
                    const real &et, PckInfo* bpc,// PckEphemeris &ephem,
                    std::vector<std::vector<real>> &xformMat);

#endif
