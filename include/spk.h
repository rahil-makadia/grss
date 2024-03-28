// code here is modified from
// https://github.com/matthewholman/assist/tree/main/src/spk
#ifndef SPK_H
#define SPK_H

#include "utilities.h"

/**
 * @brief Structure to hold a single time,pos,vel,accel record from an SPK file.
 *
 * @param spiceId SPICE ID of the body.
 * @param t Time of the state.
 * @param x X position of the state.
 * @param y Y position of the state.
 * @param z Z position of the state.
 * @param vx X velocity of the state.
 * @param vy Y velocity of the state.
 * @param vz Z velocity of the state.
 * @param ax X acceleration of the state.
 * @param ay Y acceleration of the state.
 * @param az Z acceleration of the state.
 */
struct CacheItem {
    int spiceId = -99999;
    double t;
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double ax;
    double ay;
    double az;
};

/**
 * @brief Size of each item in the cache.
 * It is set to 32, which is a bit of a buffer on the usual number of
 * 27 SpiceBodies in each PropSimulation (Sun+8planets+Moon+Pluto+16Asteriods).
 */
#define SPK_CACHE_ITEM_SIZE 32

/**
 * @brief Structure to hold a cache of SPK data.
 *
 * @param t Time of the cache.
 * @param items Array of CacheItems.
 */
struct EphemerisCache {
    double t;
    CacheItem items[SPK_CACHE_ITEM_SIZE];
};

/**
 * @brief Number of items in the cache.
 * This is the number of spk queries that are remembered.
 */
#define SPK_CACHE_SIZE 16

/**
 * @brief Length of a record in an SPK file.
 */
#define RECORD_LEN 1024

/**
 * @brief Structure to hold the data for a single body in an SPK file.
 *
 * @param code SPICE ID of the body.
 * @param cen Center of the target (usually SSB/Sun).
 * @param beg Starting epoch.
 * @param end End epoch.
 * @param res Epoch span.
 * @param one First record index.
 * @param two Final record index.
 * @param ind Span of the records.
 */
struct SpkTarget {
    int code;
    int cen;
    double beg;
    double end;
    double res;
    int *one;
    int *two;
    int ind;
};

/**
 * @brief Structure to hold the data for a single SPK file.
 *
 * @param targets Array of SpkTarget.
 * @param num Number of targets.
 * @param allocatedNum Number of allocated targets.
 * @param map Memory map of the SPK file.
 * @param len Length of the memory map.
 */
struct DafInfo {
    SpkTarget* targets;
    int num;
    int allocatedNum;
    void *map;
    size_t len;
};

/**
 * @brief Structure to hold all the data for the SPK files in a PropSimulation.
 *
 * @param mb Main body ephemeris data.
 * @param sb Small body ephemeris data.
 * @param nextIdxToWrite Next index to write to in the cache.
 * @param cache Cache of recently queried SPK data.
 */
struct Ephemeris {
    std::string mbPath;
    std::string sbPath;
    DafInfo* mb = nullptr;
    DafInfo* sb = nullptr;
    size_t nextIdxToWrite = -1;
    EphemerisCache cache[SPK_CACHE_SIZE];
};

/**
 * @brief Free the memory allocated for an DafInfo structure.
 */
void daf_free(DafInfo *pl);

/**
 * @brief Initialise a DAF file.
 */
DafInfo* daf_init(const std::string &path, const std::string &type);

/**
 * @brief Initialise an SPK file.
 */
DafInfo* spk_init(const std::string &path);

/**
 * @brief Compute pos, vel, and acc for a given body
 * at a given time using an DafInfo structure.
 */
void spk_calc(DafInfo *pl, double epoch, int spiceId, double *out_x,
             double *out_y, double *out_z, double *out_vx, double *out_vy,
             double *out_vz, double *out_ax, double *out_ay, double *out_az);

/**
 * @brief Top level function to get the state of a body at a given time
 * using the ephemeris data in a PropSimulation.
 */
void get_spk_state(const int &spiceId, const double &t0_mjd, Ephemeris &ephem,
                   double state[9]);

/**
 * @brief Using cspice to get the rotation matrix from one frame to another.
 */
void get_pck_rotMat(const std::string &from, const std::string &to,
                    const real &et, std::vector<std::vector<real>> &rotMat);

#endif
