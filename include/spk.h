// code here is modified from
// https://github.com/matthewholman/assist/tree/main/src/spk
#ifndef SPK_H
#define SPK_H

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include "utilities.h"

struct Ephemeris {
    struct spkInfo *mb;
    struct spkInfo *sb;
};

struct cacheItem {
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
};

struct ephemerisCache {
    double *t;
    double dtSign;
    struct cacheItem *items;
};

struct spkTarget {
    int code;    // Target code
    int cen;     // Centre target
    double beg;  // Begin epoch
    double end;  // End epoch
    double res;  // Epoch step
    int *one;    // Record index
    int *two;    // ... ditto
    int ind;     // Length of index
};

struct spkInfo {
    struct spkTarget *targets;
    int num;           // number of targets
    int allocatedNum;  // space allocated for this many targets
    void *map;         // memory map
    size_t len;        // map length
};

spkInfo *spk_init(const std::string &path);
int spk_free(struct spkInfo *pl);
int spk_calc(struct spkInfo *pl, double epoch, int m, double *x, double *y,
             double *z, double *out_vx, double *out_vy, double *out_vz);
void get_spk_state(const int &spiceID, const double &t0_mjd,
                   const Ephemeris &ephem, double state[6]);
#endif
