// code here is modified from
// https://github.com/matthewholman/assist/tree/main/src/spk
#ifndef SPK_H
#define SPK_H

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>
#include <iostream>
#include <vector>
#include <stdexcept>

struct CacheItem {
    int spiceID = -99999;
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

#define SPK_CACHE_ITEM_SIZE 32
struct EphemerisCache {
    double t;
    struct CacheItem items[SPK_CACHE_ITEM_SIZE];
};
#define SPK_CACHE_SIZE 16
#define RECORD_LEN 1024
struct Ephemeris {
    struct SpkInfo *mb;
    struct SpkInfo *sb;
    size_t nextIdxToWrite = -1;
    EphemerisCache cache[SPK_CACHE_SIZE];
};

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
struct SpkInfo {
    struct SpkTarget *targets;
    int num;
    int allocatedNum;
    void *map;
    size_t len;
};
SpkInfo *spk_init(const std::string &path);
void spk_calc(SpkInfo *pl, double epoch, int spiceId, double *out_x,
             double *out_y, double *out_z, double *out_vx, double *out_vy,
             double *out_vz, double *out_ax, double *out_ay, double *out_az);
void get_spk_state(const int &spiceID, const double &t0_mjd, Ephemeris &ephem,
                   double state[9]);

#endif
