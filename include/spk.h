// from https://github.com/matthewholman/assist/tree/main/src/spk
#ifndef SPK_H
#define SPK_H

#include <cstring>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "utilities.h"

struct mpos_s {
    double u[3];
    double v[3];
    double w[3];
    double jde;
};

struct spk_target {
    int code;     // Target code
    int cen;      // Centre target
    double mass;  // Mass. Set to 0 if not found in ephemeris file.
    double beg;   // Begin epoch
    double end;   // End epoch
    double res;   // Epoch step
    int *one;     // Record index
    int *two;     // ... ditto
    int ind;      // Length of index
};

struct spk_s {
    struct spk_target *targets;
    int num;            // number of targets
    int allocated_num;  // space allocated for this many targets
    void *map;          // memory map
    size_t len;         // map length
};

spk_s *assist_spk_init(const char *path);
int assist_spk_free(struct spk_s *pl);
int assist_spk_calc(struct spk_s *pl, double jde, double rel, int m, double *GM,
                    double *x, double *y, double *z, double *out_vx,
                    double *out_vy, double *out_vz);

// Order of columns in JPL Ephemeris file
enum {
    JPL_MER,  // Mercury
    JPL_VEN,  // Venus
    JPL_EMB,  // Earth
    JPL_MAR,  // Mars
    JPL_JUP,  // Jupiter
    JPL_SAT,  // Saturn
    JPL_URA,  // Uranus
    JPL_NEP,  // Neptune
    JPL_PLU,  // Pluto
    JPL_LUN,  // Moon (geocentric)
    JPL_SUN,  // the Sun
    JPL_NUT,  // nutations
    JPL_LIB,  // lunar librations
    JPL_MAN,  // lunar mantle
    JPL_TDB,  // TT-TDB (< 2 ms)

    JPL_N,  // Number of columns
};

struct jpl_s {
    double beg, end;     // begin and end times
    double inc;          // time step size
    double cau;          // definition of AU
    double cem;          // Earth/Moon mass ratio
    int32_t num;         // number of constants
    int32_t ver;         // ephemeris version
    int32_t off[JPL_N];  // indexing offset
    int32_t ncf[JPL_N];  // number of chebyshev coefficients
    int32_t niv[JPL_N];  // number of interpolation intervals
    int32_t ncm[JPL_N];  // number of components / dimension
    double mass[JPL_N];  // G*mass for all bodies
    double J2E;          // Other constant names follow JPL
    double J3E;
    double J4E;
    double J2SUN;
    double AU;
    double RE;
    double CLIGHT;
    double ASUN;
    size_t len, rec;  // file and record sizes
    void *map;        // memory mapped location
    double *con;      // constant values
    char **str;       // constant names
};

jpl_s *assist_jpl_init(const char *path);
void assist_jpl_free(struct jpl_s *jpl);
void assist_jpl_work(double *P, int ncm, int ncf, int niv, double t0, double t1,
                     double *u, double *v, double *w);
int assist_jpl_calc(struct jpl_s *pl, double jd_ref, double jd_rel, int body,
                    double *const GM, double *const x, double *const y,
                    double *const z, double *const vx, double *const vy,
                    double *const vz, double *const ax, double *const ay,
                    double *const az);

// From Weryk's code
/////// private interface :

static inline void vecpos_off(double *u, const double *v, const double w) {
    u[0] += v[0] * w;
    u[1] += v[1] * w;
    u[2] += v[2] * w;
}
static inline void vecpos_set(double *u, const double *v) {
    u[0] = v[0];
    u[1] = v[1];
    u[2] = v[2];
}
static inline void vecpos_nul(double *u) { u[0] = u[1] = u[2] = 0.0; }
static inline void vecpos_div(double *u, double v) {
    u[0] /= v;
    u[1] /= v;
    u[2] /= v;
}

#endif
