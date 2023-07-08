#include "spk.h"

int assist_spk_free(spk_s *pl) {
    if (pl == NULL) {
        return -1;
    }
    if (pl->targets) {
        for (int m = 0; m < pl->num; m++) {
            free(pl->targets[m].one);
            free(pl->targets[m].two);
        }
        free(pl->targets);
    }
    munmap(pl->map, pl->len);
    memset(pl, 0, sizeof(spk_s));
    free(pl);
    return 0;
}

// convert SPK epoch (since J2000.0) to julian day number
static double inline _jul(double eph) { return 2451545.0 + eph / 86400.0; }

#define record_length 1024
// check for any non-7bit ascii characters
static int _com(const char *record) {
    for (int n = 0; n < record_length; n++) {
        if (record[n] < 0) return 0;
    }

    return 1;
}

spk_s *assist_spk_init(const char *path) {
    // For file format information, see
    // https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/daf.html

    // Format for one summary
    struct sum_s {
        double beg;  // begin epoch, seconds since J2000.0
        double end;  // ending epoch
        int tar;     // target code
        int cen;     // centre code (10 = sun)
        int ref;     // reference frame (1 = J2000.0)
        int ver;     // type of ephemeris (2 = chebyshev)
        int one;     // initial array address
        int two;     // final array address
    };

    // File is split into records. We read one record at a time.
    union {
        char buf[record_length];
        struct {
            double next;  // The record number of the next summary record in the
                          // file. Zero if this is the final summary record.
            double
                prev;  // The record number of the previous summary record in
                       // the file. Zero if this is the initial summary record.
            double nsum;  // Number of summaries in this record
            sum_s s[25];  // Summaries (25 is the maximum)
        } summary;        // Summary record
        struct {
            char locidw[8];  // An identification word
            int nd;  // The number of double precision components in each array
                     // summary.
            int ni;  // The number of integer components in each array summary.
        } file;      // File record
    } record;

    // Try opening file.
    int fd = open(path, O_RDONLY);
    if (fd < 0) {
        return NULL;
    }

    // Read the file record.
    read(fd, &record, 1024);
    // Check if the file is a valid Double Precision Array File
    if (strncmp(record.file.locidw, "DAF/SPK", 7) != 0) {
        fprintf(stderr, "Error parsing DAF/SPK file. Incorrect header.\n");
        close(fd);
        return NULL;
    }

    // Check that the size of a summary record is equal to the size of our
    // struct.
    int nc = 8 * (record.file.nd + (record.file.ni + 1) / 2);
    if (nc != sizeof(sum_s)) {
        fprintf(stderr,
                "Error parsing DAF/SPK file. Wrong size of summary record.\n");
        close(fd);
        return NULL;
    }

    // Continue reading file until we find a non-ascii record.
    do {
        read(fd, record.buf, 1024);
    } while (_com(record.buf) > 0);

    // We are at the first summary block, validate
    if ((int64_t)record.buf[8] != 0) {
        fprintf(stderr,
                "Error parsing DAF/SPL file. Cannot find summary block.\n");
        close(fd);
        return NULL;
    }
    // std::cout << "record.summary.nsum: " << record.summary.nsum << std::endl;
    // std::cout << "record.file.nd: " << record.file.nd << std::endl;
    // std::cout << "record.file.ni: " << record.file.ni << std::endl;
    // std::cout << "nc: " << nc << std::endl;
    // okay, let's go
    spk_s *pl = (spk_s *)calloc(1, sizeof(spk_s));
    while (1) {  // Loop over records
        for (int b = 0; b < (int)record.summary.nsum;
             b++) {                             // Loop over summaries
            sum_s *sum = &record.summary.s[b];  // get current summary
            // Index in our arrays for current target
            int m = pl->num - 1;
            // New target?
            if (pl->num == 0 || sum->tar != pl->targets[m].code) {
                if (pl->num <= pl->allocated_num) {
                    pl->allocated_num += 32;  // increase space in batches of 32
                    pl->targets = (spk_target *)realloc(
                        pl->targets, pl->allocated_num * sizeof(spk_target));
                }
                m++;
                pl->targets[m].code = sum->tar;
                pl->targets[m].cen = sum->cen;
                pl->targets[m].beg = _jul(sum->beg);
                pl->targets[m].res = _jul(sum->end) - pl->targets[m].beg;
                pl->targets[m].one = (int *)calloc(32768, sizeof(int));
                pl->targets[m].two = (int *)calloc(32768, sizeof(int));
                pl->targets[m].ind = 0;
                pl->num++;
            }
            // add index for target
            pl->targets[m].one[pl->targets[m].ind] = sum->one;
            pl->targets[m].two[pl->targets[m].ind] = sum->two;
            pl->targets[m].end = _jul(sum->end);
            pl->targets[m].ind++;
        }

        // Location of next record
        long n = (long)record.summary.next - 1;
        if (n < 0) {
            // this is already the last record.
            break;
        } else {
            // Find and read next record
            lseek(fd, n * 1024, SEEK_SET);
            read(fd, record.buf, 1024);
        }
    }

    // Get file size
    struct stat sb;
    if (fstat(fd, &sb) < 0) {
        fprintf(stderr, "Error calculating size for DAF/SPL file.\n");
        return NULL;
    }
    pl->len = sb.st_size;

    // Memory map
    pl->map = mmap(NULL, pl->len, PROT_READ, MAP_SHARED, fd, 0);
    if (pl->map == NULL) {
        fprintf(stderr, "Error creating memory map.\n");
        return NULL;  // Will leak memory
    }
#if defined(MADV_RANDOM)
    if (madvise(pl->map, pl->len, MADV_RANDOM) < 0) {
        fprintf(stderr, "Error while calling madvise().\n");
        return NULL;  // Will leak memory
    }
#endif
    close(fd);
    return pl;
}

int assist_spk_calc(spk_s *pl, double jde, double rel, int m, double *GM,
                    double *out_x, double *out_y, double *out_z, double *out_vx,
                    double *out_vy, double *out_vz) {
    if (pl == NULL) {
        throw std::runtime_error("The JPL ephemeris file has not been found.");
    }

    if (m < 0 || m >= pl->num) {
        throw std::runtime_error("The requested spice ID has not been found.");
    }
    spk_target *target = &(pl->targets[m]);

    if (jde + rel < target->beg || jde + rel > target->end) {
        throw std::runtime_error(
            "The requested time is outside the coverage "
            "provided by the ephemeris file.");
    }

    *GM = target->mass;  // Note mass constants defined in DE440/441 ephemeris
                         // files. If not found mass of 0 is used.

    double T[32];
    double S[32];
    double *val, z;
    mpos_s pos;

    for (size_t n = 0; n < 3; n++) {
        pos.u[n] = pos.v[n] = 0.0;
    }

    int n, b, p, P, R;
    // find location of 'directory' describing the data records
    n = (int)((jde + rel - target->beg) / target->res);
    val = (double *)pl->map + target->two[n] - 1;

    // record size and number of coefficients per coordinate
    R = (int)val[-1];
    P = (R - 2) / 3;  // must be < 32 !!

    // pick out the precise record
    b = (int)(((jde - _jul(val[-3])) + rel) / (val[-2] / 86400.0));
    //+ sizeof(double) * b * R;
    val = (double *)pl->map + (target->one[n] - 1) + b * R;

    // scale to interpolation units
    z = ((jde - _jul(val[0])) + rel) / (val[1] / 86400.0);
    // std::cout << "n: " << n << std::endl;
    // std::cout << "b: " << b << std::endl;
    // std::cout << "p: " << p << std::endl;
    // std::cout << "P: " << P << std::endl;
    // std::cout << "R: " << R << std::endl;
    // std::cout << "z: " << z << std::endl;

    // set up Chebyshev polynomials
    T[0] = 1.0;
    T[1] = z;
    // Not used at the moment:
    S[0] = 0.0;
    S[1] = 1.0;

    for (p = 2; p < P; p++) {
        T[p] = 2.0 * z * T[p - 1] - T[p - 2];
        // Not used at the moment:
        S[p] = 2.0 * z * S[p - 1] + 2.0 * T[p - 1] - S[p - 2];
    }
    double c = (0.125 / 2 / 86400.0);
    for (n = 0; n < 3; n++) {
        b = 2 + n * P;
        // std::cout << "b: " << b << std::endl;
        // sum interpolation stuff
        for (p = 0; p < P; p++) {
            pos.u[n] += val[b + p] * T[p];
            // Not used at the moment:
            pos.v[n] += val[b + p] * S[p] * c;
        }
        // restore units to [AU] and [AU/day]
        // pos.u[n] /= 149597870.7;
        // Not used at the moment:
        // pos.v[n] /= 149597870.7 / 86400.0;
        // pos.v[n] /= val[1];
    }

    *out_x = pos.u[0];
    *out_y = pos.u[1];
    *out_z = pos.u[2];

    *out_vx = pos.v[0];
    *out_vy = pos.v[1];
    *out_vz = pos.v[2];

    return 0;
}

////////////////////////////////////////// Planet
/////////////////////////////////////////////

void assist_jpl_work(double *P, int ncm, int ncf, int niv, double t0, double t1,
                     double *u, double *v, double *w) {
    double T[24], S[24];
    double U[24];
    double t, c;
    int p, m, n, b;

    // adjust to correct interval
    t = t0 * (double)niv;
    t0 = 2.0 * fmod(t, 1.0) - 1.0;
    c = (double)(niv * 2) / t1 / 86400.0;

    b = (int)t;
    // std::cout << "t: " << t << std::endl;
    // std::cout << "t0: " << t0 << std::endl;
    // std::cout << "t1: " << t1 << std::endl;
    // std::cout << "c: " << c << std::endl;
    // std::cout << "b: " << b << std::endl;
    // std::cout << "niv: " << niv << std::endl;
    // std::cout << "ncm: " << ncm << std::endl;
    // std::cout << "ncf: " << ncf << std::endl;

    // set up Chebyshev polynomials and derivatives
    T[0] = 1.0;
    T[1] = t0;
    S[0] = 0.0;
    S[1] = 1.0;
    U[0] = 0.0;
    U[1] = 0.0;
    U[2] = 4.0;

    for (p = 2; p < ncf; p++) {
        T[p] = 2.0 * t0 * T[p - 1] - T[p - 2];
        S[p] = 2.0 * t0 * S[p - 1] + 2.0 * T[p - 1] - S[p - 2];
    }
    for (p = 3; p < ncf; p++) {
        U[p] = 2.0 * t0 * U[p - 1] + 4.0 * S[p - 1] - U[p - 2];
    }

    // compute the position/velocity
    for (m = 0; m < ncm; m++) {
        u[m] = v[m] = w[m] = 0.0;
        n = ncf * (m + b * ncm);
        // std::cout << "n: " << n << std::endl;
        for (p = 0; p < ncf; p++) {
            u[m] += T[p] * P[n + p];
            v[m] += S[p] * P[n + p] * c;
            w[m] += U[p] * P[n + p] * c * c;
        }
    }
}

static double getConstant(jpl_s *jpl, const char *name) {
    for (int p = 0; p < jpl->num; p++) {
        if (strncmp(name, jpl->str[p], 6) == 0) {
            return jpl->con[p];
        }
    }
    fprintf(stderr, "WARNING: Constant [%s] not found in ephemeris file.\n",
            name);
    return 0;
}

jpl_s *assist_jpl_init(const char *path) {
    struct stat sb;
    ssize_t ret;
    int fd;

    if ((fd = open(path, O_RDONLY)) < 0) {
        return NULL;
    }

    if (fstat(fd, &sb) < 0) {
        close(fd);
        fprintf(stderr, "Error while trying to determine filesize.\n");
        return NULL;
    }

    // skip the header and constant names for now
    if (lseek(fd, 0x0A5C, SEEK_SET) < 0) {
        close(fd);
        fprintf(stderr, "Error while seeking to header.\n");
        return NULL;
    }
    jpl_s *jpl = (jpl_s *)calloc(1, sizeof(jpl_s));
    // read header
    ret = read(fd, &jpl->beg, sizeof(double));    // Start JD
    ret += read(fd, &jpl->end, sizeof(double));   // End JD
    ret += read(fd, &jpl->inc, sizeof(double));   // Days per block
    ret += read(fd, &jpl->num, sizeof(int32_t));  // Number of constants
    ret += read(fd, &jpl->cau, sizeof(double));   // AU to km
    ret += read(fd, &jpl->cem, sizeof(double));   // Ratio between Earth/Moon

    // number of coefficients for all components
    for (int p = 0; p < JPL_N; p++) {
        jpl->ncm[p] = 3;
    }
    // exceptions:
    jpl->ncm[JPL_NUT] = 2;  // nutations
    jpl->ncm[JPL_TDB] = 1;  // TT-TDB

    for (int p = 0; p < 12; p++) {  // Columns 1-12 of Group 1050
        ret += read(fd, &jpl->off[p], sizeof(int32_t));
        ret += read(fd, &jpl->ncf[p], sizeof(int32_t));
        ret += read(fd, &jpl->niv[p], sizeof(int32_t));
    }

    ret += read(fd, &jpl->ver, sizeof(int32_t));  // Version. e.g. 440
    ret +=
        read(fd, &jpl->off[12], sizeof(int32_t));  // Columns 13 of Group 1050
    ret += read(fd, &jpl->ncf[12], sizeof(int32_t));
    ret += read(fd, &jpl->niv[12], sizeof(int32_t));

    // Get all the constant names
    jpl->str = (char **)calloc(jpl->num, sizeof(char *));

    // retrieve the names of the first 400 constants
    lseek(fd, 0x00FC, SEEK_SET);
    for (int p = 0; p < 400; p++) {  // Group 1040
        jpl->str[p] = (char *)calloc(8, sizeof(char));
        read(fd, jpl->str[p], 6);
    }

    // read the remaining constant names
    lseek(fd, 0x0B28, SEEK_SET);
    for (int p = 400; p < jpl->num; p++) {
        jpl->str[p] = (char *)calloc(8, sizeof(char));
        read(fd, jpl->str[p], 6);
    }

    for (int p = 13; p < 15; p++) {  // Columns 14 and 15 of Group 1050
        ret += read(fd, &jpl->off[p], sizeof(int32_t));
        ret += read(fd, &jpl->ncf[p], sizeof(int32_t));
        ret += read(fd, &jpl->niv[p], sizeof(int32_t));
    }

    // adjust for correct indexing (ie: zero based)
    for (int p = 0; p < JPL_N; p++) {
        jpl->off[p] -= 1;
    }

    // save file size, and determine 'kernel size' or 'block size' (=8144 bytes
    // for DE440/441)
    jpl->len = sb.st_size;
    jpl->rec = sizeof(double) * 2;

    for (int p = 0; p < JPL_N; p++) {
        jpl->rec += sizeof(double) * jpl->ncf[p] * jpl->niv[p] * jpl->ncm[p];
    }

    // memory map the file, which makes us thread-safe with kernel caching
    jpl->map = mmap(NULL, jpl->len, PROT_READ, MAP_SHARED, fd, 0);

    if (jpl->map == NULL) {
        close(fd);
        free(jpl);  // note constants leak
        fprintf(stderr, "Error while calling mmap().\n");
        return NULL;
    }

    // Read constants
    jpl->con = (double *)calloc(jpl->num, sizeof(double));
    lseek(fd, jpl->rec, SEEK_SET);  // Starts at offset of 1 block size
    for (int p = 0; p < jpl->num; p++) {
        read(fd, &jpl->con[p], sizeof(double));
        // printf("%6d  %s   %.5e\n",p,jpl->str[p],jpl->con[p]);
    }

    // Find masses
    jpl->mass[0] = getConstant(jpl, "GMS   ");    // Sun
    jpl->mass[1] = getConstant(jpl, "GM1   ");    // Mercury
    jpl->mass[2] = getConstant(jpl, "GM2   ");    // Venus
    double emrat = getConstant(jpl, "EMRAT  ");   // Earth Moon Ratio
    double gmb = getConstant(jpl, "GMB   ");      // Earth Moon combined
    jpl->mass[3] = (emrat / (1. + emrat)) * gmb;  // Earth
    jpl->mass[4] = 1. / (1 + emrat) * gmb;        // Moon
    jpl->mass[5] = getConstant(jpl, "GM4   ");    // Mars
    jpl->mass[6] = getConstant(jpl, "GM5   ");    // Jupiter
    jpl->mass[7] = getConstant(jpl, "GM6   ");    // Saturn
    jpl->mass[8] = getConstant(jpl, "GM7   ");    // Uranus
    jpl->mass[9] = getConstant(jpl, "GM8   ");    // Neptune
    jpl->mass[10] = getConstant(jpl, "GM9   ");   // Pluto

    // Other constants
    jpl->J2E = getConstant(jpl, "J2E   ");
    jpl->J3E = getConstant(jpl, "J3E   ");
    jpl->J4E = getConstant(jpl, "J4E   ");
    jpl->J2SUN = getConstant(jpl, "J2SUN ");
    jpl->AU = getConstant(jpl, "AU    ");
    jpl->RE = getConstant(jpl, "RE    ");
    jpl->CLIGHT = getConstant(jpl, "CLIGHT");
    jpl->ASUN = getConstant(jpl, "ASUN  ");

    // this file descriptor is no longer needed since we are memory mapped
    if (close(fd) < 0) {
        fprintf(stderr, "Error while closing file.\n");
    }
#if defined(MADV_RANDOM)
    if (madvise(jpl->map, jpl->len, MADV_RANDOM) < 0) {
        fprintf(stderr, "Error during madvise.\n");
    }
#endif
    return jpl;
}

void assist_jpl_free(jpl_s *jpl) {
    if (jpl == NULL) {
        return;
    }
    if (munmap(jpl->map, jpl->len) < 0) {
        fprintf(stderr, "Error during munmap().\n");
    }
    for (int p = 0; p < jpl->num; p++) {
        free(jpl->str[p]);
    }
    free(jpl->str);
    free(jpl->con);
    free(jpl);
}

/*
 *  jpl_calc
 *
 *  Caculate the position+velocity in _equatorial_ coordinates.
 *  Assumes pos is initially zero.
 */
int assist_jpl_calc(jpl_s *jpl, double jd_ref, double jd_rel, int body,
                    double *const GM, double *const out_x, double *const out_y,
                    double *const out_z, double *const out_vx,
                    double *const out_vy, double *const out_vz,
                    double *const out_ax, double *const out_ay,
                    double *const out_az) {
    double t, *z;
    uint32_t blk;

    if (jpl == NULL || jpl->map == NULL) {
        throw std::runtime_error("The JPL ephemeris file has not been found.");
    }
    // if body is less than 0 or not one of
    // 0,1,199,2,299,3,301,399,4,5,6,7,8,9,10, throw error
    if (body < 0 ||
        (body != 0 && body != 1 && body != 199 && body != 2 && body != 299 &&
         body != 3 && body != 301 && body != 399 && body != 4 && body != 5 &&
         body != 6 && body != 7 && body != 8 && body != 9 && body != 10)) {
        throw std::runtime_error("The requested spice ID has not been found.");
    }
    mpos_s pos;

    // check if covered by this file
    if (jd_ref + jd_rel < jpl->beg || jd_ref + jd_rel > jpl->end) {
        throw std::runtime_error(
            "The requested time is outside the coverage "
            "provided by the ephemeris file.");
    }

    // compute record number and 'offset' into record
    blk = (uint32_t)((jd_ref + jd_rel - jpl->beg) / jpl->inc);
    z = (double *)jpl->map + (blk + 2) * jpl->rec / sizeof(double);
    t = ((jd_ref - jpl->beg - (double)blk * jpl->inc) + jd_rel) / jpl->inc;

    // Get mass, position, velocity, and mass of body i in barycentric coords.
    switch (body) {  // The indices in the jpl-> arrays match the JPL component
                     // index for the body
        // solar system barycenter
        case 0:
            // assist_jpl_work(&z[jpl->off[0]], jpl->ncm[0], jpl->ncf[0],
            // jpl->niv[0], t, jpl->inc, pos.u, pos.v, pos.w);
            *GM = 0.0;
            break;
        // sun
        case 10:
            assist_jpl_work(&z[jpl->off[10]], jpl->ncm[10], jpl->ncf[10],
                            jpl->niv[10], t, jpl->inc, pos.u, pos.v, pos.w);
            *GM = jpl->mass[0];
            break;
        // mercury planet center or barycenter
        case 1:
        case 199:
            assist_jpl_work(&z[jpl->off[JPL_MER]], jpl->ncm[JPL_MER],
                            jpl->ncf[JPL_MER], jpl->niv[JPL_MER], t, jpl->inc,
                            pos.u, pos.v, pos.w);
            *GM = jpl->mass[1];
            break;
        // venus planet center or barycenter
        case 2:
        case 299:
            assist_jpl_work(&z[jpl->off[JPL_VEN]], jpl->ncm[JPL_VEN],
                            jpl->ncf[JPL_VEN], jpl->niv[JPL_VEN], t, jpl->inc,
                            pos.u, pos.v, pos.w);
            *GM = jpl->mass[2];
            break;
        // earth-moon barycenter
        case 3:
            assist_jpl_work(&z[jpl->off[JPL_EMB]], jpl->ncm[JPL_EMB],
                            jpl->ncf[JPL_EMB], jpl->niv[JPL_EMB], t, jpl->inc,
                            pos.u, pos.v, pos.w);
            *GM = jpl->mass[3] + jpl->mass[4];
            break;
        // earth planet center
        case 399: {
            mpos_s emb, lun;
            assist_jpl_work(&z[jpl->off[JPL_EMB]], jpl->ncm[JPL_EMB],
                            jpl->ncf[JPL_EMB], jpl->niv[JPL_EMB], t, jpl->inc,
                            emb.u, emb.v, emb.w);  // earth moon barycenter
            assist_jpl_work(&z[jpl->off[JPL_LUN]], jpl->ncm[JPL_LUN],
                            jpl->ncf[JPL_LUN], jpl->niv[JPL_LUN], t, jpl->inc,
                            lun.u, lun.v, lun.w);

            vecpos_set(pos.u, emb.u);
            vecpos_off(pos.u, lun.u, -1.0 / (1.0 + jpl->cem));

            vecpos_set(pos.v, emb.v);
            vecpos_off(pos.v, lun.v, -1.0 / (1.0 + jpl->cem));

            vecpos_set(pos.w, emb.w);
            vecpos_off(pos.w, lun.w, -1.0 / (1.0 + jpl->cem));
            *GM = jpl->mass[3];
        } break;
        // moon body center
        case 301: {
            mpos_s emb, lun;
            assist_jpl_work(&z[jpl->off[JPL_EMB]], jpl->ncm[JPL_EMB],
                            jpl->ncf[JPL_EMB], jpl->niv[JPL_EMB], t, jpl->inc,
                            emb.u, emb.v, emb.w);
            assist_jpl_work(&z[jpl->off[JPL_LUN]], jpl->ncm[JPL_LUN],
                            jpl->ncf[JPL_LUN], jpl->niv[JPL_LUN], t, jpl->inc,
                            lun.u, lun.v, lun.w);

            vecpos_set(pos.u, emb.u);
            vecpos_off(pos.u, lun.u, jpl->cem / (1.0 + jpl->cem));

            vecpos_set(pos.v, emb.v);
            vecpos_off(pos.v, lun.v, jpl->cem / (1.0 + jpl->cem));

            vecpos_set(pos.w, emb.w);
            vecpos_off(pos.w, lun.w, jpl->cem / (1.0 + jpl->cem));
            *GM = jpl->mass[4];
        } break;
        // mars barycenter
        case 4:
            assist_jpl_work(&z[jpl->off[JPL_MAR]], jpl->ncm[JPL_MAR],
                            jpl->ncf[JPL_MAR], jpl->niv[JPL_MAR], t, jpl->inc,
                            pos.u, pos.v, pos.w);
            *GM = jpl->mass[5];
            break;
        // jupiter barycenter
        case 5:
            assist_jpl_work(&z[jpl->off[JPL_JUP]], jpl->ncm[JPL_JUP],
                            jpl->ncf[JPL_JUP], jpl->niv[JPL_JUP], t, jpl->inc,
                            pos.u, pos.v, pos.w);
            *GM = jpl->mass[6];
            break;
        // saturn barycenter
        case 6:
            assist_jpl_work(&z[jpl->off[JPL_SAT]], jpl->ncm[JPL_SAT],
                            jpl->ncf[JPL_SAT], jpl->niv[JPL_SAT], t, jpl->inc,
                            pos.u, pos.v, pos.w);
            *GM = jpl->mass[7];
            break;
        // uranus barycenter
        case 7:
            assist_jpl_work(&z[jpl->off[JPL_URA]], jpl->ncm[JPL_URA],
                            jpl->ncf[JPL_URA], jpl->niv[JPL_URA], t, jpl->inc,
                            pos.u, pos.v, pos.w);
            *GM = jpl->mass[8];
            break;
        // neptune barycenter
        case 8:
            assist_jpl_work(&z[jpl->off[JPL_NEP]], jpl->ncm[JPL_NEP],
                            jpl->ncf[JPL_NEP], jpl->niv[JPL_NEP], t, jpl->inc,
                            pos.u, pos.v, pos.w);
            *GM = jpl->mass[9];
            break;
        // pluto barycenter
        case 9:
            assist_jpl_work(&z[jpl->off[JPL_PLU]], jpl->ncm[JPL_PLU],
                            jpl->ncf[JPL_PLU], jpl->niv[JPL_PLU], t, jpl->inc,
                            pos.u, pos.v, pos.w);
            *GM = jpl->mass[10];
            break;
        default:
            throw std::runtime_error(
                "The requested spice ID has not been found.");  // body not
                                                                // found
            break;
    }

    // Convert to au/day and au/day^2
    // vecpos_div(pos.u, jpl->cau);
    // vecpos_div(pos.v, jpl->cau/86400.);
    // vecpos_div(pos.w, jpl->cau/(86400.*86400.));

    *out_x = pos.u[0];
    *out_y = pos.u[1];
    *out_z = pos.u[2];
    *out_vx = pos.v[0];
    *out_vy = pos.v[1];
    *out_vz = pos.v[2];
    *out_ax = pos.w[0];
    *out_ay = pos.w[1];
    *out_az = pos.w[2];

    return 0;
}
