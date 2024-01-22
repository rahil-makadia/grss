#include "spk.h"

// convert SPK epoch (since J2000.0) to julian day number
static double inline _mjd(double et) { return 51544.5 + et / 86400.0; }

#define RECORD_LEN 1024
// check for any non-7bit ascii characters
static int _com(const char *record) {
    for (int n = 0; n < RECORD_LEN; n++) {
        if (record[n] < 0) return 0;
    }

    return 1;
}

spkInfo *spk_init(const std::string &path) {
    // For file format information, see
    // https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/daf.html

    // Format for one summary
    struct summary {
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
        char buf[RECORD_LEN];
        struct {
            double next;  // The record number of the next summary record in the
                          // file. Zero if this is the final summary record.
            double
                prev;  // The record number of the previous summary record in
                       // the file. Zero if this is the initial summary record.
            double nsum;    // Number of summaries in this record
            summary s[25];  // Summaries (25 is the maximum)
        } summaries;          // Summary record
        struct {
            char locidw[8];  // An identification word
            int nd;  // The number of double precision components in each array
                     // summary.
            int ni;  // The number of integer components in each array summary.
        } file;      // File record
    } record;

    // Try opening file.
    int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        return NULL;
    }

    // Read the file record.
    read(fd, &record, RECORD_LEN);
    // Check if the file is a valid double Precision Array File
    if (strncmp(record.file.locidw, "DAF/SPK", 7) != 0) {
        throw std::runtime_error(
            "Error parsing DAF/SPK file. Incorrect "
            "header.");
        close(fd);
        return NULL;
    }

    // Check that the size of a summary record is equal to the size of our
    // struct.
    int nc = 8 * (record.file.nd + (record.file.ni + 1) / 2);
    if (nc != sizeof(summary)) {
        throw std::runtime_error(
            "Error parsing DAF/SPK file. Wrong size of "
            "summary record.");
        close(fd);
        return NULL;
    }

    // Continue reading file until we find a non-ascii record.
    do {
        read(fd, record.buf, RECORD_LEN);
    } while (_com(record.buf) > 0);

    // We are at the first summary block, validate
    if ((int64_t)record.buf[8] != 0) {
        throw std::runtime_error(
            "Error parsing DAF/SPL file. Cannot find "
            "summary block.");
        close(fd);
        return NULL;
    }
    // std::cout << "record.summary.nsum: " << record.summary.nsum << std::endl;
    // std::cout << "record.file.nd: " << record.file.nd << std::endl;
    // std::cout << "record.file.ni: " << record.file.ni << std::endl;
    // std::cout << "nc: " << nc << std::endl;
    // okay, let's go
    spkInfo *pl = (spkInfo *)calloc(1, sizeof(spkInfo));
    while (1) {  // Loop over records
        for (int b = 0; b < (int)record.summaries.nsum;
             b++) {                               // Loop over summaries
            summary *sum = &record.summaries.s[b];  // get current summary
            // Index in our arrays for current target
            int m = pl->num - 1;
            // New target?
            if (pl->num == 0 || sum->tar != pl->targets[m].code) {
                if (pl->num <= pl->allocatedNum) {
                    pl->allocatedNum += 32;  // increase space in batches of 32
                    pl->targets = (spkTarget *)realloc(
                        pl->targets, pl->allocatedNum * sizeof(spkTarget));
                }
                m++;
                pl->targets[m].code = sum->tar;
                pl->targets[m].cen = sum->cen;
                pl->targets[m].beg = _mjd(sum->beg);
                pl->targets[m].res = _mjd(sum->end) - pl->targets[m].beg;
                pl->targets[m].one = (int *)calloc(32768, sizeof(int));
                pl->targets[m].two = (int *)calloc(32768, sizeof(int));
                pl->targets[m].ind = 0;
                pl->num++;
            }
            // add index for target
            pl->targets[m].one[pl->targets[m].ind] = sum->one;
            pl->targets[m].two[pl->targets[m].ind] = sum->two;
            pl->targets[m].end = _mjd(sum->end);
            pl->targets[m].ind++;
        }

        // Location of next record
        long n = (long)record.summaries.next - 1;
        if (n < 0) {
            // this is already the last record.
            break;
        } else {
            // Find and read next record
            lseek(fd, n * RECORD_LEN, SEEK_SET);
            read(fd, record.buf, RECORD_LEN);
        }
    }

    // Get file size
    struct stat sb;
    if (fstat(fd, &sb) < 0) {
        throw std::runtime_error("Error calculating size for DAF/SPL file.");
        return NULL;
    }
    pl->len = sb.st_size;

    // Memory map
    pl->map = mmap(NULL, pl->len, PROT_READ, MAP_SHARED, fd, 0);
    if (pl->map == NULL) {
        throw std::runtime_error("Error creating memory map.");
        return NULL;  // Will leak memory
    }
#if defined(MADV_RANDOM)
    if (madvise(pl->map, pl->len, MADV_RANDOM) < 0) {
        throw std::runtime_error("Error while calling madvise().");
        return NULL;  // Will leak memory
    }
#endif
    close(fd);
    return pl;
}

int spk_free(spkInfo *pl) {
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
    memset(pl, 0, sizeof(spkInfo));
    free(pl);
    return 0;
}

int spk_calc(spkInfo *pl, double epoch, int spiceId, double *out_x,
             double *out_y, double *out_z, double *out_vx, double *out_vy,
             double *out_vz, double *out_ax, double *out_ay, double *out_az) {
    if (pl == NULL) {
        throw std::runtime_error("The JPL ephemeris file has not been found.");
    }
    if (spiceId == 199) spiceId = 1;
    if (spiceId == 299) spiceId = 2;
    int m;
    m = spiceId;
    for (m = 0; m < pl->num; m++) {
        if (pl->targets[m].code == spiceId) {
            break;
        }
        if (m == pl->num - 1) {
            throw std::invalid_argument(
                "ERROR: Requested SPICE ID not found in SPK file");
        }
    }
    if (m < 0 || m >= pl->num) {
        throw std::runtime_error("The requested spice ID has not been found.");
    }
    spkTarget *target = &(pl->targets[m]);
    if (epoch < target->beg || epoch > target->end) {
        throw std::runtime_error(
            "The requested time is outside the coverage "
            "provided by the ephemeris file.");
    }
    double T[32];
    double S[32];
    double U[32];
    double *val, z;
    *out_x = 0.0;
    *out_y = 0.0;
    *out_z = 0.0;
    *out_vx = 0.0;
    *out_vy = 0.0;
    *out_vz = 0.0;
    *out_ax = 0.0;
    *out_ay = 0.0;
    *out_az = 0.0;
    if (target->cen == 3) {
        double xc, yc, zc, vxc, vyc, vzc, axc, ayc, azc;
        spk_calc(pl, epoch, target->cen, &xc, &yc, &zc, &vxc, &vyc, &vzc, &axc,
                 &ayc, &azc);
        *out_x = xc;
        *out_y = yc;
        *out_z = zc;
        *out_vx = vxc;
        *out_vy = vyc;
        *out_vz = vzc;
        *out_ax = axc;
        *out_ay = ayc;
        *out_az = azc;
    }
    int n, b, p, P, R;
    // find location of 'directory' describing the data records
    n = (int)((epoch - target->beg) / target->res);
    val = (double *)pl->map + target->two[n] - 1;
    // record size and number of coefficients per coordinate
    R = (int)val[-1];
    P = (R - 2) / 3;  // must be < 32 !!
    // pick out the precise record
    b = (int)((epoch - _mjd(val[-3])) / (val[-2] / 86400.0));
    val = (double *)pl->map + (target->one[n] - 1) + b * R;
    // scale to interpolation units
    z = (epoch - _mjd(val[0])) / (val[1] / 86400.0);
    // set up Chebyshev polynomials
    T[0] = 1.0;
    T[1] = z;
    S[0] = 0.0;
    S[1] = 1.0;
    U[0] = 0.0;
    U[1] = 0.0;
    U[2] = 4.0;
    for (p = 2; p < P; p++) {
        T[p] = 2.0 * z * T[p - 1] - T[p - 2];
        // Not used at the moment:
        S[p] = 2.0 * z * S[p - 1] + 2.0 * T[p - 1] - S[p - 2];
    }
    for (p = 3; p < P; p++) {
        U[p] = 2.0 * z * U[p - 1] + 4.0 * S[p - 1] - U[p - 2];
    }
    // double c = (0.125 / 2 / 86400.0);
    double t1 = 32.0;
    int niv;
    switch (spiceId) {
        case 1:
            niv = 4;
            break;
        case 2:
            niv = 2;
            break;
        case 3:
            niv = 2;
            break;
        case 4:
            niv = 1;
            break;
        case 5:
            niv = 1;
            break;
        case 6:
            niv = 1;
            break;
        case 7:
            niv = 1;
            break;
        case 8:
            niv = 1;
            break;
        case 9:
            niv = 1;
            break;
        case 10:
            niv = 2;
            break;
        case 301:
            niv = 8;
            break;
        case 399:
            niv = 8;
            break;
        default:
            niv = 1;
            break;
    }
    double c = (double)(niv * 2) / t1 / 86400.0;
    double pos[3] = {0.0, 0.0, 0.0};
    double vel[3] = {0.0, 0.0, 0.0};
    double acc[3] = {0.0, 0.0, 0.0};
    for (n = 0; n < 3; n++) {
        b = 2 + n * P;
        // sum interpolation stuff
        for (p = 0; p < P; p++) {
            pos[n] += val[b + p] * T[p] / 149597870.7;
            vel[n] += val[b + p] * S[p] * c / 149597870.7 * 86400.0;
            acc[n] += val[b + p] * U[p] * c * c / 149597870.7 * 86400.0 *
                      86400.0;
        }
    }
    *out_x += pos[0];
    *out_y += pos[1];
    *out_z += pos[2];
    *out_vx += vel[0];
    *out_vy += vel[1];
    *out_vz += vel[2];
    *out_ax += acc[0];
    *out_ay += acc[1];
    *out_az += acc[2];
    return 0;
}

void get_spk_state(const int &spiceID, const double &t0_mjd, Ephemeris &ephem,
                   double state[9]) {
    bool smallBody = spiceID > 1000000;
    spkInfo *infoToUse;
    if (smallBody) {
        infoToUse = ephem.sb;
    } else {
        infoToUse = ephem.mb;
    }
    // find what cache index corresponds to the requested SPICE ID
    int m;
    m = spiceID;
    for (m = 0; m < infoToUse->num; m++) {
        if (infoToUse->targets[m].code == spiceID) {
            break;
        }
        if (m == infoToUse->num - 1) {
            throw std::invalid_argument(
                "ERROR: Requested SPICE ID not found in SPK file");
        }
    }
    int cacheIdx = m;
    if (smallBody) {
        cacheIdx += ephem.mb->num;
    }
    // std::cout.precision(15);
    // std::cout << "cacheIdx = " << cacheIdx << ". ";
    // check if t0_mjd is in the ephem cache
    bool t0SomewhereInCache = false;
    for (size_t i = 0; i < ephem.cacheSize; i++) {
        if (ephem.cache[i].t == t0_mjd) {
            t0SomewhereInCache = true;
            if (ephem.cache[i].items[cacheIdx].t == t0_mjd &&
                ephem.cache[i].items[cacheIdx].spiceID == spiceID) {
                // std::cout << "Using cached state for " << spiceID << " at "
                // << t0_mjd << " from slot" << i << std::endl;
                state[0] = ephem.cache[i].items[cacheIdx].x;
                state[1] = ephem.cache[i].items[cacheIdx].y;
                state[2] = ephem.cache[i].items[cacheIdx].z;
                state[3] = ephem.cache[i].items[cacheIdx].vx;
                state[4] = ephem.cache[i].items[cacheIdx].vy;
                state[5] = ephem.cache[i].items[cacheIdx].vz;
                state[6] = ephem.cache[i].items[cacheIdx].ax;
                state[7] = ephem.cache[i].items[cacheIdx].ay;
                state[8] = ephem.cache[i].items[cacheIdx].az;
                return;
            }
        }
    }
    // if not, calculate it from the SPK memory map,
    double x, y, z, vx, vy, vz, ax, ay, az;
    spk_calc(infoToUse, t0_mjd, spiceID, &x, &y, &z, &vx, &vy, &vz, &ax, &ay,
             &az);
    state[0] = x;
    state[1] = y;
    state[2] = z;
    state[3] = vx;
    state[4] = vy;
    state[5] = vz;
    state[6] = ax;
    state[7] = ay;
    state[8] = az;
    if (smallBody) {
        double xSun, ySun, zSun, vxSun, vySun, vzSun, axSun, aySun, azSun;
        spk_calc(ephem.mb, t0_mjd, 10, &xSun, &ySun, &zSun, &vxSun, &vySun,
                 &vzSun, &axSun, &aySun, &azSun);
        state[0] += xSun;
        state[1] += ySun;
        state[2] += zSun;
        state[3] += vxSun;
        state[4] += vySun;
        state[5] += vzSun;
        state[6] += axSun;
        state[7] += aySun;
        state[8] += azSun;
    }
    // and add it to the cache
    if (!t0SomewhereInCache) {
        ephem.nextIdxToWrite++;
        if (ephem.nextIdxToWrite == ephem.cacheSize) {
            ephem.nextIdxToWrite = 0;
        }
    }
    // std::cout << "Adding state for " << spiceID << " at " << t0_mjd << " to
    // cache at slot" << ephem.nextIdxToWrite << std::endl;
    ephem.cache[ephem.nextIdxToWrite].t = t0_mjd;
    ephem.cache[ephem.nextIdxToWrite].items[cacheIdx].t = t0_mjd;
    ephem.cache[ephem.nextIdxToWrite].items[cacheIdx].spiceID = spiceID;
    ephem.cache[ephem.nextIdxToWrite].items[cacheIdx].x = state[0];
    ephem.cache[ephem.nextIdxToWrite].items[cacheIdx].y = state[1];
    ephem.cache[ephem.nextIdxToWrite].items[cacheIdx].z = state[2];
    ephem.cache[ephem.nextIdxToWrite].items[cacheIdx].vx = state[3];
    ephem.cache[ephem.nextIdxToWrite].items[cacheIdx].vy = state[4];
    ephem.cache[ephem.nextIdxToWrite].items[cacheIdx].vz = state[5];
    ephem.cache[ephem.nextIdxToWrite].items[cacheIdx].ax = state[6];
    ephem.cache[ephem.nextIdxToWrite].items[cacheIdx].ay = state[7];
    ephem.cache[ephem.nextIdxToWrite].items[cacheIdx].az = state[8];
}
