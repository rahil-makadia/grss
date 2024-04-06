#include "pck.h"

/**
 * @brief Calculate the Modified Julian Date from the Ephemeris Time.
 * 
 * @param[in] et Ephemeris Time (TDB seconds since J2000).
 * @return double Modified Julian Date (TDB).
 */
static double inline _mjd(double et) { return 51544.5 + et / 86400.0; }

/**
 * @param[in] bpc SpkInfo structure.
 */
void pck_free(PckInfo* bpc) {
    if (bpc == nullptr)
        return;

    if (bpc->targets){
        for (int m = 0; m < bpc->num; m++) {
            free(bpc->targets[m].one);
            free(bpc->targets[m].two);
        }
        free(bpc->targets);
    }
    munmap(bpc->map, bpc->len);
    memset(bpc, 0, sizeof(PckInfo));
    free(bpc);
}

/**
 * @param[in] path Path to the PCK file.
 * @return PckInfo* Pointer to the PckInfo structure for the PCK file.
 */
PckInfo* pck_init(const std::string &path) {
    // For file format information, see
    // https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/daf.html#The%20File%20Record
    // https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/pck.html#Binary%20PCK%20Kernel%20Format

    // Format for one summary
    struct summary {
        double beg;  // begin epoch, seconds since J2000.0
        double end;  // ending epoch
        int code;     // target frame code (3000 = ITRF93)
        int ref;     // reference frame (17 = ECLIPJ2000)
        int ver;     // type of PCK file
        int one;     // initial array address
        int two;     // final array address
    };

    // File is split into records. We read one record at a time.
    union {
        char buf[RECORD_LEN];
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
        return NULL;
    }

    // Read the file record.
    read(fd, &record, RECORD_LEN);
    // Check if the file is a valid double Precision Array File
    std::string full_file_type = "DAF/PCK";
    if (strncmp(record.file.locidw, full_file_type.c_str(), 7) != 0) {
        throw std::runtime_error(
            "Error parsing "+full_file_type+". Incorrect "
            "header.");
        close(fd);
        return NULL;
    }

    // Check that the size of a summary record is equal to the size of our
    // struct.
    int nc = 8 * (record.file.nd + (record.file.ni + 1) / 2);
    if (nc != sizeof(summary)) {
        throw std::runtime_error(
            "Error parsing "+full_file_type+". Wrong size of "
            "summary record.");
        close(fd);
        return NULL;
    }

    // Seek until the first summary record using the file record's fward pointer.
    // Record numbers start from 1 not 0 so we subtract 1 to get to the correct record.
    lseek(fd, (record.file.fward - 1) * RECORD_LEN, SEEK_SET);
    read(fd, record.buf, RECORD_LEN);

    // We are at the first summary block, validate
    if ((int64_t)record.buf[8] != 0) {
        throw std::runtime_error(
            "Error parsing "+full_file_type+". Cannot find "
            "summary block.");
        close(fd);
        return NULL;
    }
    // std::cout << "record.summaries.nsum: " << record.summaries.nsum << std::endl;
    // std::cout << "record.file.nd: " << record.file.nd << std::endl;
    // std::cout << "record.file.ni: " << record.file.ni << std::endl;
    // std::cout << "nc: " << nc << std::endl;
    // okay, here we go
    PckInfo *bpc = (PckInfo *)calloc(1, sizeof(PckInfo));
    // Loop over records
    while (true) {
        // Loop over summaries
        for (int b = 0; b < (int)record.summaries.nsum; b++) {
            // get current summary
            summary *sum = &record.summaries.s[b];
            // Index in our arrays for current target
            int m = bpc->num - 1;
            // New target?
            if (bpc->num == 0 || sum->code != bpc->targets[m].code) {
                if (bpc->num <= bpc->allocatedNum) {
                    bpc->allocatedNum += PCK_CACHE_ITEM_SIZE;  // increase space in batches of PCK_CACHE_ITEM_SIZE
                    bpc->targets = (PckTarget *)realloc(
                        bpc->targets, bpc->allocatedNum * sizeof(PckTarget));
                }
                m++;
                bpc->targets[m].code = sum->code;
                bpc->targets[m].ref = sum->ref;
                bpc->targets[m].beg = _mjd(sum->beg);
                bpc->targets[m].res = _mjd(sum->end) - bpc->targets[m].beg;
                bpc->targets[m].one = (int *)calloc(32768, sizeof(int));
                bpc->targets[m].two = (int *)calloc(32768, sizeof(int));
                bpc->targets[m].ind = 0;
                bpc->num++;
            }
            // add index for target
            bpc->targets[m].one[bpc->targets[m].ind] = sum->one;
            bpc->targets[m].two[bpc->targets[m].ind] = sum->two;
            bpc->targets[m].end = _mjd(sum->end);
            bpc->targets[m].ind++;
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
        throw std::runtime_error("Error calculating size for "+full_file_type+".");
        return NULL;
    }
    bpc->len = sb.st_size;

    // Memory map
    bpc->map = mmap(NULL, bpc->len, PROT_READ, MAP_SHARED, fd, 0);
    if (bpc->map == NULL) {
        throw std::runtime_error("Error creating memory map.");
        return NULL;  // Will leak memory
    }
    #if defined(MADV_RANDOM)
        if (madvise(bpc->map, bpc->len, MADV_RANDOM) < 0) {
            throw std::runtime_error("Error while calling madvise().");
            return NULL;  // Will leak memory
        }
    #endif
    close(fd);
    return bpc;
}

/**
 * @param[in] bpc SpkInfo structure.
 * @param[in] epoch Epoch to compute the state at (MJD ET).
 * @param[in] spiceId SPICE ID of the frame.
 * @param[in] inertialJ2000 True if the frame is inertial J2000.
 * @param[out] rotMat Rotation matrix
 * @param[out] rotMatDot Derivative of the rotation matrix
 */
void pck_calc(PckInfo *bpc, real epoch, int spiceId, bool inertialJ2000,
              real *rotMat, real *rotMatDot) {
    if (bpc == NULL) {
        throw std::runtime_error("The PCK file has not been found.");
    }
    if (spiceId != 3000) {
        throw std::runtime_error("The requested SPICE frame ID is currently not supported.");
    }
    int m;
    for (m = 0; m < bpc->num; m++) {
        if (bpc->targets[m].code == spiceId) {
            break;
        }
        if (m == bpc->num - 1) {
            throw std::invalid_argument("ERROR: Requested SPICE frame ID not found in PCK file");
        }
    }
    if (m < 0 || m >= bpc->num) {
        throw std::runtime_error("The requested SPICE frame ID has not been found.");
    }
    PckTarget *target = &(bpc->targets[m]);
    if (epoch < target->beg || epoch > target->end) {
        std::cout << "epoch: " << epoch << std::endl;
        std::cout << "target->beg: " << target->beg << std::endl;
        std::cout << "target->end: " << target->end << std::endl;
        throw std::runtime_error(
            "The requested time is outside the coverage "
            "provided by the PCK file.");
    }
    int n, b, p, P, R;
    // find location of 'directory' describing the data records
    n = (int)((epoch - target->beg) / target->res);
    double *val;
    val = (double *)bpc->map + target->two[n] - 1;
    // record size and number of coefficients per coordinate
    R = (int)val[-1];
    P = (R - 2) / 3;  // must be < 32 !!
    // pick out the precise record
    b = (int)((epoch - _mjd(val[-3])) / (val[-2] / 86400.0));
    val = (double *)bpc->map + (target->one[n] - 1) + b * R;
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
    for (p = 2; p < P; p++) {
        T[p] = 2.0 * z * T[p - 1] - T[p - 2];
        S[p] = 2.0 * z * S[p - 1] + 2.0 * T[p - 1] - S[p - 2];
    }
    double c = 1.0 / val[1]; // derivative scaling factor from SPICE/spke02.f and chbint.f
    double angle[3] = {0.0, 0.0, 0.0};
    double angleDot[3] = {0.0, 0.0, 0.0};
    for (n = 0; n < 3; n++) {
        b = 2 + n * P;
        // sum interpolation stuff
        for (p = 0; p < P; p++) {
            angle[n] += val[b + p] * T[p];
            angleDot[n] += val[b + p] * S[p] * c;
        }
    }
    // phi, delta, w 3-1-3 euler angles are in a specific order
    real euler[6];
    euler[0] = (real)angle[2];
    euler[1] = (real)angle[1];
    euler[2] = (real)angle[0];
    euler[3] = (real)angleDot[2];
    euler[4] = (real)angleDot[1];
    euler[5] = (real)angleDot[0];

    real rotMatEclip[3][3], rotMatEclipDot[3][3];
    euler313_to_rotmat(euler, *rotMatEclip, *rotMatEclipDot);
    
    if (inertialJ2000) {
        // We have ecliptic to body rotation matrix and its derivative
        // We need equatorial to body rotation matrix and its derivative
        real rotEquatEclipTranspose[3][3] = {
            {1.0, 0.0, 0.0},
            {0.0, cos(EARTH_OBLIQUITY), sin(EARTH_OBLIQUITY)},
            {0.0, -sin(EARTH_OBLIQUITY), cos(EARTH_OBLIQUITY)}
        };
        mat3_mat3_mul(*rotMatEclip, *rotEquatEclipTranspose, rotMat);
        mat3_mat3_mul(*rotMatEclipDot, *rotEquatEclipTranspose, rotMatDot);
    } else {
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                rotMat[3*i+j] = rotMatEclip[i][j];
                rotMatDot[3*i+j] = rotMatEclipDot[i][j];
            }
        }
    }
}

void euler313_to_rotmat(const real euler[6], real *rotMat, real *rotMatDot){
    // convert euler angles to rotation matrix
    // (we need w, delta, phi 3-1-3 euler angles here)
    double angZ1 = euler[2];
    double angX = euler[1];
    double angZ2 = euler[0];
    double rz2[3][3] = {
        {cos(angZ2), sin(angZ2), 0.0},
        {-sin(angZ2), cos(angZ2), 0.0},
        {0.0, 0.0, 1.0}
    };
    double rx[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, cos(angX), sin(angX)},
        {0.0, -sin(angX), cos(angX)}
    };
    double rz1[3][3] = {
        {cos(angZ1), sin(angZ1), 0.0},
        {-sin(angZ1), cos(angZ1), 0.0},
        {0.0, 0.0, 1.0}
    };
    double rotMatTemp[3][3];
    mat3_mat3_mul(*rx, *rz1, *rotMatTemp);
    mat3_mat3_mul(*rz2, *rotMatTemp, rotMat);

    // from SPICE subroutine XF2EUL entry point EUL2XF for Rdot
    double ca, sa, u, v;
    ca = cos(euler[0]);
    sa = sin(euler[0]);
    u = cos(euler[1]);
    v = sin(euler[1]);
    double solutn[3][3] = {
        {-1.0, 0.0, -u},
        {0.0, -ca, -sa*v},
        {0.0, sa, -ca*v}
    };
    double domega[3];
    for (int i = 0; i < 3; i++){
        domega[i] = 0.0;
        for (int j = 0; j < 3; j++){
            domega[i] += solutn[i][j] * euler[j+3];
        }
    }
    double rotDotTimesRotTranspose[3][3] = {
        {0.0, -domega[0], domega[2]},
        {domega[0], 0.0, -domega[1]},
        {-domega[2], domega[1], 0.0}
    };
    mat3_mat3_mul(*rotDotTimesRotTranspose, rotMat, rotMatDot);
}

/**
 * @param from Frame to rotate from.
 * @param to Frame to rotate to.
 * @param et TDB ephemeris time in seconds past J2000 to get the rotation matrix at.
 * @param rotMat Rotation matrix from 'from' to 'to'.
 */
void get_pck_rotMat(const std::string &from, const std::string &to,
                    const double &et, std::vector<std::vector<double>> &rotMat) {
    SpiceDouble rotMatSpice[6][6];
    sxform_c(from.c_str(), to.c_str(), et, rotMatSpice);
    for (size_t i = 0; i < 6; i++) {
        for (size_t j = 0; j < 6; j++) {
            rotMat[i][j] = (double)rotMatSpice[i][j];
        }
    }
}

/**
 * @param from Frame to rotate from.
 * @param to Frame to rotate to.
 * @param t0_mjd t0_mjd Epoch to compute the rotation matrix at (MJD ET).
 * @param bpc PckInfo structure.
 * @param xformMat Rotation matrix from 'from' to 'to'.
 */
void get_pck_rotMat2(const std::string &from, const std::string &to,
                    const real &t0_mjd, PckInfo* bpc,//PckEphemeris &ephem,
                    std::vector<std::vector<real>> &xformMat) {
    int spiceId = 3000;
    bool inertialJ2000 = true;
    bool bodyToInertial = true;
    if (from == "ITRF93") {
        if (to == "ECLIPJ2000") {
            inertialJ2000 = false;
        } else if (to != "J2000") {
            throw std::runtime_error("The requested frame transformation is currently not supported.");
        }
    } else if (to == "ITRF93") {
        bodyToInertial = false;
        if (from == "ECLIPJ2000") {
            inertialJ2000 = false;
        } else if (from != "J2000") {
            throw std::runtime_error("The requested frame transformation is currently not supported.");
        }
    } else {
        throw std::runtime_error("The requested frame transformation is currently not supported.");
    }
    // PckInfo *bpc = ephem.histPck;
    real rotMat[3][3], rotMatDot[3][3];
    pck_calc(bpc, t0_mjd, spiceId, inertialJ2000, *rotMat, *rotMatDot);
    // Okay, assemble the state rotation matrix
    if (bodyToInertial){
        // if bodyToInertial is true, then rotMat and rotMatDot need to be transposed
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                xformMat[i][j] = rotMat[j][i];
                xformMat[i+3][j] = rotMatDot[j][i];
                xformMat[i+3][j+3] = rotMat[j][i];
            }
        }
    } else {
        // we already built the inertial to body rotation matrices, copy them
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                xformMat[i][j] = rotMat[i][j];
                xformMat[i+3][j] = rotMatDot[i][j];
                xformMat[i+3][j+3] = rotMat[i][j];
            }
        }
    }
}
