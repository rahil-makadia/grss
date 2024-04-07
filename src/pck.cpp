#include "pck.h"

/**
 * @brief Calculate the Modified Julian Date from the Ephemeris Time.
 * 
 * @param[in] et Ephemeris Time (TDB seconds since J2000).
 * @return double Modified Julian Date (TDB).
 */
static double inline _mjd(double et) { return 51544.5 + et / 86400.0; }

/**
 * @param[in] bpc PckInfo structure.
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
        throw std::runtime_error("Error opening "+path+".");
    }

    // Read the file record.
    read(fd, &record, RECORD_LEN);
    // Check if the file is a valid double Precision Array File
    std::string full_file_type = "DAF/PCK";
    if (strncmp(record.file.locidw, full_file_type.c_str(), 7) != 0) {
        close(fd);
        throw std::runtime_error(
            "Error parsing "+full_file_type+". Incorrect "
            "header.");
    }

    // Check that the size of a summary record is equal to the size of our
    // struct.
    int nc = 8 * (record.file.nd + (record.file.ni + 1) / 2);
    if (nc != sizeof(summary)) {
        close(fd);
        throw std::runtime_error(
            "Error parsing "+full_file_type+". Wrong size of "
            "summary record.");
    }

    // Seek until the first summary record using the file record's fward pointer.
    // Record numbers start from 1 not 0 so we subtract 1 to get to the correct record.
    lseek(fd, (record.file.fward - 1) * RECORD_LEN, SEEK_SET);
    read(fd, record.buf, RECORD_LEN);

    // We are at the first summary block, validate
    if ((int64_t)record.buf[8] != 0) {
        close(fd);
        throw std::runtime_error(
            "Error parsing "+full_file_type+". Cannot find "
            "summary block.");
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
                    bpc->allocatedNum += 4;  // increase space in batches of PCK_CACHE_ITEM_SIZE(4)
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
    }
    bpc->len = sb.st_size;

    // Memory map
    bpc->map = mmap(NULL, bpc->len, PROT_READ, MAP_SHARED, fd, 0);
    if (bpc->map == NULL) {
        // this will leak memory
        throw std::runtime_error("Error creating memory map.");
    }
    #if defined(MADV_RANDOM)
        if (madvise(bpc->map, bpc->len, MADV_RANDOM) < 0) {
            // this will leak memory
            throw std::runtime_error("Error while calling madvise().");
        }
    #endif
    close(fd);
    return bpc;
}

/**
 * @param[in] bpc PckInfo structure.
 * @param[in] epoch Epoch to compute the state at (MJD ET).
 * @param[in] spiceId SPICE ID of the frame.
 * @param[out] rotMat Rotation matrix
 * @param[out] rotMatDot Derivative of the rotation matrix
 */
void pck_calc(PckInfo *bpc, real epoch, int spiceId, real *rotMat,
              real *rotMatDot) {
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
    // phi, delta, w are computed from chebyshev polynomials
    // but we need w, delta, phi for the euler313_to_rotmat function
    real euler[6];
    euler[0] = (real)angle[2];
    euler[1] = (real)angle[1];
    euler[2] = (real)angle[0];
    euler[3] = (real)angleDot[2];
    euler[4] = (real)angleDot[1];
    euler[5] = (real)angleDot[0];

    real rotMatEclip[3][3], rotMatEclipDot[3][3];
    euler313_to_rotMat(euler, *rotMatEclip, *rotMatEclipDot);
    // We have ecliptic to body rotation matrix and its derivative
    // We need equatorial to body rotation matrix and its derivative
    real rotEquatEclipTranspose[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, cos(EARTH_OBLIQUITY), sin(EARTH_OBLIQUITY)},
        {0.0, -sin(EARTH_OBLIQUITY), cos(EARTH_OBLIQUITY)}
    };
    mat3_mat3_mul(*rotMatEclip, *rotEquatEclipTranspose, rotMat);
    mat3_mat3_mul(*rotMatEclipDot, *rotEquatEclipTranspose, rotMatDot);
}

/**
 * @param[in] t0_mjd Epoch to compute the rotation matrix at (MJD ET).
 * @param[in] iauFrame IAU frame name.
 * @param[out] euler Real array of 6 elements containing the 313 Euler angles and their derivatives.
 */
void iau_to_euler(const real t0_mjd, std::string iauFrame, real *euler){
    real secInDay = 86400.0;
    real secInCentury = 36525.0 * secInDay;
    real D = t0_mjd - 51544.5;
    real T = D / 36525.0;
    real cra0, cra1, cra2;
    cra0 = cra1 = cra2 = 0.0;
    real cdec0, cdec1, cdec2;
    cdec0 = cdec1 = cdec2 = 0.0;
    real cw0, cw1, cw2;
    cw0 = cw1 = cw2 = 0.0;
    // numbers are from NAIF's pck00011.tpc file
    if (iauFrame == "IAU_SUN"){
        // BODY10_POLE_RA         = (  286.13       0.          0. )
        // BODY10_POLE_DEC        = (   63.87       0.          0. )
        // BODY10_PM              = (   84.176     14.18440     0. )
        cra0 = 286.13;
        cdec0 = 63.87;
        cw0 = 84.176;
        cw1 = 14.1844;
    } else if (iauFrame == "IAU_MERCURY"){
        // BODY199_POLE_RA          = (  281.0103   -0.0328     0. )
        // BODY199_POLE_DEC         = (   61.4155   -0.0049     0. )
        // BODY199_PM               = (  329.5988    6.1385108  0. )
        std::cout << "Note: the IAU Mercury model has nutation and precision "
                     "terms that are not included in GRSS."
                  << std::endl;
        cra0 = 281.0103;
        cra1 = -0.0328;
        cdec0 = 61.4155;
        cdec1 = -0.0049;
        cw0 = 329.5988;
        cw1 = 6.1385108;
    } else if (iauFrame == "IAU_VENUS"){
        // BODY299_POLE_RA          = (  272.76       0.          0. )
        // BODY299_POLE_DEC         = (   67.16       0.          0. )
        // BODY299_PM               = (  160.20      -1.4813688   0. )
        cra0 = 272.76;
        cdec0 = 67.16;
        cw0 = 160.20;
        cw1 = -1.4813688;
    } else if (iauFrame == "IAU_EARTH"){
        // BODY399_POLE_RA        = (    0.      -0.641         0. )
        // BODY399_POLE_DEC       = (   90.      -0.557         0. )
        // BODY399_PM             = (  190.147  360.9856235     0. )
        cra1 = -0.641;
        cdec0 = 90.0;
        cdec1 = -0.557;
        cw0 = 190.147;
        cw1 = 360.9856235;
    } else if (iauFrame == "IAU_MOON") {
        // BODY301_POLE_RA      = (  269.9949        0.0031        0.      )
        // BODY301_POLE_DEC     = (   66.5392        0.0130        0.      )
        // BODY301_PM           = (   38.3213       13.17635815   -1.4D-12 )
        std::cout << "Note: the IAU Moon model has nutation and precision "
                     "terms that are not included in GRSS."
                  << std::endl;
        cra0 = 269.9949;
        cra1 = 0.0031;
        cdec0 = 66.5392;
        cdec1 = 0.0130;
        cw0 = 38.3213;
        cw1 = 13.17635815;
        cw2 = -1.4e-12;
    } else if (iauFrame == "IAU_MARS"){
        // BODY499_POLE_RA          = (  317.269202  -0.10927547        0.  )
        // BODY499_POLE_DEC         = (   54.432516  -0.05827105        0.  )
        // BODY499_PM               = (  176.049863  +350.891982443297  0.  )
        std::cout << "Note: the IAU Mars model has nutation and precision "
                     "terms that are not included in GRSS."
                  << std::endl;
        cra0 = 317.269202;
        cra1 = -0.10927547;
        cdec0 = 54.432516;
        cdec1 = -0.05827105;
        cw0 = 176.049863;
        cw1 = 350.891982443297;
    } else if (iauFrame == "IAU_JUPITER"){
        // BODY599_POLE_RA        = (   268.056595     -0.006499       0. )
        // BODY599_POLE_DEC       = (    64.495303      0.002413       0. )
        // BODY599_PM             = (   284.95        870.5360000      0. )
        std::cout << "Note: the IAU Jupiter model has nutation and precision "
                     "terms that are not included in GRSS."
                  << std::endl;
        cra0 = 268.056595;
        cra1 = -0.006499;
        cdec0 = 64.495303;
        cdec1 = 0.002413;
        cw0 = 284.95;
        cw1 = 870.5360000;
    } else if (iauFrame == "IAU_SATURN"){
        // BODY699_POLE_RA        = (   40.589       -0.036       0. )
        // BODY699_POLE_DEC       = (   83.537       -0.004       0. )
        // BODY699_PM             = (   38.90       810.7939024   0. )
        std::cout << "Note: the IAU Saturn model has nutation and precision "
                     "terms that are not included in GRSS."
                  << std::endl;
        cra0 = 40.589;
        cra1 = -0.036;
        cdec0 = 83.537;
        cdec1 = -0.004;
        cw0 = 38.90;
        cw1 = 810.7939024;
    } else if (iauFrame == "IAU_URANUS"){
        // BODY799_POLE_RA        = (  257.311     0.         0.  )
        // BODY799_POLE_DEC       = (  -15.175     0.         0.  )
        // BODY799_PM             = (  203.81   -501.1600928  0.  )
        std::cout << "Note: the IAU Uranus model has nutation and precision "
                     "terms that are not included in GRSS."
                  << std::endl;
        cra0 = 257.311;
        cdec0 = -15.175;
        cw0 = 203.81;
        cw1 = -501.1600928;
    } else if (iauFrame == "IAU_NEPTUNE"){
        // BODY899_POLE_RA        = (  299.36     0.         0. )
        // BODY899_POLE_DEC       = (   43.46     0.         0. )
        // BODY899_PM             = (  249.978  541.1397757  0. )
        std::cout << "Note: the IAU Neptune model has nutation and precision "
                     "terms that are not included in GRSS."
                  << std::endl;
        cra0 = 299.36;
        cdec0 = 43.46;
        cw0 = 249.978;
        cw1 = 541.1397757;
    } else if (iauFrame == "IAU_PLUTO"){
        // BODY999_POLE_RA        = (  132.993   0.          0. )
        // BODY999_POLE_DEC       = (   -6.163   0.          0. )
        // BODY999_PM             = (  302.695   56.3625225  0. )
        cra0 = 132.993;
        cdec0 = -6.163;
        cw0 = 302.695;
        cw1 = 56.3625225;
    } else {
        throw std::runtime_error("iau_to_euler: The IAU frame is not supported.");
    }
    real ra, dec, w, raDot, decDot, wDot;
    ra = cra0 + T * (cra1 + T * cra2);
    dec = cdec0 + T * (cdec1 + T * cdec2);
    w = cw0 + D * (cw1 + D * cw2);
    raDot = (cra1 + 2.0 * T * cra2) / secInCentury;
    decDot = (cdec1 + 2.0 * T * cdec2) / secInCentury;
    wDot = (cw1 + 2.0 * D * cw2) / secInDay;
    w *= DEG2RAD;
    ra *= DEG2RAD;
    dec *= DEG2RAD;
    wDot *= DEG2RAD;
    raDot *= DEG2RAD;
    decDot *= DEG2RAD;
    euler[0] = fmod(w, 2*PI); // w
    euler[1] = PI/2 - dec; // delta
    euler[2] = ra + PI/2; // phi
    euler[3] = wDot; // wDot
    euler[4] = -decDot; // deltaDot
    euler[5] = raDot; // phiDot
}

/**
 * @param[in] euler Real array of 6 elements containing the 313 Euler angles and their derivatives.
 * @param[out] rotMat Real array of 9 elements to store the rotation matrix.
 * @param[out] rotMatDot Real array of 9 elements to store the rotation matrix derivative.
 */
void euler313_to_rotMat(const real euler[6], real *rotMat, real *rotMatDot){
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
 * @param[in] from Frame to rotate from.
 * @param[in] to Frame to rotate to.
 * @param[in] t0_mjd t0_mjd Epoch to compute the rotation matrix at (MJD ET).
 * @param[in] ephem PckEphemeris structure.
 * @param[out] xformMat Rotation matrix from 'from' to 'to'.
 */
void get_pck_rotMat(const std::string &from, const std::string &to,
                    const real &t0_mjd, PckEphemeris &ephem,
                    std::vector<std::vector<real>> &xformMat) {
    // trivial case
    if (from == to){
        for (int i = 0; i < 6; i++){
            for (int j = 0; j < 6; j++){
                xformMat[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
        return;
    }
    bool bodyToInertial = true;
    std::vector<std::string> fromTo = {from, to};
    std::vector<std::string> validBodyFrames = {"ITRF93", "IAU_SUN",
                                                "IAU_MERCURY", "IAU_VENUS",
                                                "IAU_EARTH", "IAU_MOON",
                                                "IAU_MARS", "IAU_JUPITER",
                                                "IAU_SATURN", "IAU_URANUS",
                                                "IAU_NEPTUNE", "IAU_PLUTO"};
    // make sure either from or to frame is J2000
    int bodyFrameIdx = -1;
    if (from == "J2000"){
        bodyToInertial = false;
        bodyFrameIdx = 1;
    } else if (to == "J2000"){
        bodyFrameIdx = 0;
    } else {
        throw std::runtime_error("get_pck_rotMat: The inertial frame is not J2000.");
    }
    bool nonBinaryBodyFrame = false;
    // and the other frame is from a list of valid frames
    if (std::find(validBodyFrames.begin(), validBodyFrames.end(), fromTo[bodyFrameIdx]) == validBodyFrames.end()){
        std::cout << "Body frame: " << fromTo[bodyFrameIdx] << std::endl;
        throw std::runtime_error("get_pck_rotMat: The body frame is not valid.");
    }

    if (fromTo[bodyFrameIdx].substr(0, 4) == "IAU_"){
        nonBinaryBodyFrame = true;
    }
    real rotMat[3][3], rotMatDot[3][3];
    if (nonBinaryBodyFrame){
        real euler[6];
        iau_to_euler(t0_mjd, fromTo[bodyFrameIdx], euler);
        euler313_to_rotMat(euler, *rotMat, *rotMatDot);
    } else {
        if (ephem.histPck == nullptr || ephem.latestPck == nullptr || ephem.predictPck == nullptr){
        throw std::invalid_argument(
            "get_pck_rotMat: PCK kernels are not loaded. Memory map "
            "the ephemeris using PropSimulation.map_ephemeris() method first.");
        }
        // pick the correct PCK file memory map based on the epoch
        PckInfo *bpc;
        if (t0_mjd < ephem.histPck->targets[0].beg || t0_mjd > ephem.predictPck->targets[0].end){
            throw std::runtime_error("get_pck_rotMat: The epoch is outside the range of the binary PCK files.");
        }
        if (t0_mjd <= ephem.histPck->targets[0].end){
            bpc = ephem.histPck;
            // std::cout << "Using histPck" << std::endl;
        } else if (t0_mjd <= ephem.latestPck->targets[0].end){
            bpc = ephem.latestPck;
            // std::cout << "Using latestPck" << std::endl;
        } else {
            bpc = ephem.predictPck;
            // std::cout << "Using predictPck" << std::endl;
        }
        int BinaryFrameSpiceId;
        if (fromTo[bodyFrameIdx] == "ITRF93"){
            BinaryFrameSpiceId = 3000;
        } else {
            throw std::runtime_error("get_pck_rotMat: The binary body frame is not recognized.");
        }
        pck_calc(bpc, t0_mjd, BinaryFrameSpiceId, *rotMat, *rotMatDot);
    }
    // Okay, assemble the state rotation matrix
    if (bodyToInertial){
        // if bodyToInertial is true, then rotMat and rotMatDot need to be transposed
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                xformMat[i][j] = xformMat[i+3][j+3] = rotMat[j][i];
                xformMat[i+3][j] = rotMatDot[j][i];
            }
        }
    } else {
        // we already built the inertial to body rotation matrices, copy them
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                xformMat[i][j] = xformMat[i+3][j+3] = rotMat[i][j];
                xformMat[i+3][j] = rotMatDot[i][j];
            }
        }
    }
}
