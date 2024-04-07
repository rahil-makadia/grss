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
    std::vector<real> raNutPrec, decNutPrec, wNutPrec, nutPrecIntercept, nutPrecSlope;
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
        // BODY199_NUT_PREC_RA  = ( 0. 0. 0. 0. 0. )
        // BODY199_NUT_PREC_DEC = ( 0. 0. 0. 0. 0. )
        // BODY199_NUT_PREC_PM  = (    0.01067257
        //                            -0.00112309
        //                            -0.00011040
        //                            -0.00002539
        //                            -0.00000571  )
        // BODY1_NUT_PREC_ANGLES  = ( 174.7910857  0.14947253587500003E+06
        //                            349.5821714  0.29894507175000006E+06
        //                            164.3732571  0.44841760762500006E+06
        //                            339.1643429  0.59789014350000012E+06
        //                            153.9554286  0.74736267937499995E+06 )
        cra0 = 281.0103;
        cra1 = -0.0328;
        cdec0 = 61.4155;
        cdec1 = -0.0049;
        cw0 = 329.5988;
        cw1 = 6.1385108;
        wNutPrec = {0.01067257, -0.00112309, -0.00011040, -0.00002539,
                    -0.00000571};
        nutPrecIntercept = {174.7910857, 349.5821714, 164.3732571, 339.1643429,
                            153.9554286};
        nutPrecSlope = {0.14947253587500003E+06, 0.29894507175000006E+06,
                        0.44841760762500006E+06, 0.59789014350000012E+06,
                        0.74736267937499995E+06};
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
        // BODY301_NUT_PREC_RA  = (   -3.8787   -0.1204   0.0700   -0.0172
        //                             0.0       0.0072   0.0       0.0
        //                             0.0      -0.0052   0.0       0.0
        //                             0.0043                              )

        // BODY301_NUT_PREC_DEC = (   1.5419     0.0239  -0.0278    0.0068
        //                            0.0       -0.0029   0.0009    0.0
        //                            0.0        0.0008   0.0       0.0
        //                           -0.0009                               )

        // BODY301_NUT_PREC_PM  = (   3.5610     0.1208  -0.0642    0.0158
        //                            0.0252    -0.0066  -0.0047   -0.0046
        //                            0.0028     0.0052   0.0040    0.0019
        //                           -0.0044                               )
        // BODY3_NUT_PREC_ANGLES  = (  125.045         -1935.5364525000
        //                             250.089         -3871.0729050000
        //                             260.008        475263.3328725000
        //                             176.625        487269.6299850000
        //                             357.529         35999.0509575000
        //                             311.589        964468.4993100000
        //                             134.963        477198.8693250000
        //                             276.617         12006.3007650000
        //                              34.226         63863.5132425000
        //                              15.134         -5806.6093575000
        //                             119.743           131.8406400000
        //                             239.961          6003.1503825000
        //                              25.053        473327.7964200000 )
        cra0 = 269.9949;
        cra1 = 0.0031;
        cdec0 = 66.5392;
        cdec1 = 0.0130;
        cw0 = 38.3213;
        cw1 = 13.17635815;
        cw2 = -1.4e-12;
        raNutPrec = {-3.8787, -0.1204, 0.0700,  -0.0172, 0.0, 0.0072, 0.0,
                     0.0,     0.0,     -0.0052, 0.0,     0.0, 0.0043};
        decNutPrec = {1.5419, 0.0239, -0.0278, 0.0068, 0.0, -0.0029, 0.0009,
                      0.0,    0.0,    0.0008,  0.0,    0.0, -0.0009};
        wNutPrec = {3.5610,  0.1208, -0.0642, 0.0158, 0.0252, -0.0066, -0.0047,
                    -0.0046, 0.0028, 0.0052,  0.0040, 0.0019, -0.0044};
        nutPrecIntercept = {125.045, 250.089, 260.008, 176.625, 357.529,
                            311.589, 134.963, 276.617, 34.226,  15.134,
                            119.743, 239.961, 25.053};
        nutPrecSlope = {-1935.5364525000,  -3871.0729050000, 475263.3328725000,
                        487269.6299850000, 35999.0509575000, 964468.4993100000,
                        477198.8693250000, 12006.3007650000, 63863.5132425000,
                        -5806.6093575000,  131.8406400000,   6003.1503825000,
                        473327.7964200000};
    } else if (iauFrame == "IAU_MARS"){
        // BODY499_POLE_RA          = (  317.269202  -0.10927547        0.  )
        // BODY499_POLE_DEC         = (   54.432516  -0.05827105        0.  )
        // BODY499_PM               = (  176.049863  +350.891982443297  0.  )
        // BODY499_NUT_PREC_RA      = (  0     0     0     0     0
        //                               0     0     0     0     0
        //                               0.000068
        //                               0.000238
        //                               0.000052
        //                               0.000009
        //                               0.419057                  )
        // BODY499_NUT_PREC_DEC     = (  0     0     0     0     0
        //                               0     0     0     0     0
        //                               0     0     0     0     0
        //                               0.000051
        //                               0.000141
        //                               0.000031
        //                               0.000005
        //                               1.591274                  )
        // BODY499_NUT_PREC_PM      = (  0     0     0     0     0
        //                               0     0     0     0     0
        //                               0     0     0     0     0
        //                               0     0     0     0     0
        //                               0.000145
        //                               0.000157
        //                               0.000040
        //                               0.000001
        //                               0.000001
        //                               0.584542                  )
        // BODY4_NUT_PREC_ANGLES  = (

        //     190.72646643      15917.10818695   0
        //     21.46892470       31834.27934054   0
        //     332.86082793      19139.89694742   0
        //     394.93256437      38280.79631835   0
        //     189.63271560   41215158.18420050   12.711923222

        //     121.46893664        660.22803474   0
        //     231.05028581        660.99123540   0
        //     251.37314025       1320.50145245   0
        //     217.98635955      38279.96125550   0
        //     196.19729402      19139.83628608   0

        //     198.991226        19139.4819985    0
        //     226.292679        38280.8511281    0
        //     249.663391        57420.7251593    0
        //     266.183510        76560.6367950    0
        //     79.398797             0.5042615    0

        //     122.433576        19139.9407476    0
        //     43.058401         38280.8753272    0
        //     57.663379         57420.7517205    0
        //     79.476401         76560.6495004    0
        //     166.325722            0.5042615    0

        //     129.071773        19140.0328244    0
        //     36.352167         38281.0473591    0
        //     56.668646         57420.9295360    0
        //     67.364003         76560.2552215    0
        //     104.792680        95700.4387578    0
        //     95.391654             0.5042615    0 )
        cra0 = 317.269202;
        cra1 = -0.10927547;
        cdec0 = 54.432516;
        cdec1 = -0.05827105;
        cw0 = 176.049863;
        cw1 = 350.891982443297;
        raNutPrec = {0.0, 0.0, 0.0, 0.0, 0.0, 
                     0.0, 0.0, 0.0, 0.0, 0.0,
                        0.000068, 0.000238, 0.000052, 0.000009, 0.419057};
        decNutPrec = {0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0, 0.0,
                        0.000051, 0.000141, 0.000031, 0.000005, 1.591274};
        wNutPrec = {0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0,
                    0.000145, 0.000157, 0.000040, 0.000001, 0.000001, 0.584542};
        nutPrecIntercept = {190.72646643, 21.46892470, 332.86082793, 394.93256437, 189.63271560,
                            121.46893664, 231.05028581, 251.37314025, 217.98635955, 196.19729402,
                            198.991226, 226.292679, 249.663391, 266.183510, 79.398797,
                            122.433576, 43.058401, 57.663379, 79.476401, 166.325722,
                            129.071773, 36.352167, 56.668646, 67.364003, 104.792680, 95.391654};
        nutPrecSlope = {15917.10818695, 31834.27934054, 19139.89694742, 38280.79631835, 41215158.18420050,
                        660.22803474, 660.99123540, 1320.50145245, 38279.96125550, 19139.83628608,
                        19139.4819985, 38280.8511281, 57420.7251593, 76560.6367950, 0.5042615,
                        19139.9407476, 38280.8753272, 57420.7517205, 76560.6495004, 0.5042615,
                        19140.0328244, 38281.0473591, 57420.9295360, 76560.2552215, 95700.4387578, 0.5042615};
    } else if (iauFrame == "IAU_JUPITER"){
        // BODY599_POLE_RA        = (   268.056595     -0.006499       0. )
        // BODY599_POLE_DEC       = (    64.495303      0.002413       0. )
        // BODY599_PM             = (   284.95        870.5360000      0. )
        // BODY599_NUT_PREC_RA  = ( 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.  0.000117
        //                                                         0.000938
        //                                                         0.001432
        //                                                         0.000030
        //                                                         0.002150 )
        // BODY599_NUT_PREC_DEC = ( 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.  0.000050
        //                                                         0.000404
        //                                                         0.000617
        //                                                        -0.000013
        //                                                         0.000926 )
        // BODY599_NUT_PREC_PM  = ( 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.  0.0
        //                                                         0.0
        //                                                         0.0
        //                                                         0.0
        //                                                         0.0  )
        // BODY5_NUT_PREC_ANGLES  = (    73.32      91472.9
        //                               24.62      45137.2
        //                              283.90       4850.7
        //                              355.80       1191.3
        //                              119.90        262.1
        //                              229.80         64.3
        //                              352.25       2382.6
        //                              113.35       6070.0
        //                              146.64     182945.8
        //                               49.24      90274.4
        //                               99.360714   4850.4046
        //                              175.895369   1191.9605
        //                              300.323162    262.5475
        //                              114.012305   6070.2476
        //                               49.511251     64.3000  )
        cra0 = 268.056595;
        cra1 = -0.006499;
        cdec0 = 64.495303;
        cdec1 = 0.002413;
        cw0 = 284.95;
        cw1 = 870.5360000;
        raNutPrec = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.000117, 0.000938, 0.001432, 0.000030, 0.002150};
        decNutPrec = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.000050, 0.000404, 0.000617, -0.000013, 0.000926};
        nutPrecIntercept = {73.32, 24.62, 283.90, 355.80, 119.90, 229.80, 352.25, 113.35,
                            146.64, 49.24, 99.360714, 175.895369, 300.323162, 114.012305,
                            49.511251};
        nutPrecSlope = {91472.9, 45137.2, 4850.7, 1191.3, 262.1, 64.3, 2382.6, 6070.0,
                        182945.8, 90274.4, 4850.4046, 1191.9605, 262.5475, 6070.2476,
                        64.3000};
    } else if (iauFrame == "IAU_SATURN"){
        // BODY699_POLE_RA        = (   40.589       -0.036       0. )
        // BODY699_POLE_DEC       = (   83.537       -0.004       0. )
        // BODY699_PM             = (   38.90       810.7939024   0. )
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
        cra0 = 257.311;
        cdec0 = -15.175;
        cw0 = 203.81;
        cw1 = -501.1600928;
    } else if (iauFrame == "IAU_NEPTUNE"){
        // BODY899_POLE_RA        = (  299.36     0.         0. )
        // BODY899_POLE_DEC       = (   43.46     0.         0. )
        // BODY899_PM             = (  249.978  541.1397757  0. )
        // BODY899_NUT_PREC_RA    = (  0.70 0. 0. 0. 0. 0. 0. 0. )
        // BODY899_NUT_PREC_DEC   = ( -0.51 0. 0. 0. 0. 0. 0. 0. )
        // BODY899_NUT_PREC_PM    = ( -0.48 0. 0. 0. 0. 0. 0. 0. )
        // BODY8_NUT_PREC_ANGLES = (   357.85         52.316
        //                             323.92      62606.6
        //                             220.51      55064.2
        //                             354.27      46564.5
        //                             75.31      26109.4
        //                             35.36      14325.4
        //                             142.61       2824.6
        //                             177.85         52.316
        //                             647.840    125213.200
        //                             355.700       104.632
        //                             533.550       156.948
        //                             711.400       209.264
        //                             889.250       261.580
        //                             1067.100       313.896
        //                             1244.950       366.212
        //                             1422.800       418.528
        //                             1600.650       470.844   )
        cra0 = 299.36;
        cdec0 = 43.46;
        cw0 = 249.978;
        cw1 = 541.1397757;
        raNutPrec = {0.70};
        decNutPrec = {-0.51};
        wNutPrec = {-0.48};
        nutPrecIntercept = {357.85, 323.92, 220.51, 354.27, 75.31, 35.36, 142.61, 177.85,
                            647.840, 355.700, 533.550, 711.400, 889.250, 1067.100, 1244.950,
                            1422.800, 1600.650};
        nutPrecSlope = {52.316, 62606.6, 55064.2, 46564.5, 26109.4, 14325.4, 2824.6, 52.316,
                        125213.200, 104.632, 156.948, 209.264, 261.580, 313.896, 366.212,
                        418.528, 470.844};
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
    // maxIter is the maximum of the lengths of raNutPrec, decNutPrec, wNutPrec
    size_t maxIter = std::max({raNutPrec.size(), decNutPrec.size(), wNutPrec.size()});
    for (size_t i = 0; i < maxIter; i++){
        const real theta = (nutPrecIntercept[i] + nutPrecSlope[i] * T) * DEG2RAD;
        const real cTheta = cos(theta);
        const real sTheta = sin(theta);
        const real thetaDot = nutPrecSlope[i] * DEG2RAD / secInCentury;
        if (i < raNutPrec.size()){
            ra += raNutPrec[i] * sTheta;
            raDot += raNutPrec[i] * cTheta * thetaDot;
        }
        if (i < decNutPrec.size()){
            dec += decNutPrec[i] * cTheta;
            decDot += decNutPrec[i] * -sTheta * thetaDot;
        }
        if (i < wNutPrec.size()){
            w += wNutPrec[i] * sTheta;
            wDot += wNutPrec[i] * cTheta * thetaDot;
        }
    }
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
    real angZ1 = euler[2];
    real cAngZ1 = cos(angZ1);
    real sAngZ1 = sin(angZ1);
    real angX = euler[1];
    real cAngX = cos(angX);
    real sAngX = sin(angX);
    real angZ2 = euler[0];
    real cAngZ2 = cos(angZ2);
    real sAngZ2 = sin(angZ2);
    rotMat[0] = cAngZ1*cAngZ2 - cAngX*sAngZ1*sAngZ2;
    rotMat[1] = cAngZ2*sAngZ1 + cAngX*cAngZ1*sAngZ2;
    rotMat[2] = sAngX*sAngZ2;
    rotMat[3] = -cAngZ1*sAngZ2 - cAngX*cAngZ2*sAngZ1;
    rotMat[4] = cAngX*cAngZ1*cAngZ2 - sAngZ1*sAngZ2;
    rotMat[5] = cAngZ2*sAngX;
    rotMat[6] = sAngX*sAngZ1;
    rotMat[7] = -cAngZ1*sAngX;
    rotMat[8] = cAngX;

    // from SPICE subroutine XF2EUL entry point EUL2XF for Rdot
    real ca, sa, u, v;
    ca = cos(euler[0]);
    sa = sin(euler[0]);
    u = cos(euler[1]);
    v = sin(euler[1]);
    real solutn[3][3] = {
        {-1.0, 0.0, -u},
        {0.0, -ca, -sa*v},
        {0.0, sa, -ca*v}
    };
    real domega[3];
    for (int i = 0; i < 3; i++){
        domega[i] = 0.0;
        for (int j = 0; j < 3; j++){
            domega[i] += solutn[i][j] * euler[j+3];
        }
    }
    real rotDotTimesRotTranspose[3][3] = {
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
