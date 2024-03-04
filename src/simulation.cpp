#include "simulation.h"

void get_observer_state(const real &tObsMjd,
                        const std::vector<real> &observerInfo,
                        propSimulation *propSim, const bool tObsInUTC,
                        std::vector<real> &observerState) {
    SpiceInt baseBody = observerInfo[0];
    if (observerInfo[0] == 500) baseBody = 399;
    if (baseBody == 0) {
        observerState[0] = 0.0L;
        observerState[1] = 0.0L;
        observerState[2] = 0.0L;
        observerState[3] = 0.0L;
        observerState[4] = 0.0L;
        observerState[5] = 0.0L;
        return;
    }
    real tObsET;
    real tObsMjdTDB;
    mjd_to_et(tObsMjd, tObsET);
    if (tObsInUTC) {
        // std::cout << tObsMjd << " MJD UTC" << tObsET << " s UTC -> ";
        SpiceDouble et_minus_utc;
        real sec_past_j2000_utc = tObsET;
        deltet_c(sec_past_j2000_utc, "UTC", &et_minus_utc);
        tObsET += et_minus_utc;
        et_to_mjd(tObsET, tObsMjdTDB);
        // std::cout << tObsET << " s ET " << tObsMjdTDB << " MJD TDB" <<
        // std::endl;
    } else {
        tObsMjdTDB = tObsMjd;
    }
    double baseBodyState[9];
    get_spk_state(baseBody, tObsMjdTDB, propSim->ephem, baseBodyState);
    ConstSpiceChar *baseBodyFrame;
    switch ((int)observerInfo[0]) {
        case 10:
            baseBodyFrame = "IAU_SUN";
            break;
        case 1:
        case 199:
            baseBodyFrame = "IAU_MERCURY";
            break;
        case 2:
        case 299:
            baseBodyFrame = "IAU_VENUS";
            break;
        case 399:
            baseBodyFrame = "ITRF93";
            // High precision frame is not defined before 1972 JAN 01 00:00:42.183 TDB
            if (tObsMjdTDB < 41317.0004882291666666666L) {
                baseBodyFrame = "IAU_EARTH";
            }
            break;
        case 500:
            observerState[0] = (real) baseBodyState[0] + observerInfo[1]/propSim->consts.du2m;
            observerState[1] = (real) baseBodyState[1] + observerInfo[2]/propSim->consts.du2m;
            observerState[2] = (real) baseBodyState[2] + observerInfo[3]/propSim->consts.du2m;
            observerState[3] = (real) baseBodyState[3] + observerInfo[4]/propSim->consts.duptu2mps;
            observerState[4] = (real) baseBodyState[4] + observerInfo[5]/propSim->consts.duptu2mps;
            observerState[5] = (real) baseBodyState[5] + observerInfo[6]/propSim->consts.duptu2mps;
            return;
            break;
        case 499:
            baseBodyFrame = "IAU_MARS";
            break;
        case 599:
            baseBodyFrame = "IAU_JUPITER";
            break;
        case 699:
            baseBodyFrame = "IAU_SATURN";
            break;
        case 799:
            baseBodyFrame = "IAU_URANUS";
            break;
        case 899:
            baseBodyFrame = "IAU_NEPTUNE";
            break;
        case 999:
            baseBodyFrame = "IAU_PLUTO";
            break;
        default:
            std::cout << "Given base body: " << baseBody << std::endl;
            throw std::invalid_argument("Given base body not supported");
            break;
    }
    real lon = observerInfo[1];
    real lat = observerInfo[2];
    real rho = observerInfo[3];
    ConstSpiceDouble bodyFixedX = rho * cos(lat) * cos(lon) / 1.0e3L;
    ConstSpiceDouble bodyFixedY = rho * cos(lat) * sin(lon) / 1.0e3L;
    ConstSpiceDouble bodyFixedZ = rho * sin(lat) / 1.0e3L;
    ConstSpiceDouble bodyFixedState[6] = {bodyFixedX, bodyFixedY, bodyFixedZ,
                                          0.0,        0.0,        0.0};
    SpiceDouble observerStateInertial[6];
    SpiceDouble rotMat[6][6];
    sxform_c(baseBodyFrame, "J2000", tObsET, rotMat);
    mxvg_c(rotMat, bodyFixedState, 6, 6, observerStateInertial);
    observerStateInertial[0] *= (real)1.0e3L / propSim->consts.du2m;
    observerStateInertial[1] *= (real)1.0e3L / propSim->consts.du2m;
    observerStateInertial[2] *= (real)1.0e3L / propSim->consts.du2m;
    observerStateInertial[3] *= (real)1.0e3L / propSim->consts.duptu2mps;
    observerStateInertial[4] *= (real)1.0e3L / propSim->consts.duptu2mps;
    observerStateInertial[5] *= (real)1.0e3L / propSim->consts.duptu2mps;
    observerState[0] = baseBodyState[0] + observerStateInertial[0];
    observerState[1] = baseBodyState[1] + observerStateInertial[1];
    observerState[2] = baseBodyState[2] + observerStateInertial[2];
    observerState[3] = baseBodyState[3] + observerStateInertial[3];
    observerState[4] = baseBodyState[4] + observerStateInertial[4];
    observerState[5] = baseBodyState[5] + observerStateInertial[5];
}

void Body::set_J2(real J2, real poleRA, real poleDec) {
    this->J2 = J2;
    if (this->J2 != 0.0L) {
        this->isJ2 = true;
    } else {
        this->isJ2 = false;
    }
    this->poleRA = poleRA * DEG2RAD;
    this->poleDec = poleDec * DEG2RAD;
}

SpiceBody::SpiceBody(std::string name, int spiceId, real t0, real mass,
                     real radius) {
    this->name = name;
    this->spiceId = spiceId;
    if (this->spiceId > 1000000) {
        this->caTol = 0.05;
    }
    this->t0 = t0;
    this->mass = mass;
    this->radius = radius;
    this->isNongrav = false;
    this->isPPN = false;
    this->isMajor = false;
    // this->pos = {0.0L, 0.0L, 0.0L};
    this->pos[0] = 0.0L;
    this->pos[1] = 0.0L;
    this->pos[2] = 0.0L;
    // this->vel = {0.0L, 0.0L, 0.0L};
    this->vel[0] = 0.0L;
    this->vel[1] = 0.0L;
    this->vel[2] = 0.0L;
    // this->acc = {0.0L, 0.0L, 0.0L};
    this->acc[0] = 0.0L;
    this->acc[1] = 0.0L;
    this->acc[2] = 0.0L;
}

IntegBody::IntegBody(std::string name, real t0, real mass, real radius,
                     std::vector<real> cometaryState,
                     NongravParamaters ngParams) {
    this->name = name;
    this->t0 = t0;
    this->mass = mass;
    this->radius = radius;
    this->caTol = 0.0;
    std::vector<real> cartesianStateEclip(6);
    std::vector<real> cartesianPos(3);
    std::vector<real> cartesianVel(3);
    this->isCometary = true;
    this->initState = cometaryState;

    cometary_to_cartesian(t0, cometaryState, cartesianStateEclip);
    // rotate to eme2000
    std::vector<std::vector<real>> eclipToEquatorial(3, std::vector<real>(3));
    rot_mat_x(EARTH_OBLIQUITY, eclipToEquatorial);
    mat_vec_mul(eclipToEquatorial,
                {cartesianStateEclip[0], cartesianStateEclip[1],
                 cartesianStateEclip[2]},
                cartesianPos);
    mat_vec_mul(eclipToEquatorial,
                {cartesianStateEclip[3], cartesianStateEclip[4],
                 cartesianStateEclip[5]},
                cartesianVel);
    this->pos[0] = cartesianPos[0];
    this->pos[1] = cartesianPos[1];
    this->pos[2] = cartesianPos[2];
    this->vel[0] = cartesianVel[0];
    this->vel[1] = cartesianVel[1];
    this->vel[2] = cartesianVel[2];
    this->acc[0] = 0.0L;
    this->acc[1] = 0.0L;
    this->acc[2] = 0.0L;
    this->isNongrav = false;
    if (ngParams.a1 != 0.0L || ngParams.a2 != 0.0L || ngParams.a3 != 0.0L) {
        this->ngParams.a1 = ngParams.a1;
        this->ngParams.a2 = ngParams.a2;
        this->ngParams.a3 = ngParams.a3;
        this->ngParams.alpha = ngParams.alpha;
        this->ngParams.k = ngParams.k;
        this->ngParams.m = ngParams.m;
        this->ngParams.n = ngParams.n;
        this->ngParams.r0_au = ngParams.r0_au;
        this->isNongrav = true;
    }
    this->isPPN = false;
    this->isMajor = false;
}

IntegBody::IntegBody(std::string name, real t0, real mass, real radius,
                     std::vector<real> pos, std::vector<real> vel,
                     NongravParamaters ngParams) {
    this->name = name;
    this->t0 = t0;
    this->mass = mass;
    this->radius = radius;
    this->caTol = 0.0;
    this->isCometary = false;
    this->initState = {pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]};
    this->pos[0] = pos[0];
    this->pos[1] = pos[1];
    this->pos[2] = pos[2];
    this->vel[0] = vel[0];
    this->vel[1] = vel[1];
    this->vel[2] = vel[2];
    this->acc[0] = 0.0L;
    this->acc[1] = 0.0L;
    this->acc[2] = 0.0L;
    this->isNongrav = false;
    if (ngParams.a1 != 0.0L || ngParams.a2 != 0.0L || ngParams.a3 != 0.0L) {
        this->ngParams.a1 = ngParams.a1;
        this->ngParams.a2 = ngParams.a2;
        this->ngParams.a3 = ngParams.a3;
        this->ngParams.alpha = ngParams.alpha;
        this->ngParams.k = ngParams.k;
        this->ngParams.m = ngParams.m;
        this->ngParams.n = ngParams.n;
        this->ngParams.r0_au = ngParams.r0_au;
        this->isNongrav = true;
    }
    this->isPPN = false;
    this->isMajor = false;
}

void IntegBody::prepare_stm(){
    int stmSize = 36;
    size_t numParams = 0;
    if (this->isNongrav) {
        if (ngParams.a1 != 0.0L) {
            stmSize += 6;
            numParams++;
        }
        if (ngParams.a2 != 0.0L) {
            stmSize += 6;
            numParams++;
        }
        if (ngParams.a3 != 0.0L) {
            stmSize += 6;
            numParams++;
        }
    }
    this->stm = std::vector<real>(stmSize, 0.0L);
    for (size_t i = 0; i < 6; i++) {
        this->stm[6 * i + i] = 1.0L;
    }
    this->n2Derivs += (size_t) stmSize/2;
    this->propStm = true;
    if (this->isCometary){
        std::vector<std::vector<real>> partialsEclip(6, std::vector<real>(6, 0.0L));
        get_cartesian_partials(this->t0, this->initState, "com2cart", partialsEclip);
        std::vector<std::vector<real>> bigRotMat(6, std::vector<real>(6, 0.0L));
        bigRotMat[0][0] = bigRotMat[3][3] = 1.0L;
        bigRotMat[1][1] = bigRotMat[4][4] = bigRotMat[2][2] = bigRotMat[5][5] = cos(EARTH_OBLIQUITY);
        bigRotMat[1][2] = bigRotMat[4][5] = -sin(EARTH_OBLIQUITY);
        bigRotMat[2][1] = bigRotMat[5][4] = sin(EARTH_OBLIQUITY);
        std::vector<std::vector<real>> partials(6, std::vector<real>(6, 0.0L));
        mat_mat_mul(bigRotMat, partialsEclip, partials);
        this->dCartdState = partials;
        for (size_t i = 0; i < 36; i++) {
            this->stm[i] = partials[i/6][i%6];
        }
        for (size_t i = 0; i < 6; i++) {
            for (size_t j = 0; j < numParams; j++) {
                this->dCartdState[i].push_back(0.0L);
            }
        }
    } else {
        this->dCartdState = std::vector<std::vector<real>>(6, std::vector<real>(6+numParams, 0.0L));
        for (size_t i = 0; i < 6; i++) {
            this->dCartdState[i][i] = 1.0L;
        }
    }
    for (size_t i = 0; i < numParams; i++) {
        this->dCartdState.push_back(std::vector<real>(6+numParams, 0.0L));
    }
    for (size_t i = 6; i < 6+numParams; i++) {
        this->dCartdState[i][i] = 1.0L;
    }
}

void ImpulseEvent::apply(const real& t, std::vector<real>& xInteg,
                         const real& propDir) {
    if (t != this->t) {
        throw std::runtime_error(
            "ImpulseEvent::apply: Integration time does "
            "not match event time. Cannot apply impulse.");
    }
    size_t velStartIdx = 6 * this->bodyIndex + 3;
    for (size_t i = 0; i < 3; i++) {
        xInteg[velStartIdx + i] += propDir * this->multiplier * this->deltaV[i];
    }
}

propSimulation::propSimulation(std::string name, real t0,
                               const int defaultSpiceBodies,
                               std::string DEkernelPath) {
    this->name = name;
    this->DEkernelPath = DEkernelPath;
    this->integParams.t0 = t0;

    this->integParams.nInteg = 0;
    this->integParams.nSpice = 0;
    this->integParams.nTotal = 0;
    this->integParams.n2Derivs = 0;
    this->integParams.timestepCounter = 0;

    std::string mapKernelPath =
        DEkernelPath.substr(0, DEkernelPath.find_last_of("/\\")) + "/";
    switch (defaultSpiceBodies) {
        case 0: {
            std::string kernel_sb = mapKernelPath + "sb441-n16s.bsp";
            std::string kernel_mb = mapKernelPath + "de440.bsp";
            spkInfo* mbInfo = spk_init(kernel_mb);
            spkInfo* sbInfo = spk_init(kernel_sb);
            this->ephem.mb = mbInfo;
            this->ephem.sb = sbInfo;
            break;
        }
        // DE430 or DE431
        case 430:
        case 431: {
            std::string kernel_sb = mapKernelPath + "sb431-n16s.bsp";
            std::string kernel_mb = mapKernelPath + "de430.bsp";
            spkInfo* mbInfo = spk_init(kernel_mb);
            spkInfo* sbInfo = spk_init(kernel_sb);
            this->ephem.mb = mbInfo;
            this->ephem.sb = sbInfo;
            real G = 6.6743e-11L /
                (149597870700.0L * 149597870700.0L * 149597870700.0L) *
                86400.0L * 86400.0L;  // default kg au^3 / day^2
            // add planets and planetary bodies from DE431 header
            // (https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de431_tech-comments.txt)
            SpiceBody Sun("Sun", 10, this->integParams.t0,
                          2.959122082855911e-4L / G, 6.96e8L);
            SpiceBody MercuryBarycenter("Mercury Barycenter", 1,
                                        this->integParams.t0,
                                        4.91248045036476e-11L / G, 2440.53e3L);
            SpiceBody VenusBarycenter("Venus Barycenter", 2,
                                      this->integParams.t0,
                                      7.24345233264412e-10L / G, 6051.8e3L);
            SpiceBody Earth("Earth", 399, this->integParams.t0,
                            8.887692445125634e-10L / G, 6378136.3L);
            SpiceBody Moon("Moon", 301, this->integParams.t0,
                           1.093189450742374e-11L / G, 1738.1e3L);
            SpiceBody MarsBarycenter("Mars Barycenter", 4, this->integParams.t0,
                                     9.54954869555077e-11L / G, 3396.19e3L);
            SpiceBody JupiterBarycenter("Jupiter Barycenter", 5,
                                        this->integParams.t0,
                                        2.82534584083387e-07L / G, 0.0L);
            SpiceBody SaturnBarycenter("Saturn Barycenter", 6,
                                       this->integParams.t0,
                                       8.45970607324503e-08L / G, 0.0L);
            SpiceBody UranusBarycenter("Uranus Barycenter", 7,
                                       this->integParams.t0,
                                       1.29202482578296e-08L / G, 0.0L);
            SpiceBody NeptuneBarycenter("Neptune Barycenter", 8,
                                        this->integParams.t0,
                                        1.52435734788511e-08L / G, 0.0L);
            SpiceBody PlutoBarycenter("Pluto Barycenter", 9,
                                      this->integParams.t0,
                                      2.17844105197418e-12L / G, 0.0L);
            Sun.isPPN = true;
            Sun.isMajor = true;
            MercuryBarycenter.isPPN = true;
            MercuryBarycenter.isMajor = true;
            VenusBarycenter.isPPN = true;
            VenusBarycenter.isMajor = true;
            Earth.isPPN = true;
            Earth.isMajor = true;
            Moon.isPPN = true;
            Moon.isMajor = true;
            MarsBarycenter.isPPN = true;
            MarsBarycenter.isMajor = true;
            JupiterBarycenter.isPPN = true;
            JupiterBarycenter.isMajor = true;
            SaturnBarycenter.isPPN = true;
            SaturnBarycenter.isMajor = true;
            UranusBarycenter.isPPN = true;
            UranusBarycenter.isMajor = true;
            NeptuneBarycenter.isPPN = true;
            NeptuneBarycenter.isMajor = true;
            PlutoBarycenter.isPPN = true;
            PlutoBarycenter.isMajor = true;
            Sun.set_J2(
                2.1106088532726840e-7L, 286.13L,
                63.87L);  // https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf
            Earth.set_J2(
                0.00108262545L, 0.0L,
                90.0L);  // https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf
            Sun.caTol = 0.25;
            MercuryBarycenter.caTol = 0.1;
            VenusBarycenter.caTol = 0.1;
            Earth.caTol = 0.1;
            Moon.caTol = 0.05;
            MarsBarycenter.caTol = 0.1;
            JupiterBarycenter.caTol = 0.25;
            SaturnBarycenter.caTol = 0.25;
            UranusBarycenter.caTol = 0.25;
            NeptuneBarycenter.caTol = 0.25;
            PlutoBarycenter.caTol = 0.1;
            add_spice_body(Sun);
            add_spice_body(MercuryBarycenter);
            add_spice_body(VenusBarycenter);
            add_spice_body(Earth);
            add_spice_body(Moon);
            add_spice_body(MarsBarycenter);
            add_spice_body(JupiterBarycenter);
            add_spice_body(SaturnBarycenter);
            add_spice_body(UranusBarycenter);
            add_spice_body(NeptuneBarycenter);
            add_spice_body(PlutoBarycenter);

            // add DE431 big16 asteroids from JPL sb431-big16s.bsp, mass
            // parameters from DE431 header
            // (https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de431_tech-comments.txt)
            SpiceBody Ceres("Ceres", 2000001, this->integParams.t0,
                            1.400476556172344e-13L / G, 0.0L);
            SpiceBody Vesta("Vesta", 2000004, this->integParams.t0,
                            3.85475018780881e-14L / G, 0.0L);
            SpiceBody Pallas("Pallas", 2000002, this->integParams.t0,
                             3.104448198938713e-14L / G, 0.0L);
            SpiceBody Hygiea("Hygiea", 2000010, this->integParams.t0,
                             1.235800787294125e-14L / G, 0.0L);
            SpiceBody Euphrosyne("Euphrosyne", 2000031, this->integParams.t0,
                                 6.343280473648602e-15L / G, 0.0L);
            SpiceBody Interamnia("Interamnia", 2000704, this->integParams.t0,
                                 5.256168678493662e-15L / G, 0.0L);
            SpiceBody Davida("Davida", 2000511, this->integParams.t0,
                             5.198126979457498e-15L / G, 0.0L);
            SpiceBody Eunomia("Eunomia", 2000015, this->integParams.t0,
                              4.678307418350905e-15L / G, 0.0L);
            SpiceBody Juno("Juno", 2000003, this->integParams.t0,
                           3.617538317147937e-15L / G, 0.0L);
            SpiceBody Psyche("Psyche", 2000016, this->integParams.t0,
                             3.411586826193812e-15L / G, 0.0L);
            SpiceBody Cybele("Cybele", 2000065, this->integParams.t0,
                             3.180659282652541e-15L / G, 0.0L);
            SpiceBody Thisbe("Thisbe", 2000088, this->integParams.t0,
                             2.577114127311047e-15L / G, 0.0L);
            SpiceBody Doris("Doris", 2000048, this->integParams.t0,
                            2.531091726015068e-15L / G, 0.0L);
            SpiceBody Europa("Europa", 2000052, this->integParams.t0,
                             2.476788101255867e-15L / G, 0.0L);
            SpiceBody Patientia("Patientia", 2000451, this->integParams.t0,
                                2.295559390637462e-15L / G, 0.0L);
            SpiceBody Sylvia("Sylvia", 2000087, this->integParams.t0,
                             2.199295173574073e-15L / G, 0.0L);
            add_spice_body(Ceres);
            add_spice_body(Vesta);
            add_spice_body(Pallas);
            add_spice_body(Hygiea);
            add_spice_body(Euphrosyne);
            add_spice_body(Interamnia);
            add_spice_body(Davida);
            add_spice_body(Eunomia);
            add_spice_body(Juno);
            add_spice_body(Psyche);
            add_spice_body(Cybele);
            add_spice_body(Thisbe);
            add_spice_body(Doris);
            add_spice_body(Europa);
            add_spice_body(Patientia);
            add_spice_body(Sylvia);
            break;
        }
        // DE440 or DE441
        case 440:
        case 441: {
            std::string kernel_sb = mapKernelPath + "sb441-n16s.bsp";
            std::string kernel_mb = mapKernelPath + "de440.bsp";
            spkInfo* mbInfo = spk_init(kernel_mb);
            spkInfo* sbInfo = spk_init(kernel_sb);
            this->ephem.mb = mbInfo;
            this->ephem.sb = sbInfo;
            real G = 6.6743e-11L /
                (149597870700.0L * 149597870700.0L * 149597870700.0L) *
                86400.0L * 86400.0L;  // default kg au^3 / day^2
            // add planets and planetary bodies from DE441 header
            // (https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_tech-comments.txt)
            SpiceBody Sun("Sun", 10, this->integParams.t0,
                          2.9591220828411956e-04L / G, 6.96e8L);
            SpiceBody MercuryBarycenter(
                "Mercury Barycenter", 1, this->integParams.t0,
                4.9125001948893182e-11L / G, 2440.53e3L);
            SpiceBody VenusBarycenter("Venus Barycenter", 2,
                                      this->integParams.t0,
                                      7.2434523326441187e-10L / G, 6051.8e3L);
            SpiceBody Earth("Earth", 399, this->integParams.t0,
                            8.8876924467071033e-10L / G, 6378136.6L);
            SpiceBody Moon("Moon", 301, this->integParams.t0,
                           1.0931894624024351e-11L / G, 1738.1e3L);
            SpiceBody MarsBarycenter("Mars Barycenter", 4, this->integParams.t0,
                                     9.5495488297258119e-11L / G, 3396.19e3L);
            SpiceBody JupiterBarycenter("Jupiter Barycenter", 5,
                                        this->integParams.t0,
                                        2.8253458252257917e-07L / G, 0.0L);
            SpiceBody SaturnBarycenter("Saturn Barycenter", 6,
                                       this->integParams.t0,
                                       8.4597059933762903e-08L / G, 0.0L);
            SpiceBody UranusBarycenter("Uranus Barycenter", 7,
                                       this->integParams.t0,
                                       1.2920265649682399e-08L / G, 0.0L);
            SpiceBody NeptuneBarycenter("Neptune Barycenter", 8,
                                        this->integParams.t0,
                                        1.5243573478851939e-08L / G, 0.0L);
            SpiceBody PlutoBarycenter("Pluto Barycenter", 9,
                                      this->integParams.t0,
                                      2.1750964648933581e-12L / G, 0.0L);
            Sun.isPPN = true;
            Sun.isMajor = true;
            MercuryBarycenter.isPPN = true;
            MercuryBarycenter.isMajor = true;
            VenusBarycenter.isPPN = true;
            VenusBarycenter.isMajor = true;
            Earth.isPPN = true;
            Earth.isMajor = true;
            Moon.isPPN = true;
            Moon.isMajor = true;
            MarsBarycenter.isPPN = true;
            MarsBarycenter.isMajor = true;
            JupiterBarycenter.isPPN = true;
            JupiterBarycenter.isMajor = true;
            SaturnBarycenter.isPPN = true;
            SaturnBarycenter.isMajor = true;
            UranusBarycenter.isPPN = true;
            UranusBarycenter.isMajor = true;
            NeptuneBarycenter.isPPN = true;
            NeptuneBarycenter.isMajor = true;
            PlutoBarycenter.isPPN = true;
            PlutoBarycenter.isMajor = true;
            Sun.set_J2(
                2.1961391516529825e-07L, 286.13L,
                63.87L);  // https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf
            Earth.set_J2(
                1.0826253900000000e-03L, 0.0L,
                90.0L);  // https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf
            Sun.caTol = 0.25;
            MercuryBarycenter.caTol = 0.1;
            VenusBarycenter.caTol = 0.1;
            Earth.caTol = 0.1;
            Moon.caTol = 0.05;
            MarsBarycenter.caTol = 0.1;
            JupiterBarycenter.caTol = 0.25;
            SaturnBarycenter.caTol = 0.25;
            UranusBarycenter.caTol = 0.25;
            NeptuneBarycenter.caTol = 0.25;
            PlutoBarycenter.caTol = 0.1;
            add_spice_body(Sun);
            add_spice_body(MercuryBarycenter);
            add_spice_body(VenusBarycenter);
            add_spice_body(Earth);
            add_spice_body(Moon);
            add_spice_body(MarsBarycenter);
            add_spice_body(JupiterBarycenter);
            add_spice_body(SaturnBarycenter);
            add_spice_body(UranusBarycenter);
            add_spice_body(NeptuneBarycenter);
            add_spice_body(PlutoBarycenter);

            // add DE441 big16 asteroids from JPL SSD IOM 392R-21-005
            // (ftp://ssd.jpl.nasa.gov/pub/eph/small_bodies/asteroids_de441/SB441_IOM392R-21-005_perturbers.pdf)
            SpiceBody Ceres("Ceres", 2000001, this->integParams.t0,
                            1.3964518123081070e-13L / G, 0.0L);
            SpiceBody Vesta("Vesta", 2000004, this->integParams.t0,
                            3.8548000225257904e-14L / G, 0.0L);
            SpiceBody Pallas("Pallas", 2000002, this->integParams.t0,
                             3.0471146330043200e-14L / G, 0.0L);
            SpiceBody Hygiea("Hygiea", 2000010, this->integParams.t0,
                             1.2542530761640810e-14L / G, 0.0L);
            SpiceBody Davida("Davida", 2000511, this->integParams.t0,
                             8.6836253492286545e-15L / G, 0.0L);
            SpiceBody Interamnia("Interamnia", 2000704, this->integParams.t0,
                                 6.3110343420878887e-15L / G, 0.0L);
            SpiceBody Europa("Europa", 2000052, this->integParams.t0,
                             5.9824315264869841e-15L / G, 0.0L);
            SpiceBody Sylvia("Sylvia", 2000087, this->integParams.t0,
                             4.8345606546105521e-15L / G, 0.0L);
            SpiceBody Eunomia("Eunomia", 2000015, this->integParams.t0,
                              4.5107799051436795e-15L / G, 0.0L);
            SpiceBody Juno("Juno", 2000003, this->integParams.t0,
                           4.2823439677995011e-15L / G, 0.0L);
            SpiceBody Psyche("Psyche", 2000016, this->integParams.t0,
                             3.5445002842488978e-15L / G, 0.0L);
            SpiceBody Camilla("Camilla", 2000107, this->integParams.t0,
                              3.2191392075878588e-15L / G, 0.0L);
            SpiceBody Thisbe("Thisbe", 2000088, this->integParams.t0,
                             2.6529436610356353e-15L / G, 0.0L);
            SpiceBody Iris("Iris", 2000007, this->integParams.t0,
                           2.5416014973471498e-15L / G, 0.0L);
            SpiceBody Euphrosyne("Euphrosyne", 2000031, this->integParams.t0,
                                 2.4067012218937576e-15L / G, 0.0L);
            SpiceBody Cybele("Cybele", 2000065, this->integParams.t0,
                             2.0917175955133682e-15L / G, 0.0L);
            add_spice_body(Ceres);
            add_spice_body(Vesta);
            add_spice_body(Pallas);
            add_spice_body(Hygiea);
            add_spice_body(Davida);
            add_spice_body(Interamnia);
            add_spice_body(Europa);
            add_spice_body(Sylvia);
            add_spice_body(Eunomia);
            add_spice_body(Juno);
            add_spice_body(Psyche);
            add_spice_body(Camilla);
            add_spice_body(Thisbe);
            add_spice_body(Iris);
            add_spice_body(Euphrosyne);
            add_spice_body(Cybele);
            break;
        }
        default:
            throw std::invalid_argument(
                "The defaultSpiceBodies argument is only defined for no "
                "default "
                "SPICE bodies (case 0) or DE431 (case 431) or DE441 (case "
                "441).");
            break;
    }
}

propSimulation::propSimulation(std::string name, const propSimulation& simRef) {
    this->name = name;
    this->DEkernelPath = simRef.DEkernelPath;
    this->ephem = simRef.ephem;
    this->consts = simRef.consts;
    this->integParams = simRef.integParams;
    this->spiceBodies = simRef.spiceBodies;
    // this->integBodies = simRef.integBodies;
    // this->events = simRef.events;
    this->tEvalUTC = simRef.tEvalUTC;
    this->evalApparentState = simRef.evalApparentState;
    this->evalMeasurements = simRef.evalMeasurements;
    this->convergedLightTime = simRef.convergedLightTime;
    this->observerInfo = simRef.observerInfo;
    this->xObserver = simRef.xObserver;
    this->tEvalMargin = simRef.tEvalMargin;
    this->tEval = simRef.tEval;
    this->radarObserver = simRef.radarObserver;
}

void propSimulation::prepare_for_evaluation(
    std::vector<real>& tEval, std::vector<std::vector<real>>& observerInfo) {
    const bool forwardProp = this->integParams.t0 <= this->integParams.tf;
    const bool backwardProp = this->integParams.t0 >= this->integParams.tf;
    if (forwardProp && backwardProp) {
        throw std::invalid_argument(
            "The initial and final times must be different.");
    }
    // sort tEval into ascending order or descending order based on the
    // integration direction (not during orbit fits)
    if (observerInfo.size() == 0) {
    std::vector<size_t> tEvalSortedIdx(tEval.size());
        sort_vector(tEval, forwardProp, tEvalSortedIdx);
    }
    if (forwardProp) {
        int removeCounter = 0;
        while (tEval[0] < this->integParams.t0 - this->tEvalMargin) {
            // remove any tEval values that are before the integration start
            // time
            tEval.erase(tEval.begin());
            if (observerInfo.size() != 0)
                observerInfo.erase(observerInfo.begin());
            removeCounter++;
        }
        while (tEval.back() > this->integParams.tf + this->tEvalMargin) {
            // remove any tEval values that are after the integration end time
            tEval.pop_back();
            if (observerInfo.size() != 0) observerInfo.pop_back();
            removeCounter++;
        }
        if (removeCounter > 0) {
            std::cout << "WARNING: " << removeCounter
                      << " tEval and observerInfo value(s) were removed "
                         "because they were outside the interpolation range, "
                         "i.e., integration range with a margin of "
                      << this->tEvalMargin << " day(s)." << std::endl;
        }
    } else if (backwardProp) {
        int removeCounter = 0;
        while (tEval[0] > this->integParams.t0 + this->tEvalMargin) {
            // remove any tEval values that are after the integration start time
            tEval.erase(tEval.begin());
            if (observerInfo.size() != 0)
                observerInfo.erase(observerInfo.begin());
            removeCounter++;
        }
        while (tEval.back() < this->integParams.tf - this->tEvalMargin) {
            // remove any tEval values that are before the integration end time
            tEval.pop_back();
            if (observerInfo.size() != 0) observerInfo.pop_back();
            removeCounter++;
        }
        if (removeCounter > 0) {
            std::cout << "WARNING: " << removeCounter
                      << " tEval and observerInfo value(s) were removed "
                         "because they were outside the interpolation range, "
                         "i.e., integration range with a margin of "
                      << this->tEvalMargin << " day(s)." << std::endl;
        }
    }

    if (observerInfo.size() == 0) {
        observerInfo = std::vector<std::vector<real>>(
            tEval.size(), std::vector<real>(4, 0.0L));
    }

    if (this->observerInfo.size() == 0) {
        this->observerInfo = observerInfo;
    } else if (this->observerInfo.size() != 0) {
        for (size_t i = 0; i < observerInfo.size(); i++) {
            this->observerInfo.push_back(observerInfo[i]);
        }
    }

    std::vector<std::vector<real>> xObserver = std::vector<std::vector<real>>(
        tEval.size(), std::vector<real>(6, 0.0L));
    std::vector<int> radarObserver = std::vector<int>(tEval.size(), 0);
    if (this->tEvalUTC) {
        furnsh_c(this->DEkernelPath.c_str());
    }
    if (tEval.size() != 0) {
        for (size_t i = 0; i < tEval.size(); i++) {
            if (observerInfo[i].size() == 4 || observerInfo[i].size() == 7) {
                radarObserver[i] = 0;
            } else if (observerInfo[i].size() == 9) {
                radarObserver[i] = 1;
            } else if (observerInfo[i].size() == 10) {
                radarObserver[i] = 2;
            } else {
                throw std::invalid_argument(
                    "The observerInfo vector must have 4/7 (optical), 9 (radar "
                    "delay), or 10 elements (radar doppler).");
            }
            get_observer_state(tEval[i], observerInfo[i], this,
                               this->tEvalUTC, xObserver[i]);
            // std::cout << xObserver[i][0] << " " << xObserver[i][1] << " " <<
            // xObserver[i][2] << " " << xObserver[i][3] << " " <<
            // xObserver[i][4] << " " << xObserver[i][5] << std::endl;
        }
    }
    if (this->tEvalUTC) {
        unload_c(this->DEkernelPath.c_str());
    }

    if (this->tEval.size() == 0) {
        this->tEval = tEval;
        this->xObserver = xObserver;
        this->radarObserver = radarObserver;
    } else if (this->tEval.size() != 0) {
        for (size_t i = 0; i < tEval.size(); i++) {
            this->tEval.push_back(tEval[i]);
            this->xObserver.push_back(xObserver[i]);
            this->radarObserver.push_back(radarObserver[i]);
        }
    }
}

void propSimulation::add_spice_body(SpiceBody body) {
    // check if body already exists. if so, throw error
    for (size_t i = 0; i < this->spiceBodies.size(); i++) {
        if (this->spiceBodies[i].name == body.name) {
            throw std::invalid_argument("SPICE Body with name " + body.name +
                                        " already exists in simulation " +
                                        this->name);
        }
    }
    body.radius /= this->consts.du2m;
    this->spiceBodies.push_back(body);
    this->integParams.nSpice++;
    this->integParams.nTotal++;
}

std::vector<real> propSimulation::get_spiceBody_state(const real t, const std::string &bodyName) {
    int spiceID = -1;
    for (size_t i = 0; i < this->spiceBodies.size(); i++){
        if (this->spiceBodies[i].name == bodyName){
            spiceID = this->spiceBodies[i].spiceId;
            break;
        }
    }
    if (spiceID == -1){
        throw std::invalid_argument("SPICE Body with name " + bodyName +
                                        " does not exist in simulation " +
                                        this->name);
    }
    double spiceState[9];
    get_spk_state(spiceID, t, this->ephem, spiceState);
    std::vector<real> state = {spiceState[0], spiceState[1], spiceState[2],
                               spiceState[3], spiceState[4], spiceState[5]};
    return state;
}

void propSimulation::add_integ_body(IntegBody body) {
    // check if body already exists. if so, throw error
    for (size_t i = 0; i < this->integBodies.size(); i++) {
        if (this->integBodies[i].name == body.name) {
            throw std::invalid_argument(
                "Integration body with name " + body.name +
                " already exists in simulation " + this->name);
        }
    }
    if (body.t0 != this->integParams.t0) {
        throw std::invalid_argument(
            "Integration body " + body.name + " has initial time MJD " +
            std::to_string(body.t0) +
            " TDB which is different from the simulation initial time: MJD " +
            std::to_string(this->integParams.t0) + " TDB.");
    }
    if (body.isCometary) {
        double sunState[9];
        get_spk_state(10, body.t0, this->ephem, sunState);
        body.pos[0] += sunState[0];
        body.pos[1] += sunState[1];
        body.pos[2] += sunState[2];
        body.vel[0] += sunState[3];
        body.vel[1] += sunState[4];
        body.vel[2] += sunState[5];
    }
    body.initState = {body.pos[0], body.pos[1], body.pos[2],
                      body.vel[0], body.vel[1], body.vel[2]};
    body.radius /= this->consts.du2m;
    this->integBodies.push_back(body);
    this->integParams.nInteg++;
    this->integParams.nTotal++;
    this->integParams.n2Derivs += body.n2Derivs;
}

void propSimulation::remove_body(std::string name) {
    for (size_t i = 0; i < this->spiceBodies.size(); i++) {
        if (this->spiceBodies[i].name == name) {
            this->spiceBodies.erase(this->spiceBodies.begin() + i);
            this->integParams.nSpice--;
            this->integParams.nTotal--;
            return;
        }
    }
    for (size_t i = 0; i < this->integBodies.size(); i++) {
        if (this->integBodies[i].name == name) {
            this->integBodies.erase(this->integBodies.begin() + i);
            this->integParams.nInteg--;
            this->integParams.nTotal--;
            return;
        }
    }
    std::cout << "Error: Body " << name << " not found." << std::endl;
}

void propSimulation::add_event(IntegBody body, real tEvent,
                               std::vector<real> deltaV, real multiplier) {
    // check if tEvent is valid
    const bool forwardProp = this->integParams.tf > this->integParams.t0;
    const bool backwardProp = this->integParams.tf < this->integParams.t0;
    if ((forwardProp &&
         (tEvent < this->integParams.t0 || tEvent >= this->integParams.tf)) ||
        (backwardProp &&
         (tEvent > this->integParams.t0 || tEvent <= this->integParams.tf))) {
        throw std::invalid_argument("Event time " + std::to_string(tEvent) +
                                    " is not within simulation time bounds.");
    }
    // check if body exists
    bool bodyExists = false;
    size_t bodyIndex;
    for (size_t i = 0; i < this->integParams.nInteg; i++) {
        if (this->integBodies[i].name == body.name) {
            bodyExists = true;
            bodyIndex = i;
            break;
        }
    }
    if (!bodyExists) {
        throw std::invalid_argument("Integration body with name " + body.name +
                                    " does not exist in simulation " +
                                    this->name);
    }
    ImpulseEvent event;
    event.t = tEvent;
    event.deltaV = deltaV;
    event.multiplier = multiplier;
    event.bodyName = body.name;
    event.bodyIndex = bodyIndex;
    // add event to this->events sorted by time
    if (this->events.size() == 0) {
        this->events.push_back(event);
    } else {
        for (size_t i = 0; i < this->events.size(); i++) {
            if (event.t < this->events[i].t) {
                this->events.insert(this->events.begin() + i, event);
                break;
            } else if (i == this->events.size() - 1) {
                this->events.push_back(event);
                break;
            }
        }
    }
}

void propSimulation::set_sim_constants(real du2m, real tu2s, real G,
                                       real clight) {
    this->consts.du2m = du2m;
    this->consts.tu2s = tu2s;
    this->consts.duptu2mps = du2m / tu2s;
    this->consts.G = G;
    this->consts.clight = clight;
    this->consts.j2000Jd = 2451545.0;
    this->consts.JdMinusMjd = 2400000.5;
}

void propSimulation::set_integration_parameters(
    real tf, std::vector<real> tEval, bool tEvalUTC, bool evalApparentState,
    bool convergedLightTime, std::vector<std::vector<real>> observerInfo,
    bool adaptiveTimestep, real dt0, real dtMax, real dtMin,
    real dtChangeFactor, real tolInteg, real tolPC) {
    this->integParams.tf = tf;
    this->tEvalUTC = tEvalUTC;
    this->evalApparentState = evalApparentState;
    this->convergedLightTime = convergedLightTime;
    if (tEval.size() != 0) {
        prepare_for_evaluation(tEval, observerInfo);
    }
    this->integParams.dt0 = dt0;
    this->integParams.dtMax = dtMax;
    this->integParams.dtMin = dtMin;
    this->integParams.dtChangeFactor = dtChangeFactor;
    this->integParams.adaptiveTimestep = adaptiveTimestep;
    this->integParams.tolPC = tolPC;
    this->integParams.tolInteg = tolInteg;
}

std::vector<real> propSimulation::get_sim_constants() {
    std::vector<real> constants = {
        this->consts.du2m,      this->consts.tu2s,   this->consts.duptu2mps,
        this->consts.G,         this->consts.clight, this->consts.j2000Jd,
        this->consts.JdMinusMjd};
    return constants;
}

std::vector<real> propSimulation::get_integration_parameters() {
    std::vector<real> integration_parameters = {
        (real)this->integParams.nInteg,
        (real)this->integParams.nSpice,
        (real)this->integParams.nTotal,
        this->integParams.t0,
        this->integParams.tf,
        (real)this->integParams.adaptiveTimestep,
        this->integParams.dt0,
        this->integParams.dtMax,
        this->integParams.dtMin,
        this->integParams.dtChangeFactor,
        this->integParams.tolInteg,
        this->integParams.tolPC};
    return integration_parameters;
}

void propSimulation::preprocess() {
    if (!this->isPreprocessed) {
        this->t = this->integParams.t0;
        for (size_t i = 0; i < this->integParams.nInteg; i++) {
            for (size_t j = 0; j < 3; j++) {
                this->xInteg.push_back(integBodies[i].pos[j]);
            }
            for (size_t j = 0; j < 3; j++) {
                this->xInteg.push_back(integBodies[i].vel[j]);
            }
            if (integBodies[i].propStm) {
                for (size_t j = 0; j < integBodies[i].stm.size(); j++) {
                    this->xInteg.push_back(integBodies[i].stm[j]);
                }
            }
        }
        this->interpParams.tStack.push_back(t);
        this->interpParams.xIntegStack.push_back(xInteg);
        bool backwardProp = this->integParams.t0 > this->integParams.tf;
        if (backwardProp) {
            std::reverse(this->events.begin(), this->events.end());
        }
        this->isPreprocessed = true;
    }
}

void propSimulation::extend(real tf, std::vector<real> tEvalNew,
                            std::vector<std::vector<real>> xObserverNew) {
    std::cout << "WARNING: The extend() method is under development and may "
                 "not work properly with the interpolation routines."
              << std::endl;

    // empty existing vectors from previous integration
    this->caParams.clear();
    this->interpParams.tStack.clear();
    this->interpParams.xIntegStack.clear();
    this->interpParams.bStack.clear();
    this->interpParams.accIntegStack.clear();
    this->interpIdx = 0;
    this->xObserver.clear();
    this->observerInfo.clear();
    this->tEval.clear();
    this->radarObserver.clear();
    this->lightTimeEval.clear();
    this->xIntegEval.clear();
    this->opticalObs.clear();
    this->opticalPartials.clear();
    this->radarObs.clear();
    this->radarPartials.clear();

    this->integParams.t0 = this->t;
    this->interpParams.tStack.push_back(this->t);
    this->interpParams.xIntegStack.push_back(this->xInteg);
    this->set_integration_parameters(tf, tEvalNew, this->tEvalUTC,
                                     this->evalApparentState,
                                     this->convergedLightTime, xObserverNew);
    this->integrate();
}
