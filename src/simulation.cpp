#include "simulation.h"

void Body::set_J2(real J2, real poleRA, real poleDec) {
    this->J2 = J2;
    if (this->J2 != 0.0L) {
        this->isJ2 = true;
    } else {
        this->isJ2 = false;
    }
    this->poleRA = poleRA*DEG2RAD;
    this->poleDec = poleDec*DEG2RAD;
}

SpiceBody::SpiceBody(std::string DEkernelPath, std::string name, int spiceId,
                     real t0, real mass, real radius, Constants consts) {
    this->name = std::to_string(spiceId) + " " + name;
    this->spiceId = spiceId;
    this->t0 = t0;
    this->mass = mass;
    this->radius = radius / consts.du2m;
    this->isNongrav = false;
    this->isPPN = false;
    this->isMajor = false;
    if (this->isSpice) {
        double state[6];
        double lt;
        furnsh_c(DEkernelPath.c_str());
        get_spice_state_lt(this->spiceId, this->t0, consts, state, lt);
        unload_c(DEkernelPath.c_str());
        this->pos = {state[0], state[1], state[2]};
        this->vel = {state[3], state[4], state[5]};
    }
}

IntegBody::IntegBody(std::string DEkernelPath, std::string name, real t0,
                     real mass, real radius, std::vector<real> cometaryState,
                     std::vector<std::vector<real>> covariance,
                     NongravParamaters ngParams, Constants consts) {
    this->name = name;
    this->t0 = t0;
    this->mass = mass;
    this->radius = radius / consts.du2m;
    std::vector<real> cartesianStateEclip(6);
    std::vector<real> cartesianPos(3);
    std::vector<real> cartesianVel(3);

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
    // shift heliocentric to barycentric
    double sunState[6];
    double lt;
    furnsh_c(DEkernelPath.c_str());
    get_spice_state_lt(10, t0, consts, sunState, lt);
    unload_c(DEkernelPath.c_str());
    for (size_t i = 0; i < 3; i++) {
        cartesianPos[i] += sunState[i];
        cartesianVel[i] += sunState[i + 3];
    }
    this->pos = cartesianPos;
    this->vel = cartesianVel;
    this->covariance = covariance;
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
                     std::vector<std::vector<real>> covariance,
                     NongravParamaters ngParams, Constants consts) {
    this->name = name;
    this->t0 = t0;
    this->mass = mass;
    this->radius = radius / consts.du2m;
    this->pos = pos;
    this->vel = vel;
    this->covariance = covariance;
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

void ImpulseEvent::apply(const real &t, std::vector<real> &xInteg,
                         const real &propDir) {
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
    this->integParams.timestepCounter = 0;

    switch (defaultSpiceBodies) {
        case 0: {
            break;
        }
        case 431: {
            real G = 6.6743e-11L /
                (149597870700.0L * 149597870700.0L * 149597870700.0L) *
                86400.0L * 86400.0L;  // default kg au^3 / day^2
            // add planets and planetary bodies from DE431 header
            // (https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de431_tech-comments.txt)
            SpiceBody Sun(DEkernelPath, "Sun", 10, this->integParams.t0,
                          2.959122082855911e-4L / G, 6.96e8L, this->consts);
            SpiceBody MercuryBarycenter(
                DEkernelPath, "Mercury Barycenter", 1, this->integParams.t0,
                4.91248045036476e-11L / G, 0.0L, this->consts);
            SpiceBody VenusBarycenter(
                DEkernelPath, "Venus Barycenter", 2, this->integParams.t0,
                7.24345233264412e-10L / G, 0.0L, this->consts);
            SpiceBody Earth(DEkernelPath, "Earth", 399, this->integParams.t0,
                            8.887692445125634e-10L / G, 6378136.3L,
                            this->consts);
            SpiceBody Moon(DEkernelPath, "Moon", 301, this->integParams.t0,
                           1.093189450742374e-11L / G, 0.0L, this->consts);
            SpiceBody MarsBarycenter(
                DEkernelPath, "Mars Barycenter", 4, this->integParams.t0,
                9.54954869555077e-11L / G, 0.0L, this->consts);
            SpiceBody JupiterBarycenter(
                DEkernelPath, "Jupiter Barycenter", 5, this->integParams.t0,
                2.82534584083387e-07L / G, 0.0L, this->consts);
            SpiceBody SaturnBarycenter(
                DEkernelPath, "Saturn Barycenter", 6, this->integParams.t0,
                8.45970607324503e-08L / G, 0.0L, this->consts);
            SpiceBody UranusBarycenter(
                DEkernelPath, "Uranus Barycenter", 7, this->integParams.t0,
                1.29202482578296e-08L / G, 0.0L, this->consts);
            SpiceBody NeptuneBarycenter(
                DEkernelPath, "Neptune Barycenter", 8, this->integParams.t0,
                1.52435734788511e-08L / G, 0.0L, this->consts);
            SpiceBody PlutoBarycenter(
                DEkernelPath, "Pluto Barycenter", 9, this->integParams.t0,
                2.17844105197418e-12L / G, 0.0L, this->consts);
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
                2.1106088532726840e-7L,
                286.13L, 63.87L);  // https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf
            Earth.set_J2(
                0.00108262545L,
                0.0L, 90.0L);  // https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf
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
            SpiceBody Ceres(DEkernelPath, "Ceres", 2000001,
                            this->integParams.t0, 1.400476556172344e-13L / G,
                            0.0L, this->consts);
            SpiceBody Vesta(DEkernelPath, "Vesta", 2000004,
                            this->integParams.t0, 3.85475018780881e-14L / G,
                            0.0L, this->consts);
            SpiceBody Pallas(DEkernelPath, "Pallas", 2000002,
                             this->integParams.t0, 3.104448198938713e-14L / G,
                             0.0L, this->consts);
            SpiceBody Hygiea(DEkernelPath, "Hygiea", 2000010,
                             this->integParams.t0, 1.235800787294125e-14L / G,
                             0.0L, this->consts);
            SpiceBody Euphrosyne(
                DEkernelPath, "Euphrosyne", 2000031, this->integParams.t0,
                6.343280473648602e-15L / G, 0.0L, this->consts);
            SpiceBody Interamnia(
                DEkernelPath, "Interamnia", 2000704, this->integParams.t0,
                5.256168678493662e-15L / G, 0.0L, this->consts);
            SpiceBody Davida(DEkernelPath, "Davida", 2000511,
                             this->integParams.t0, 5.198126979457498e-15L / G,
                             0.0L, this->consts);
            SpiceBody Eunomia(DEkernelPath, "Eunomia", 2000015,
                              this->integParams.t0, 4.678307418350905e-15L / G,
                              0.0L, this->consts);
            SpiceBody Juno(DEkernelPath, "Juno", 2000003, this->integParams.t0,
                           3.617538317147937e-15L / G, 0.0L, this->consts);
            SpiceBody Psyche(DEkernelPath, "Psyche", 2000016,
                             this->integParams.t0, 3.411586826193812e-15L / G,
                             0.0L, this->consts);
            SpiceBody Cybele(DEkernelPath, "Cybele", 2000065,
                             this->integParams.t0, 3.180659282652541e-15L / G,
                             0.0L, this->consts);
            SpiceBody Thisbe(DEkernelPath, "Thisbe", 2000088,
                             this->integParams.t0, 2.577114127311047e-15L / G,
                             0.0L, this->consts);
            SpiceBody Doris(DEkernelPath, "Doris", 2000048,
                            this->integParams.t0, 2.531091726015068e-15L / G,
                            0.0L, this->consts);
            SpiceBody Europa(DEkernelPath, "Europa", 2000052,
                             this->integParams.t0, 2.476788101255867e-15L / G,
                             0.0L, this->consts);
            SpiceBody Patientia(DEkernelPath, "Patientia", 2000451,
                                this->integParams.t0,
                                2.295559390637462e-15L / G, 0.0L, this->consts);
            SpiceBody Sylvia(DEkernelPath, "Sylvia", 2000087,
                             this->integParams.t0, 2.199295173574073e-15L / G,
                             0.0L, this->consts);
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
        case 441: {
            real G = 6.6743e-11L /
                (149597870700.0L * 149597870700.0L * 149597870700.0L) *
                86400.0L * 86400.0L;  // default kg au^3 / day^2
            // add planets and planetary bodies from DE441 header
            // (https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_tech-comments.txt)
            SpiceBody Sun(DEkernelPath, "Sun", 10, this->integParams.t0,
                          2.9591220828411956e-04L / G, 6.96e8L, this->consts);
            SpiceBody MercuryBarycenter(
                DEkernelPath, "Mercury Barycenter", 1, this->integParams.t0,
                4.9125001948893182e-11L / G, 0.0L, this->consts);
            SpiceBody VenusBarycenter(
                DEkernelPath, "Venus Barycenter", 2, this->integParams.t0,
                7.2434523326441187e-10L / G, 0.0L, this->consts);
            SpiceBody Earth(DEkernelPath, "Earth", 399, this->integParams.t0,
                            8.8876924467071033e-10L / G, 6378136.6L,
                            this->consts);
            SpiceBody Moon(DEkernelPath, "Moon", 301, this->integParams.t0,
                           1.0931894624024351e-11L / G, 0.0L, this->consts);
            SpiceBody MarsBarycenter(
                DEkernelPath, "Mars Barycenter", 4, this->integParams.t0,
                9.5495488297258119e-11L / G, 0.0L, this->consts);
            SpiceBody JupiterBarycenter(
                DEkernelPath, "Jupiter Barycenter", 5, this->integParams.t0,
                2.8253458252257917e-07L / G, 0.0L, this->consts);
            SpiceBody SaturnBarycenter(
                DEkernelPath, "Saturn Barycenter", 6, this->integParams.t0,
                8.4597059933762903e-08L / G, 0.0L, this->consts);
            SpiceBody UranusBarycenter(
                DEkernelPath, "Uranus Barycenter", 7, this->integParams.t0,
                1.2920265649682399e-08L / G, 0.0L, this->consts);
            SpiceBody NeptuneBarycenter(
                DEkernelPath, "Neptune Barycenter", 8, this->integParams.t0,
                1.5243573478851939e-08L / G, 0.0L, this->consts);
            SpiceBody PlutoBarycenter(
                DEkernelPath, "Pluto Barycenter", 9, this->integParams.t0,
                2.1750964648933581e-12L / G, 0.0L, this->consts);
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
                2.1961391516529825e-07L,
                286.13L, 63.87L);  // https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf
            // MercuryBarycenter.set_J2(50.3e-6L, 0.034L*DEG2RAD); //
            // https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html
            // VenusBarycenter.set_J2(4.458e-6L, 177.36L*DEG2RAD); //
            // https://nssdc.gsfc.nasa.gov/planetary/factsheet/venusfact.html
            Earth.set_J2(
                1.0826253900000000e-03L,
                0.0L, 90.0L);  // https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf
            // Moon.set_J2(2.0321568464952570e-4L, 5.145L*DEG2RAD); //
            // https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
            // MarsBarycenter.set_J2(1960.45e-6L, 25.19L*DEG2RAD); //
            // https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
            // JupiterBarycenter.set_J2(14736.0e-6L, 3.13L*DEG2RAD); //
            // https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
            // SaturnBarycenter.set_J2(16298.0e-6L, 26.73L*DEG2RAD); //
            // https://nssdc.gsfc.nasa.gov/planetary/factsheet/saturnfact.html
            // UranusBarycenter.set_J2(3343.43e-6L, 97.77L*DEG2RAD); //
            // https://nssdc.gsfc.nasa.gov/planetary/factsheet/uranusfact.html
            // NeptuneBarycenter.set_J2(3411.0e-6L, 28.32L*DEG2RAD); //
            // https://nssdc.gsfc.nasa.gov/planetary/factsheet/neptunefact.html
            // PlutoBarycenter.set_J2(0.0L, 122.53L); //
            // https://nssdc.gsfc.nasa.gov/planetary/factsheet/plutofact.html
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
            SpiceBody Ceres(DEkernelPath, "Ceres", 2000001,
                            this->integParams.t0, 1.3964518123081070e-13L / G,
                            0.0L, this->consts);
            SpiceBody Vesta(DEkernelPath, "Vesta", 2000004,
                            this->integParams.t0, 3.8548000225257904e-14L / G,
                            0.0L, this->consts);
            SpiceBody Pallas(DEkernelPath, "Pallas", 2000002,
                             this->integParams.t0, 3.0471146330043200e-14L / G,
                             0.0L, this->consts);
            SpiceBody Hygiea(DEkernelPath, "Hygiea", 2000010,
                             this->integParams.t0, 1.2542530761640810e-14L / G,
                             0.0L, this->consts);
            SpiceBody Davida(DEkernelPath, "Davida", 2000511,
                             this->integParams.t0, 8.6836253492286545e-15L / G,
                             0.0L, this->consts);
            SpiceBody Interamnia(
                DEkernelPath, "Interamnia", 2000704, this->integParams.t0,
                6.3110343420878887e-15L / G, 0.0L, this->consts);
            SpiceBody Europa(DEkernelPath, "Europa", 2000052,
                             this->integParams.t0, 5.9824315264869841e-15L / G,
                             0.0L, this->consts);
            SpiceBody Sylvia(DEkernelPath, "Sylvia", 2000087,
                             this->integParams.t0, 4.8345606546105521e-15L / G,
                             0.0L, this->consts);
            SpiceBody Eunomia(DEkernelPath, "Eunomia", 2000015,
                              this->integParams.t0, 4.5107799051436795e-15L / G,
                              0.0L, this->consts);
            SpiceBody Juno(DEkernelPath, "Juno", 2000003, this->integParams.t0,
                           4.2823439677995011e-15L / G, 0.0L, this->consts);
            SpiceBody Psyche(DEkernelPath, "Psyche", 2000016,
                             this->integParams.t0, 3.5445002842488978e-15L / G,
                             0.0L, this->consts);
            SpiceBody Camilla(DEkernelPath, "Camilla", 2000107,
                              this->integParams.t0, 3.2191392075878588e-15L / G,
                              0.0L, this->consts);
            SpiceBody Thisbe(DEkernelPath, "Thisbe", 2000088,
                             this->integParams.t0, 2.6529436610356353e-15L / G,
                             0.0L, this->consts);
            SpiceBody Iris(DEkernelPath, "Iris", 2000007, this->integParams.t0,
                           2.5416014973471498e-15L / G, 0.0L, this->consts);
            SpiceBody Euphrosyne(
                DEkernelPath, "Euphrosyne", 2000031, this->integParams.t0,
                2.4067012218937576e-15L / G, 0.0L, this->consts);
            SpiceBody Cybele(DEkernelPath, "Cybele", 2000065,
                             this->integParams.t0, 2.0917175955133682e-15L / G,
                             0.0L, this->consts);
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

propSimulation::propSimulation(std::string name, const propSimulation &simRef) {
    this->name = name;
    this->DEkernelPath = simRef.DEkernelPath;
    this->integParams = simRef.integParams;
    this->consts = simRef.consts;
    this->integBodies = simRef.integBodies;
    this->spiceBodies = simRef.spiceBodies;
}

void propSimulation::prepare_for_evaluation(
    std::vector<real> &tEval, std::vector<std::vector<real>> &observerInfo) {
    bool forwardProp = this->integParams.t0 <= this->integParams.tf;
    bool backwardProp = this->integParams.t0 >= this->integParams.tf;
    if (forwardProp && backwardProp) {
        throw std::invalid_argument(
            "The initial and final times must be different.");
    }
    // sort tEval and corresponding observerInfo into ascending order or
    // descending order based on the integration direction
    sort_vector(tEval, forwardProp);
    if (observerInfo.size() != 0) {
        if (observerInfo.size() != tEval.size()) {
            throw std::invalid_argument(
                "The number of tEval values and the number of observerInfo "
                "vectors must be equal.");
        }
        sort_vector_by_another(observerInfo, tEval, forwardProp);
        if (backwardProp) {
            std::reverse(observerInfo.begin(), observerInfo.end());
        }
    }
    if (forwardProp) {
        int removeCounter = 0;
        while (tEval[0] < this->integParams.t0 - this->tEvalMargin) {
            // remove any tEval values that are before the integration start
            // time
            tEval.erase(tEval.begin());
            if (observerInfo.size() != 0) observerInfo.erase(observerInfo.begin());
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
            if (observerInfo.size() != 0) observerInfo.erase(observerInfo.begin());
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
    furnsh_c(this->DEkernelPath.c_str());
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
            get_observer_state(tEval[i], observerInfo[i], this->consts,
                               this->tEvalUTC, xObserver[i]);
            // std::cout << xObserver[i][0] << " " << xObserver[i][1] << " " <<
            // xObserver[i][2] << " " << xObserver[i][3] << " " <<
            // xObserver[i][4] << " " << xObserver[i][5] << std::endl;
        }
    }
    unload_c(this->DEkernelPath.c_str());

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

void propSimulation::add_spice_body(std::string DEkernelPath, std::string name,
                                    int spiceId, real t0, real mass,
                                    real radius, Constants consts) {
    // check if body already exists. if so, throw error
    for (size_t i = 0; i < this->spiceBodies.size(); i++) {
        if (this->spiceBodies[i].name == name) {
            throw std::invalid_argument("SPICE Body with name " + name +
                                        " already exists in simulation " +
                                        this->name);
        }
    }
    SpiceBody body(DEkernelPath, name, spiceId, t0, mass, radius, consts);
    this->spiceBodies.push_back(body);
    this->integParams.nSpice++;
    this->integParams.nTotal++;
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
    this->spiceBodies.push_back(body);
    this->integParams.nSpice++;
    this->integParams.nTotal++;
}

void propSimulation::add_integ_body(std::string DEkernelPath, std::string name,
                                    real t0, real mass, real radius,
                                    std::vector<real> cometaryState,
                                    std::vector<std::vector<real>> covariance,
                                    NongravParamaters ngParams,
                                    Constants consts) {
    // check if body already exists. if so, throw error
    for (size_t i = 0; i < this->integBodies.size(); i++) {
        if (this->integBodies[i].name == name) {
            throw std::invalid_argument("Integration body with name " + name +
                                        " already exists in simulation " +
                                        this->name);
        }
    }
    IntegBody body(DEkernelPath, name, t0, mass, radius, cometaryState,
                   covariance, ngParams, consts);
    this->integBodies.push_back(body);
    this->integParams.nInteg++;
    this->integParams.nTotal++;
}

void propSimulation::add_integ_body(std::string name, real t0, real mass,
                                    real radius, std::vector<real> pos,
                                    std::vector<real> vel,
                                    std::vector<std::vector<real>> covariance,
                                    NongravParamaters ngParams,
                                    Constants consts) {
    // check if body already exists. if so, throw error
    for (size_t i = 0; i < this->integBodies.size(); i++) {
        if (this->integBodies[i].name == name) {
            throw std::invalid_argument("Integration body with name " + name +
                                        " already exists in simulation " +
                                        this->name);
        }
    }
    IntegBody body(name, t0, mass, radius, pos, vel, covariance, ngParams,
                   consts);
    this->integBodies.push_back(body);
    this->integParams.nInteg++;
    this->integParams.nTotal++;
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
    this->integBodies.push_back(body);
    this->integParams.nInteg++;
    this->integParams.nTotal++;
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
    bool forwardProp = this->integParams.tf > this->integParams.t0;
    bool backwardProp = this->integParams.tf < this->integParams.t0;
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

void propSimulation::set_sim_constants(real du2m, real tu2sec, real G,
                                       real clight) {
    this->consts.du2m = du2m;
    this->consts.tu2sec = tu2sec;
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

    bool backwardProp = this->integParams.t0 > this->integParams.tf;
    if (backwardProp) {
        std::reverse(this->events.begin(), this->events.end());
    }
}

std::vector<real> propSimulation::get_sim_constants() {
    // std::cout << "The conversion from distance units to meters is: " <<
    // consts.du2m << std::endl; std::cout << "The conversion from time units to
    // seconds is: " << consts.tu2sec << std::endl; std::cout << "The universal
    // gravitational constant in nondimensional units is: " << consts.G <<
    // std::endl; std::cout << "The speed of light in nondimensional units is: "
    // << consts.clight << std::endl; std::cout << "The Julian date of J2000 is:
    // " << consts.j2000Jd << std::endl; std::cout << "The difference between
    // the Julian date and the Modified Julian date is: " << consts.JdMinusMjd
    // << std::endl;

    std::vector<real> constants = {
        this->consts.du2m,   this->consts.tu2sec,  this->consts.G,
        this->consts.clight, this->consts.j2000Jd, this->consts.JdMinusMjd};
    return constants;
}

std::vector<real> propSimulation::get_integration_parameters() {
    // std::cout << "The number of integrated bodies is: " << integParams.nInteg
    // << std::endl; std::cout << "The number of bodies whose states are taken
    // from SPICE kernels is: " << integParams.nSpice << std::endl; std::cout <<
    // "The total number of bodies in this simulation is: " <<
    // integParams.nTotal << std::endl; std::cout << "The initial time of the
    // simulation is: " << integParams.t0 << std::endl; std::cout << "The final
    // time of the simulation is: " << integParams.tf << std::endl; std::cout <<
    // "The initial timestep of the simulation is: " << integParams.dt0 <<
    // std::endl; std::cout << "The maximum timestep of the simulation is: " <<
    // integParams.dtMax << std::endl; std::cout << "The factor by which the
    // timestep is changed is: " << integParams.dtChangeFactor << std::endl;
    // std::cout << "The adaptive timestep flag is: " <<
    // integParams.adaptiveTimestep << std::endl; std::cout << "The tolerance
    // for the position and velocity error is: " << integParams.tolPC <<
    // std::endl; std::cout << "The tolerance for the integration error is: " <<
    // integParams.tolInteg << std::endl;

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
    this->t = this->integParams.t0;
    for (size_t i = 0; i < this->integParams.nInteg; i++) {
        for (size_t j = 0; j < 3; j++) {
            this->xInteg.push_back(integBodies[i].pos[j]);
        }
        for (size_t j = 0; j < 3; j++) {
            this->xInteg.push_back(integBodies[i].vel[j]);
        }
        this->forceParams.masses.push_back(integBodies[i].mass);
        this->forceParams.radii.push_back(integBodies[i].radius);
        this->forceParams.spiceIdList.push_back(-99999);
        this->forceParams.isPPNList.push_back(integBodies[i].isPPN);
        this->forceParams.isJ2List.push_back(integBodies[i].isJ2);
        this->forceParams.J2List.push_back(integBodies[i].J2);
        this->forceParams.poleRAList.push_back(integBodies[i].poleRA);
        this->forceParams.poleDecList.push_back(integBodies[i].poleDec);
        this->forceParams.isNongravList.push_back(integBodies[i].isNongrav);
        this->forceParams.ngParamsList.push_back(integBodies[i].ngParams);
        this->forceParams.isMajorList.push_back(integBodies[i].isMajor);
        this->forceParams.isThrustingList.push_back(integBodies[i].isThrusting);
    }
    for (size_t i = 0; i < this->integParams.nSpice; i++) {
        this->forceParams.masses.push_back(spiceBodies[i].mass);
        this->forceParams.radii.push_back(spiceBodies[i].radius);
        this->forceParams.spiceIdList.push_back(spiceBodies[i].spiceId);
        this->forceParams.isPPNList.push_back(spiceBodies[i].isPPN);
        this->forceParams.isJ2List.push_back(spiceBodies[i].isJ2);
        this->forceParams.J2List.push_back(spiceBodies[i].J2);
        this->forceParams.poleRAList.push_back(spiceBodies[i].poleRA);
        this->forceParams.poleDecList.push_back(spiceBodies[i].poleDec);
        this->forceParams.isMajorList.push_back(spiceBodies[i].isMajor);
    }
    this->tStep.push_back(t);
    this->xIntegStep.push_back(xInteg);
}

void propSimulation::extend(real tf, std::vector<real> tEvalNew,
                            std::vector<std::vector<real>> xObserverNew) {
    // std::reverse(this->xObserver.begin(), this->xObserver.end());
    // std::reverse(this->observerInfo.begin(), this->observerInfo.end());
    // std::reverse(this->tEval.begin(), this->tEval.end());
    // std::reverse(this->radarObserver.begin(), this->radarObserver.end());
    // std::reverse(this->lightTimeEval.begin(), this->lightTimeEval.end());
    // std::reverse(this->xIntegEval.begin(), this->xIntegEval.end());
    // std::reverse(this->radarObsEval.begin(), this->radarObsEval.end());

    std::cout << "WARNING: The extend() method is under development and may "
                 "not work properly with the interpolation routines."
              << std::endl;

    // empty existing vectors from previous integration
    this->tEval.clear();
    this->xIntegEval.clear();
    this->observerInfo.clear();
    this->radarObserver.clear();
    this->xObserver.clear();
    this->lightTimeEval.clear();
    this->radarObsEval.clear();

    this->integParams.t0 = this->t;
    this->set_integration_parameters(tf, tEvalNew, this->tEvalUTC,
                                     this->evalApparentState,
                                     this->convergedLightTime, xObserverNew);
    this->integrate();
}
