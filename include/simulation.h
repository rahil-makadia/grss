#ifndef SIMULATION_H
#define SIMULATION_H

#include "elements.h"

void get_observer_state(const real &tObsMjd,
                        const std::vector<real> &observerInfo,
                        propSimulation *propSim, const bool tObsInUTC,
                        std::vector<real> &observerState);

struct Body {
   private:
   public:
    real t0;
    real mass;
    real radius;
    real J2 = 0.0L;
    real poleRA = 0.0L;
    real poleDec = 90.0L;
    std::string name;
    int spiceId;
    real pos[3], vel[3], acc[3];
    bool isPPN = false;
    bool isJ2 = false;
    bool isNongrav = false;
    bool isMajor = false;
    real caTol = 0.1;
    void set_J2(real J2, real poleRA, real poleDec);
};

class SpiceBody : public Body {
   private:
   public:
    bool isSpice = true;
    // constructor
    SpiceBody(std::string name, int spiceID, real t0, real mass, real radius);
};

class IntegBody : public Body {
   private:
   public:
    int spiceId = -99999;
    bool isCometary = false;
    std::vector<real> initState;
    bool isInteg = true;
    bool isThrusting = false;
    NongravParamaters ngParams;
    size_t n2Derivs = 3;
    bool propStm = false;
    std::vector<real> stm;
    std::vector<std::vector<real>> dCartdState;
    // constructors
    IntegBody(std::string name, real t0, real mass, real radius,
              std::vector<real> cometaryState,
              NongravParamaters ngParams);
    IntegBody(std::string name, real t0, real mass, real radius,
              std::vector<real> pos, std::vector<real> vel,
              NongravParamaters ngParams);
    void prepare_stm();
};

class Event {
   private:
   public:
    real t;
    std::string bodyName;
    size_t bodyIndex;
};

class ImpulseEvent : public Event {
   private:
   public:
    std::vector<real> deltaV = {0.0L, 0.0L, 0.0L};
    real multiplier = 1.0L;
    void apply(const real &t, std::vector<real> &xInteg, const real &propDir);
};

class propSimulation {
   private:
    void prepare_for_evaluation(std::vector<real> &tEval,
                                std::vector<std::vector<real>> &observerInfo);
    // preprocessor
    bool isPreprocessed = false;
    void preprocess();

   public:
    // name and path to DE kernels
    std::string name;
    std::string DEkernelPath;
    // constructor and copy constructor
    propSimulation(std::string name, real t0, const int defaultSpiceBodies,
                   std::string DEkernelPath);
    propSimulation(std::string name, const propSimulation &simRef);

    // memory-mapped ephemeris
    Ephemeris ephem;
    std::vector<real> get_spiceBody_state(const real t, const std::string &bodyName);

    // constants
    Constants consts;

    // integration parameters
    IntegrationParameters integParams;

    // bodies and events
    std::vector<SpiceBody> spiceBodies;
    std::vector<IntegBody> integBodies;
    std::vector<ImpulseEvent> events;

    // close approach parameters
    std::vector<CloseApproachParameters> caParams;
    std::vector<ImpactParameters> impactParams;

    // preprocessor variables
    real t;
    std::vector<real> xInteg;

    // interpolation parameters
    InterpolationParameters interpParams;
    size_t interpIdx = 0;
    bool tEvalUTC = false;
    bool evalApparentState = false;
    bool evalMeasurements = false;
    bool convergedLightTime = false;
    std::vector<std::vector<real>> xObserver;
    std::vector<std::vector<real>> observerInfo;
    real tEvalMargin = 0.0L;
    std::vector<real> tEval;
    std::vector<int> radarObserver;
    std::vector<std::vector<real>> lightTimeEval;
    std::vector<std::vector<real>> xIntegEval;
    std::vector<std::vector<real>> opticalObs;
    std::vector<std::vector<real>> opticalPartials;
    std::vector<std::vector<real>> radarObs;
    std::vector<std::vector<real>> radarPartials;
    std::vector<real> interpolate(const real t);

    // add/remove bodies and add events
    void add_spice_body(SpiceBody body);
    void add_integ_body(IntegBody body);
    void remove_body(std::string name);
    void add_event(IntegBody body, real tEvent, std::vector<real> deltaV,
                   real multiplier = 1.0L);

    // setters
    void set_sim_constants(
        real du2m = 149597870700.0L, real tu2s = 86400.0L,
        real G = 6.6743e-11L /
            (149597870700.0L * 149597870700.0L * 149597870700.0L) * 86400.0L *
            86400.0L,
        real clight = 299792458.0L / 149597870700.0L * 86400.0L);
    void set_integration_parameters(
        real tf, std::vector<real> tEval = std::vector<real>(),
        bool tEvalUTC = false, bool evalApparentState = false,
        bool convergedLightTime = false,
        std::vector<std::vector<real>> observerInfo =
            std::vector<std::vector<real>>(),
        bool adaptiveTimestep = true, real dt0 = 0.0L, real dtMax = 21.0L,
        real dtMin = 5.0e-3L, real dtChangeFactor = 0.25L,
        real tolInteg = 1.0e-9L, real tolPC = 1.0e-16L);

    // getters
    std::vector<real> get_sim_constants();
    std::vector<real> get_integration_parameters();

    // integrator
    void integrate();

    // integration extension
    void extend(real tf, std::vector<real> tEvalNew = std::vector<real>(),
                std::vector<std::vector<real>> xObserverNew =
                    std::vector<std::vector<real>>());
};

#endif
