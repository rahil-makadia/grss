#ifndef SIMULATION_H
#define SIMULATION_H

#include "elements.h"

// forward declarations for main PropSimulation class
class PropSimulation;

void get_observer_state(const real &tObsMjd,
                        const std::vector<real> &observerInfo,
                        PropSimulation *propSim, const bool tObsInUTC,
                        std::vector<real> &observerState);

struct Constants {
    real du2m = 149597870700.0L;  // default au to m
    real tu2s = 86400.0;          // default day to sec
    real duptu2mps = du2m/tu2s;   // default au/day to m/s
    real G = 6.6743e-11L / (du2m * du2m * du2m) * tu2s *
        tu2s;                                  // default kg au^3 / day^2
    real clight = 299792458.0L / du2m * tu2s;  // default au/day
    real j2000Jd = 2451545.0;
    real JdMinusMjd = 2400000.5;
};

struct IntegrationParameters {
    size_t nInteg;
    size_t nSpice;
    size_t nTotal;
    size_t n2Derivs;
    real t0;
    real tf;
    real dt0;
    real dtMax;
    real dtMin;
    real dtChangeFactor;
    bool adaptiveTimestep;
    size_t timestepCounter;
    real tolPC;
    real tolInteg;
};

class Body {
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

struct NongravParameters {
    // from https://ssd.jpl.nasa.gov/horizons/manual.html
    real a1 = 0.0L;
    real a2 = 0.0L;
    real a3 = 0.0L;
    real alpha = 0.1112620426L;
    real k = 4.6142L;
    real m = 2.15L;
    real n = 5.093L;
    real r0_au = 2.808L;
};

class IntegBody : public Body {
   private:
   public:
    int spiceId = -99999;
    bool isCometary = false;
    std::vector<real> initState;
    bool isInteg = true;
    bool isThrusting = false;
    NongravParameters ngParams;
    size_t n2Derivs = 3;
    bool propStm = false;
    std::vector<real> stm;
    std::vector<std::vector<real>> dCartdState;
    // constructors
    IntegBody(std::string name, real t0, real mass, real radius,
              std::vector<real> cometaryState,
              NongravParameters ngParams);
    IntegBody(std::string name, real t0, real mass, real radius,
              std::vector<real> pos, std::vector<real> vel,
              NongravParameters ngParams);
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

struct BPlaneParameters {
    real x;
    real y;
    real z;
    std::vector<real> dx = std::vector<real>(6, 0.0L);
    std::vector<real> dy = std::vector<real>(6, 0.0L);
};

class CloseApproachParameters {
   private:
   public:
    real t;
    std::vector<real> xRel = std::vector<real>(6, 0.0L);
    real tMap;
    std::vector<real> xRelMap = std::vector<real>(6, 0.0L);
    real dist;
    real vel;
    real vInf;
    std::string flybyBody;
    int flybyBodyIdx;
    std::string centralBody;
    int centralBodyIdx;
    int centralBodySpiceId;
    bool impact;
    real tPeri;
    real tLin;
    std::vector<real> bVec = std::vector<real>(3, 0.0L);
    real bMag;
    real gravFocusFactor;
    BPlaneParameters kizner;
    BPlaneParameters opik;
    BPlaneParameters scaled;
    BPlaneParameters mtp;
    std::vector<real> dTLinMinusT = std::vector<real>(6, 0.0L);
    std::vector<real> dt = std::vector<real>(6, 0.0L);
    void get_ca_parameters(PropSimulation *propSim, const real &tMap);
    void print_summary(int prec=8);
};

class ImpactParameters : public CloseApproachParameters {
   private:
   public:
    std::vector<real> xRelBodyFixed = std::vector<real>(6, std::numeric_limits<real>::quiet_NaN());
    real lon;
    real lat;
    real alt;
    void get_impact_parameters(PropSimulation *propSim);
    void print_summary(int prec=8);
};

struct InterpolationParameters {
    std::vector<real> tStack;
    std::vector<std::vector<real>> xIntegStack;
    std::vector<std::vector<std::vector<real>>> bStack;
    std::vector<std::vector<real>> accIntegStack;
};

class PropSimulation {
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
    PropSimulation(std::string name, real t0, const int defaultSpiceBodies,
                   std::string DEkernelPath);
    PropSimulation(std::string name, const PropSimulation &simRef);

    // memory-mapped ephemeris
    Ephemeris ephem;
    std::vector<real> get_spiceBody_state(const real t, const std::string &bodyName);

    // constants
    Constants consts;

    // integration parameters
    IntegrationParameters integParams;
    bool parallelMode = false;

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
        real tolInteg = 1.0e-11L, real tolPC = 1.0e-16L);

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
