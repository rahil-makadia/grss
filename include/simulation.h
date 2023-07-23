#ifndef SIMULATION_H
#define SIMULATION_H

#include "force.h"

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
    std::vector<real> pos;
    std::vector<real> vel;
    bool isPPN = false;
    bool isJ2 = false;
    bool isNongrav = false;
    bool isMajor = false;
    void set_J2(real J2, real poleRA, real poleDec);
};

class SpiceBody : public Body {
   private:
   public:
    int spiceId;
    bool isSpice = true;
    // constructor
    SpiceBody(std::string DEkernelPath, std::string name, int spiceID, real t0,
              real mass, real radius, Constants consts);
};

class IntegBody : public Body {
   private:
   public:
    bool isInteg = true;
    bool isThrusting = false;
    std::vector<std::vector<real>> covariance;
    NongravParamaters ngParams;
    // constructors
    IntegBody(std::string DEkernelPath, std::string name, real t0, real mass,
              real radius, std::vector<real> cometaryState,
              std::vector<std::vector<real>> covariance,
              NongravParamaters ngParams, Constants consts);
    IntegBody(std::string name, real t0, real mass, real radius,
              std::vector<real> pos, std::vector<real> vel,
              std::vector<std::vector<real>> covariance,
              NongravParamaters ngParams, Constants consts);
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

   public:
    // name and path to DE kernels
    std::string name;
    std::string DEkernelPath;
    // constructor and copy constructor
    propSimulation(std::string name, real t0, const int defaultSpiceBodies,
                   std::string DEkernelPath);
    propSimulation(std::string name, const propSimulation &simRef);

    // constants
    Constants consts;

    // integration parameters
    IntegrationParameters integParams;
    std::vector<real> tStep;
    std::vector<std::vector<real>> xIntegStep;

    // bodies and events
    std::vector<SpiceBody> spiceBodies;
    std::vector<IntegBody> integBodies;
    std::vector<ImpulseEvent> events;

    // preprocessor variables
    real t;
    std::vector<real> xInteg;
    ForceParameters forceParams;

    // interpolator variables
    size_t interpIdx = 0;
    bool tEvalUTC = false;
    bool evalApparentState = false;
    bool convergedLightTime = false;
    std::vector<std::vector<real>> xObserver;
    std::vector<std::vector<real>> observerInfo;
    real tEvalMargin = 0.0L;
    std::vector<real> tEval;
    std::vector<int> radarObserver;
    std::vector<std::vector<real>> lightTimeEval;
    std::vector<std::vector<real>> xIntegEval;
    std::vector<std::vector<real>> radarObsEval;

    // add/remove bodies and add events
    void add_spice_body(std::string DEkernelPath, std::string name, int spiceID,
                        real t0, real mass, real radius, Constants consts);
    void add_spice_body(SpiceBody body);
    void add_integ_body(std::string DEkernelPath, std::string name, real t0,
                        real mass, real radius, std::vector<real> cometaryState,
                        std::vector<std::vector<real>> covariance,
                        NongravParamaters ngParams, Constants consts);
    void add_integ_body(std::string name, real t0, real mass, real radius,
                        std::vector<real> pos, std::vector<real> vel,
                        std::vector<std::vector<real>> covariance,
                        NongravParamaters ngParams, Constants consts);
    void add_integ_body(IntegBody body);
    void remove_body(std::string name);
    void add_event(IntegBody body, real tEvent, std::vector<real> deltaV,
                   real multiplier = 1.0L);

    // setters
    void set_sim_constants(
        real du2m = 149597870700.0L, real tu2sec = 86400.0L,
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
        bool adaptiveTimestep = true, real dt0 = 0.0L, real dtMax = 6.0L,
        real dtMin = 5.0e-3L, real dtChangeFactor = 0.25L,
        real tolInteg = 1.0e-6L, real tolPC = 1.0e-16L);

    // getters
    std::vector<real> get_sim_constants();
    std::vector<real> get_integration_parameters();

    // preprocessor
    void preprocess();

    // integrator
    void integrate();

    // integration extension
    void extend(real tf, std::vector<real> tEvalNew = std::vector<real>(),
                std::vector<std::vector<real>> xObserverNew =
                    std::vector<std::vector<real>>());
};

#endif
