#ifndef SIMULATION_H
#define SIMULATION_H

#include "elements.h"

// forward declaration for main PropSimulation class
class PropSimulation;

// forward declaration for reconstruct_stm function (defined in stm.cpp)
/**
 * @brief Reconstruct the STM matrix from the flattened STM vector.
 */
std::vector<std::vector<real>> reconstruct_stm(const std::vector<real> &stm);

/**
 * @brief Get the name of body-fixed frame for a given SPICE ID.
 */
void get_baseBodyFrame(const int &spiceId, const real &tMjdTDB,
                       std::string &baseBodyFrame);

/**
 * @brief Get the observer state for a given time.
 */
void get_observer_state(const real &tObsMjd,
                        const std::vector<real> &observerInfo,
                        PropSimulation *propSim, const bool tObsInUTC,
                        std::vector<real> &observerState);

/**
 * @brief Structure to hold constants used in a PropSimulation.
 * 
 * @param du2m Distance unit conversion factor (default: AU to meters).
 * @param tu2s Time unit conversion factor (default: days to seconds).
 * @param duptu2mps Distance unit per time unit conversion factor (default: AU/day to m/s).
 * @param G Gravitational constant (default: kg AU^3/day^2).
 * @param clight Speed of light (default: AU/day).
 * @param j2000Jd Julian date of J2000 epoch.
 * @param JdMinusMjd Julian date minus modified Julian date.
 */
struct Constants {
    real du2m = 149597870700.0L;
    real tu2s = 86400.0;
    real duptu2mps = du2m / tu2s;
    real G = 6.6743e-11L / (du2m * du2m * du2m) * tu2s * tu2s;
    real clight = 299792458.0L / du2m * tu2s;
    real j2000Jd = 2451545.0;
    real JdMinusMjd = 2400000.5;
};

/**
 * @brief Structure to hold integration parameters used in a PropSimulation.
 * 
 * @param nInteg Number of integrated bodies.
 * @param nSpice Number of SPICE bodies (queried from PropSimulation Ephemeris object).
 * @param nTotal Total number of bodies (nInteg + nSpice).
 * @param n2Derivs Number of second derivatives (usually 3 for each body unless STM is propagated).
 * @param t0 Initial time.
 * @param tf Final time.
 * @param dt0 Initial timestep.
 * @param dtMax Maximum timestep.
 * @param dtMin Minimum timestep.
 * @param dtChangeFactor Maximum factor by which to change timestep.
 * @param adaptiveTimestep Flag to use adaptive timestep.
 * @param timestepCounter Timestep counter.
 * @param tolPC Tolerance for Gauss-Radau predictor-corrector loop.
 * @param tolInteg Tolerance for integrator.
 */
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

/**
 * @brief Parent class for all bodies in a PropSimulation.
 * 
 * @param t0 Initial time.
 * @param mass Mass of the body.
 * @param radius Radius of the body.
 * @param J2 J2 coefficient for oblateness.
 * @param poleRA Right ascension of the pole.
 * @param poleDec Declination of the pole.
 * @param name Name of the body.
 * @param spiceId SPICE ID of the body.
 * @param pos Position of the body.
 * @param vel Velocity of the body.
 * @param acc Acceleration of the body.
 * @param isPPN Flag to indicate if the body is a PPN body (for EIH relativity loop 1).
 * @param isJ2 Flag to indicate if the body has J2 oblateness.
 * @param isNongrav Flag to indicate if the body has nongravitational forces.
 * @param isMajor Flag to indicate if the body is a major body (for EIH relativity loop 2).
 * @param caTol Distance tolerance for close approach detection (default: 0.1 AU).
 */
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
    /**
     * @brief Set the J2 for the body.
     */
    void set_J2(real J2, real poleRA, real poleDec);
};

/**
 * @brief Class for SPICE bodies in a PropSimulation.
 * 
 * @param isSpice Flag to indicate if the body is a SPICE body.
 */
class SpiceBody : public Body {
   private:
   public:
    bool isSpice = true;
    /**
     * @brief Construct a new SpiceBody object
     */
    SpiceBody(std::string name, int spiceId, real t0, real mass, real radius);
};

/**
 * @brief Structure to hold parameters for asteroi/comet nongravitational forces.
 * 
 * @param a1 Parameter for radial nongravitational forces.
 * @param a2 Parameter for transverse nongravitational forces.
 * @param a3 Parameter for normal nongravitational forces.
 * @param a1Est Flag to indicate if a1 is estimated (for orbit determination).
 * @param a2Est Flag to indicate if a2 is estimated (for orbit determination).
 * @param a3Est Flag to indicate if a3 is estimated (for orbit determination).
 * @param alpha Normalizing factor for nongravitational forces.
 * @param k Exponent for nongravitational forces.
 * @param m Exponent for nongravitational forces.
 * @param n Exponent for nongravitational forces.
 * @param r0_au Normalizing distance for nongravitational forces.
 */
struct NongravParameters {
    // from https://ssd.jpl.nasa.gov/horizons/manual.html
    real a1 = 0.0L;
    real a2 = 0.0L;
    real a3 = 0.0L;
    bool a1Est = false;
    bool a2Est = false;
    bool a3Est = false;
    real alpha = 0.1112620426L;
    real k = 4.6142L;
    real m = 2.15L;
    real n = 5.093L;
    real r0_au = 2.808L;
};

/**
 * @brief Class for integrated bodies in a PropSimulation.
 * 
 * @param spiceId SPICE ID of the body.
 * @param isCometary Flag to indicate if the body is initialized with a cometary state (as opposed to Cartesian state).
 * @param initState Initial state of the body.
 * @param isInteg Flag to indicate if the body is integrated.
 * @param isThrusting Flag to indicate if the body is thrusting.
 * @param ngParams Nongravitational parameters for the body.
 * @param n2Derivs Number of second derivatives for the body.
 * @param propStm Flag to indicate if the state transition matrix is propagated.
 * @param stm State transition matrix.
 * @param dCartdState Derivatives of initial Cartesian state with respect to initial state.
 */
class IntegBody : public Body {
   private:
   public:
    int spiceId = -99999;
    bool isCometary = false;
    std::vector<real> initState;
    std::vector<real> initCart;
    bool isInteg = true;
    bool isThrusting = false;
    NongravParameters ngParams;
    size_t n2Derivs = 3;
    bool propStm = false;
    std::vector<real> stm;
    std::vector<std::vector<real>> dCartdState;
    /**
     * @brief Construct a new IntegBody object.
     */
    IntegBody(std::string name, real t0, real mass, real radius,
              std::vector<real> cometaryState,
              NongravParameters ngParams);
    /**
     * @brief Construct a new IntegBody object.
     */
    IntegBody(std::string name, real t0, real mass, real radius,
              std::vector<real> pos, std::vector<real> vel,
              NongravParameters ngParams);
    /**
     * @brief Prepare the state transition matrix for propagation.
     */
    void prepare_stm();
};

/**
 * @brief Class for events in a PropSimulation.
 * 
 * @param t Time of the event.
 * @param bodyName Name of the IntegBody for the event
 * @param bodyIndex Index of the IntegBody in the PropSimulation.
 */
class Event {
   private:
   public:
    real t;
    std::string bodyName;
    size_t bodyIndex;
};

/**
 * @brief Class for impulse events in a PropSimulation.
 * 
 * @param deltaV Delta-V for the impulse.
 * @param multiplier Multiplier for the Delta-V.
 */
class ImpulseEvent : public Event {
   private:
   public:
    std::vector<real> deltaV = {0.0L, 0.0L, 0.0L};
    real multiplier = 1.0L;
    /**
     * @brief Apply the impulse event to the body.
     * 
     * @param[in] t Time of the event.
     * @param[inout] xInteg State of the body.
     * @param[in] propDir Direction of propagation.
     */
    void apply(const real &t, std::vector<real> &xInteg, const real &propDir);
};

/** 
 * @brief Structure to hold parameters for a close approach B-plane.
 * 
 * @param x X-coordinate of the B-plane.
 * @param y Y-coordinate of the B-plane.
 * @param z Z-coordinate of the B-plane (usually 0).
 * @param dx Derivatives of x with respect to the state at the close approach.
 * @param dy Derivatives of y with respect to the state at the close approach.
*/
struct BPlaneParameters {
    real x;
    real y;
    real z;
    std::vector<real> dx = std::vector<real>(6, std::numeric_limits<real>::quiet_NaN());
    std::vector<real> dy = std::vector<real>(6, std::numeric_limits<real>::quiet_NaN());
};

/**
 * @brief Class for close approach parameters in a PropSimulation.
 * 
 * @param t Time of the close approach.
 * @param xRel Relative state at the close approach.
 * @param tMap Time of the B-plane map.
 * @param xRelMap Relative state at the B-plane map time.
 * @param dist Distance at the close approach.
 * @param vel Relative velocity at the close approach.
 * @param vInf Hyperbolic excess velocity of the close approach.
 * @param flybyBody Name of the flyby body.
 * @param flybyBodyIdx Index of the flyby body in the PropSimulation.
 * @param centralBody Name of the central body.
 * @param centralBodyIdx Index of the central body in the PropSimulation.
 * @param centralBodySpiceId SPICE ID of the central body.
 * @param impact Flag to indicate if the close approach is an impact.
 * @param tPeri Time of periapsis.
 * @param tLin Linearized time of intersection on the B-plane.
 * @param bVec B-vector of the B-plane.
 * @param bMag Magnitude of the B-vector.
 * @param gravFocusFactor Gravitational focus factor to account for nonlinear effects from central body.
 * @param kizner Kizner B-plane parameters.
 * @param opik Öpik B-plane parameters.
 * @param scaled Scaled Kizner B-plane parameters.
 * @param mtp Modified Target Plane parameters.
 * @param dTLinMinusT Partial derivatives of the (linearized intersection minus map) time with respect to the state at the close approach.
 * @param dt Partial derivatives of the time of closest approach with respect to the state at the close approach.
 */
class CloseApproachParameters {
   private:
   public:
    real t;
    std::vector<real> xRel;
    real tMap;
    std::vector<real> xRelMap;
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
    /**
     * @brief Get the close approach parameters.
     */
    void get_ca_parameters(PropSimulation *propSim, const real &tMap);
    /**
     * @brief Print a summary of the close approach parameters.
     */
    void print_summary(int prec=8);
};

/**
 * @brief Class for impact parameters in a PropSimulation.
 * 
 * @param xRelBodyFixed Relative state at the impact time in body-fixed frame.
 * @param lon Longitude of the impact [rad].
 * @param lat Latitude of the impact [rad].
 * @param alt Altitude of the impact [km].
 */
class ImpactParameters : public CloseApproachParameters {
   private:
   public:
    std::vector<real> xRelBodyFixed = std::vector<real>(6, std::numeric_limits<real>::quiet_NaN());
    real lon;
    real lat;
    real alt;
    /**
     * @brief Get the impact parameters.
     */
    void get_impact_parameters(PropSimulation *propSim);
    /**
     * @brief Print a summary of the impact parameters.
     */
    void print_summary(int prec=8);
};

/**
 * @brief Structure to hold interpolation parameters for a PropSimulation.
 * 
 * @param tStack Stack of integrator epochs.
 * @param xIntegStack Stack of states at integrator epochs.
 * @param bStack Stack of interpolation coefficients at integrator epochs.
 * @param accIntegStack Stack of accelerations at integrator epochs.
 */
struct InterpolationParameters {
    std::vector<real> tStack;
    std::vector<std::vector<real>> xIntegStack;
    std::vector<std::vector<std::vector<real>>> bStack;
    std::vector<std::vector<real>> accIntegStack;
};

/**
 * @brief Class for a PropSimulation.
 * 
 * @param name Name of the simulation.
 * @param DEkernelPath Path to DE kernels.
 * @param ephem Ephemeris object for the simulation.
 * @param consts Constants object for the simulation.
 * @param integParams Integration parameters for the simulation.
 * @param spiceBodies Vector of SpiceBody objects in the simulation.
 * @param integBodies Vector of IntegBody objects in the simulation.
 * @param events Vector of ImpulseEvent objects in the simulation.
 * @param caParams Vector of CloseApproachParameters objects in the simulation.
 * @param impactParams Vector of ImpactParameters objects in the simulation.
 * @param isPreprocessed Flag to indicate if the simulation is preprocessed.
 * @param t Current time of the simulation.
 * @param xInteg Current state of the simulation.
 * @param interpParams Interpolation parameters for the simulation.
 * @param interpIdx Index of the interpolation parameters.
 * @param tEvalUTC Flag to indicate if the evaluation times are in UTC.
 * @param evalApparentState Flag to indicate if the apparent state is evaluated.
 * @param evalMeasurements Flag to indicate if measurements are evaluated.
 * @param convergedLightTime Flag to indicate if converged light time needs to be computed.
 * @param xObserver Observer states array.
 * @param observerInfo Observer information array.
 * @param tEvalMargin Margin for interpolation outside integration arc.
 * @param tEval Vector of times at which to evaluate the integrated state.
 * @param radarObserver Vector of radar observer flags (0=none, 1=delay, 2=doppler).
 * @param lightTimeEval Vector of computed light times.
 * @param xIntegEval Vector of integrated states at evaluation times.
 * @param opticalObs Vector of optical observations.
 * @param opticalPartials Vector of optical observation partials.
 * @param radarObs Vector of radar observations.
 * @param radarPartials Vector of radar observation partials.
 */
class PropSimulation {
   private:
    /**
     * @brief Prepare the simulation for evaluation.
     */
    void prepare_for_evaluation(std::vector<real> &tEval,
                                std::vector<std::vector<real>> &observerInfo);
    // preprocessor
    bool isPreprocessed = false;
    /**
     * @brief Preprocess the simulation before integration.
     */
    void preprocess();
   public:
    std::string name;
    std::string DEkernelPath;
    SpkEphemeris spkEphem;
    PckEphemeris pckEphem;
    /**
     * @brief Construct a new PropSimulation object.
     */
    PropSimulation(std::string name, real t0, const int defaultSpiceBodies,
                   std::string DEkernelPath);
    /**
     * @brief Construct a new PropSimulation object from a reference simulation.
     */
    PropSimulation(std::string name, const PropSimulation &simRef);
    /**
     * @brief Memory map the ephemeris files.
     */
    void map_ephemeris();
    /**
     * @brief Unmap the ephemeris files.
     */
    void unmap_ephemeris();
    /**
     * @brief Get the state of a SpiceBody in the PropSimulation at a given time.
     */
    std::vector<real> get_spiceBody_state(const real t, const std::string &bodyName);
    Constants consts;
    IntegrationParameters integParams;
    std::vector<SpiceBody> spiceBodies;
    std::vector<IntegBody> integBodies;
    std::vector<ImpulseEvent> events;
    std::vector<CloseApproachParameters> caParams;
    std::vector<ImpactParameters> impactParams;
    real t;
    std::vector<real> xInteg;
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
    /**
     * @brief Get the integrated state from the simulation at a given time.
     */
    std::vector<real> interpolate(const real t);
    /**
     * @brief Add a SpiceBody to the simulation.
     */
    void add_spice_body(SpiceBody body);
    /**
     * @brief Add an IntegBody to the simulation.
     */
    void add_integ_body(IntegBody body);
    /**
     * @brief Remove a body from the simulation.
     */
    void remove_body(std::string name);
    /**
     * @brief Add an impulse event to the simulation.
     */
    void add_event(IntegBody body, real tEvent, std::vector<real> deltaV,
                   real multiplier = 1.0L);
    /**
     * @brief Set the values of the PropSimulation Constants object.
     */
    void set_sim_constants(
        real du2m = 149597870700.0L, real tu2s = 86400.0L,
        real G = 6.6743e-11L /
            (149597870700.0L * 149597870700.0L * 149597870700.0L) * 86400.0L *
            86400.0L,
        real clight = 299792458.0L / 149597870700.0L * 86400.0L);
    /**
     * @brief Set the values of the PropSimulation IntegrationParameters object.
     */
    void set_integration_parameters(
        real tf, std::vector<real> tEval = std::vector<real>(),
        bool tEvalUTC = false, bool evalApparentState = false,
        bool convergedLightTime = false,
        std::vector<std::vector<real>> observerInfo =
            std::vector<std::vector<real>>(),
        bool adaptiveTimestep = true, real dt0 = 0.0L, real dtMax = 21.0L,
        real dtMin = 5.0e-3L, real dtChangeFactor = 0.25L,
        real tolInteg = 1.0e-11L, real tolPC = 1.0e-16L);
    /**
     * @brief Get the constants used in the simulation.
     */
    std::vector<real> get_sim_constants();
    /**
     * @brief Get the integration parameters used in the simulation.
     */
    std::vector<real> get_integration_parameters();
    /**
     * @brief Integrate the simulation to the final time.
     */
    void integrate();
    /**
     * @brief Extend the simulation to a new final time.
     */
    void extend(real tf, std::vector<real> tEvalNew = std::vector<real>(),
                std::vector<std::vector<real>> xObserverNew =
                    std::vector<std::vector<real>>());
    /**
     * @brief Save the simulation to a file.
     */
    void save(std::string filename);
};

#endif
