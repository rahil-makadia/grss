#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "force.h"
#include "gr15.h"
#include "grss.h"
#include "interpolate.h"
#include "simulation.h"
#include "utilities.h"

namespace py = pybind11;

PYBIND11_MODULE(prop_simulation, m) {
    m.doc() = "Simulation classes for the C++ GRSS orbit propagation code";

    // from utilities.h
    py::class_<Constants>(m, "Constants", R"mydelimiter(
        The Constants class contains physical constants and conversion factors
        used in the GRSS orbit propagation code.
        )mydelimiter")
        .def(py::init<>())
        .def_readwrite("du2m", &Constants::du2m, R"mydelimiter(
            Conversion factor from distance units to meters.
            )mydelimiter")
        .def_readwrite("tu2sec", &Constants::tu2sec, R"mydelimiter(
            Conversion factor from time units to seconds.
            )mydelimiter")
        .def_readwrite("G", &Constants::G, R"mydelimiter(
            Gravitational constant.
            )mydelimiter")
        .def_readwrite("clight", &Constants::clight, R"mydelimiter(
            Speed of light in a vacuum.
            )mydelimiter")
        .def_readwrite("j2000Jd", &Constants::j2000Jd, R"mydelimiter(
            Julian date of J2000 epoch.
            )mydelimiter")
        .def_readwrite("JdMinusMjd", &Constants::JdMinusMjd, R"mydelimiter(
            Difference between Julian date and modified Julian date.
            )mydelimiter");

    py::class_<IntegrationParameters>(m, "IntegrationParameters", R"mydelimiter(
        The IntegrationParameters class contains parameters used in the numerical
        integration of the equations of motion.
        )mydelimiter")
        .def(py::init<>())
        .def_readwrite("nInteg", &IntegrationParameters::nInteg, R"mydelimiter(
            Number of integrated bodies.
            )mydelimiter")
        .def_readwrite("nSpice", &IntegrationParameters::nSpice, R"mydelimiter(
            Number of bodies with SPICE ephemerides.
            )mydelimiter")
        .def_readwrite("nTotal", &IntegrationParameters::nTotal, R"mydelimiter(
            Total number of bodies. nTotal = nInteg + nSpice.
            )mydelimiter")
        .def_readwrite("t0", &IntegrationParameters::t0, R"mydelimiter(
            Initial time of integration (MJD TDB).
            )mydelimiter")
        .def_readwrite("tf", &IntegrationParameters::tf, R"mydelimiter(
            Final time of integration (MJD TDB).
            )mydelimiter")
        .def_readwrite("dt0", &IntegrationParameters::dt0, R"mydelimiter(
            Initial time step.
            )mydelimiter")
        .def_readwrite("dtMax", &IntegrationParameters::dtMax, R"mydelimiter(
            Maximum time step.
            )mydelimiter")
        .def_readwrite("dtMin", &IntegrationParameters::dtMin, R"mydelimiter(
            Minimum time step.
            )mydelimiter")
        .def_readwrite("dtChangeFactor", &IntegrationParameters::dtChangeFactor,
                       R"mydelimiter(
            Factor by which to limit the change in time step.
            )mydelimiter")
        .def_readwrite("adaptiveTimestep",
                       &IntegrationParameters::adaptiveTimestep, R"mydelimiter(
                       Flag to use adaptive time step.
                       )mydelimiter")
        .def_readonly("timestepCounter",
                      &IntegrationParameters::timestepCounter, R"mydelimiter(
                      Counter for number of time steps.
                      )mydelimiter")
        .def_readwrite("tolPC", &IntegrationParameters::tolPC, R"mydelimiter(
                        Tolerance for predictor-corrector within IAS15.
                        )mydelimiter")
        .def_readwrite("tolInteg", &IntegrationParameters::tolInteg,
                       R"mydelimiter(
                        Tolerance for integration.
                        )mydelimiter");

    py::class_<NongravParamaters>(m, "NongravParamaters", R"mydelimiter(
        The NongravParamaters class contains constants used for calculating the
        non-gravitational accelerations on integrated bodies.
        )mydelimiter")
        .def(py::init<>())
        .def_readwrite("a1", &NongravParamaters::a1, R"mydelimiter(
            Radial non-gravitational parameter a1.
            )mydelimiter")
        .def_readwrite("a2", &NongravParamaters::a2, R"mydelimiter(
            Transverse non-gravitational parameter a2.
            )mydelimiter")
        .def_readwrite("a3", &NongravParamaters::a3, R"mydelimiter(
            Normal non-gravitational parameter a3.
            )mydelimiter")
        .def_readwrite("alpha", &NongravParamaters::alpha, R"mydelimiter(
            Non-gravitational parameter alpha from Marsden et al. (1973).
            )mydelimiter")
        .def_readwrite("k", &NongravParamaters::k, R"mydelimiter(
            Non-gravitational parameter k from Marsden et al. (1973).
            )mydelimiter")
        .def_readwrite("m", &NongravParamaters::m, R"mydelimiter(
            Non-gravitational parameter m from Marsden et al. (1973).
            )mydelimiter")
        .def_readwrite("n", &NongravParamaters::n, R"mydelimiter(
            Non-gravitational parameter n from Marsden et al. (1973).
            )mydelimiter")
        .def_readwrite("r0_au", &NongravParamaters::r0_au, R"mydelimiter(
            Non-gravitational parameter r0 in AU from Marsden et al. (1973).
            )mydelimiter");

    m.def(
        "cometary_to_cartesian",
        [](real epochMjd, std::vector<real> cometaryState, real GM) {
            std::vector<real> cartesianState(6);
            cometary_to_cartesian(epochMjd, cometaryState, cartesianState, GM);
            return cartesianState;
        },
        py::arg("epochMjd"), py::arg("cometaryState"),
        py::arg("GM") = 2.9591220828559115e-4L, R"mydelimiter(
        Convert cometary state to cartesian state.

        Parameters
        ----------
        epochMjd : real
            Epoch in modified Julian date.
        cometaryState : list of real
            Cometary state vector.
        GM : real, optional
            Gravitational parameter of the central body, by default 0.00029591220828559115L.
        
        Returns
        -------
        cartesianState : list of real
            Cartesian state vector.
        )mydelimiter");
    m.def(
        "cartesian_to_cometary",
        [](real epochMjd, std::vector<real> cartesianState, real GM) {
            std::vector<real> cometaryState(6);
            cartesian_to_cometary(epochMjd, cartesianState, cometaryState, GM);
            return cometaryState;
        },
        py::arg("epochMjd"), py::arg("cartesianState"),
        py::arg("GM") = 2.9591220828559115e-4L, R"mydelimiter(
        Convert cartesian state to cometary state.

        Parameters
        ----------
        epochMjd : real
            Epoch in modified Julian date.
        cartesianState : list of real
            Cartesian state vector.
        GM : real, optional
            Gravitational parameter of the central body, by default 0.00029591220828559115L.

        Returns
        -------
        cometaryState : list of real
            Cometary state vector.
        )mydelimiter");

    // // from force.h
    py::class_<ForceParameters>(m, "ForceParameters", R"mydelimiter(
        The ForceParameters class contains constants used for calculating the
        forces on integrated bodies.
        )mydelimiter")
        .def(py::init<>())
        .def_readwrite("masses", &ForceParameters::masses, R"mydelimiter(
            Masses of the bodies.
            )mydelimiter")
        .def_readwrite("radii", &ForceParameters::radii, R"mydelimiter(
            Radii of the bodies.
            )mydelimiter")
        .def_readwrite("spiceIdList", &ForceParameters::spiceIdList,
                       R"mydelimiter(
            SPICE IDs of the bodies.
            )mydelimiter")
        .def_readwrite("ngParamsList", &ForceParameters::ngParamsList,
                       R"mydelimiter(
            Non-gravitational parameters of the bodies.
            )mydelimiter")
        .def_readwrite("isPPNList", &ForceParameters::isPPNList, R"mydelimiter(
            Whether the bodies are PPN bodies.
            )mydelimiter")
        .def_readwrite("isJ2List", &ForceParameters::isJ2List, R"mydelimiter(
            Whether the bodies are J2 bodies.
            )mydelimiter")
        .def_readwrite("J2List", &ForceParameters::J2List, R"mydelimiter(
            J2 parameters of the bodies.
            )mydelimiter")
        .def_readwrite("poleRAList", &ForceParameters::poleRAList,
                       R"mydelimiter(
            Right ascension of the poles of the bodies.
            )mydelimiter")
        .def_readwrite("poleDecList", &ForceParameters::poleDecList,
                       R"mydelimiter(
            Declination of the poles of the bodies.
            )mydelimiter")
        .def_readwrite("isNongravList", &ForceParameters::isNongravList,
                       R"mydelimiter(
            Whether the bodies have non-gravitational accelerations.
            )mydelimiter")
        .def_readwrite("isMajorList", &ForceParameters::isMajorList,
                       R"mydelimiter(
            Whether the bodies are major bodies (used for EIH PPN).
            )mydelimiter")
        .def_readwrite("isThrustingList", &ForceParameters::isThrustingList,
                       R"mydelimiter(
            Whether the bodies are thrusting.
            )mydelimiter");

    // from simulation.h
    py::class_<Body>(m, "Body", R"mydelimiter(
        The Body class contains the properties of an integrated or SPICE body.
        )mydelimiter")
        .def(py::init<>())
        .def_readwrite("t0", &Body::t0, R"mydelimiter(
            Initial MJD TDB time of the body. Same as the initial time of the propagator.
            )mydelimiter")
        .def_readwrite("mass", &Body::mass, R"mydelimiter(
            Mass of the body.
            )mydelimiter")
        .def_readwrite("radius", &Body::radius, R"mydelimiter(
            Radius of the body.
            )mydelimiter")
        .def_readwrite("J2", &Body::J2, R"mydelimiter(
            J2 parameter of the body.
            )mydelimiter")
        .def_readwrite("poleRA", &Body::poleRA, R"mydelimiter(
            Right ascension of the pole of the body.
            )mydelimiter")
        .def_readwrite("poleDec", &Body::poleDec, R"mydelimiter(
            Declination of the pole of the body.
            )mydelimiter")
        .def_readwrite("name", &Body::name, R"mydelimiter(
            Name of the body.
            )mydelimiter")
        .def_readwrite("pos", &Body::pos, R"mydelimiter(
            Position of the body at the initial time.
            )mydelimiter")
        .def_readwrite("vel", &Body::vel, R"mydelimiter(
            Velocity of the body at the initial time.
            )mydelimiter")
        .def_readwrite("isPPN", &Body::isPPN, R"mydelimiter(
            Whether the body is a PPN body.
            )mydelimiter")
        .def_readwrite("isMajor", &Body::isMajor, R"mydelimiter(
            Whether the body is a major body (used for EIH PPN).
            )mydelimiter")
        .def_readwrite("isJ2", &Body::isJ2, R"mydelimiter(
            Whether the body is a J2 body.
            )mydelimiter")
        .def_readwrite("isNongrav", &Body::isNongrav, R"mydelimiter(
            Whether the body has non-gravitational accelerations.
            )mydelimiter")
        .def("set_J2", &Body::set_J2, py::arg("J2"),
             py::arg("poleRA"), py::arg("poleDec"), R"mydelimiter(
            Set the J2 parameter of the body.

            Parameters
            ----------
            J2 : real
                J2 parameter of the body.
            poleRA : real
                Right ascension of the pole of the body.
            poleDec : real
                Declination of the pole of the body.

            Returns
            -------
            None : NoneType
                None.
            )mydelimiter");

    py::class_<SpiceBody, Body>(m, "SpiceBody", R"mydelimiter(
        The SpiceBody class contains the properties of a SPICE body.
        )mydelimiter")
        .def(py::init<std::string, std::string, int, real, real, real,
                      Constants>(),
             py::arg("DEkernelPath"), py::arg("name"), py::arg("spiceId"),
             py::arg("t0"), py::arg("mass"), py::arg("radius"),
             py::arg("constants"), R"mydelimiter(
            Constructor for the SpiceBody class.

            DEkernelPath : str
                Path to the SPICE DE kernel.
            name : str
                Name of the body.
            spiceId : int
                SPICE ID of the body.
            t0 : real
                Initial MJD TDB time of the body. Same as the initial time of the propagator.
            mass : real
                Mass of the body.
            radius : real
                Radius of the body.
            constants : propSimulation.Constants
                Constants of the simulation.
            )mydelimiter")
        .def_readwrite("spiceId", &SpiceBody::spiceId, R"mydelimiter(
            SPICE ID of the body.
            )mydelimiter")
        .def_readwrite("isSpice", &SpiceBody::isSpice, R"mydelimiter(
            Whether the body is a SPICE body. Always True.
            )mydelimiter");

    py::class_<IntegBody, Body>(m, "IntegBody", R"mydelimiter(
        The IntegBody class contains the properties of an integrated body.
        )mydelimiter")
        .def(py::init<std::string, std::string, real, real, real,
                      std::vector<real>, std::vector<std::vector<real>>,
                      NongravParamaters, Constants>(),
             py::arg("DEkernelPath"), py::arg("name"), py::arg("t0"),
             py::arg("mass"), py::arg("radius"), py::arg("cometaryState"),
             py::arg("covariance"), py::arg("ngParams"), py::arg("constants"),
             R"mydelimiter(
            Constructor for the IntegBody class.

            DEkernelPath : str
                Path to the SPICE DE kernel.
            name : str
                Name of the body.
            t0 : real
                Initial MJD TDB time of the body. Same as the initial time of the propagator.
            mass : real
                Mass of the body.
            radius : real
                Radius of the body.
            cometaryState : list of real
                Initial Heliocentric Ecliptic Cometary state of the body.
            covariance : list of list of real
                Covariance of the body's initial state.
            ngParams : propSimulation.NongravParamaters
                Non-gravitational parameters of the body.
            constants : propSimulation.Constants
                Constants of the simulation.
            )mydelimiter")

        .def(py::init<std::string, real, real, real, std::vector<real>,
                      std::vector<real>, std::vector<std::vector<real>>,
                      NongravParamaters, Constants>(),
             py::arg("name"), py::arg("t0"), py::arg("mass"), py::arg("radius"),
             py::arg("pos"), py::arg("vel"), py::arg("covariance"),
             py::arg("ngParams"), py::arg("constants"), R"mydelimiter(
            Constructor for the IntegBody class.

            name : str
                Name of the body.
            t0 : real
                Initial MJD TDB time of the body. Same as the initial time of the propagator.
            mass : real
                Mass of the body.
            radius : real
                Radius of the body.
            pos : list of real
                Initial barycentric Cartesian position of the body.
            vel : list of real
                Initial barycentric Cartesian velocity of the body.
            covariance : list of list of real
                Covariance of the body's initial state.
            ngParams : propSimulation.NongravParamaters
                Non-gravitational parameters of the body.
            constants : propSimulation.Constants
                Constants of the simulation.
            )mydelimiter")
        .def_readwrite("isInteg", &IntegBody::isInteg, R"mydelimiter(
            Whether the body is an integrated body. Always True.
            )mydelimiter")
        .def_readwrite("isThrusting", &IntegBody::isThrusting, R"mydelimiter(
            Whether the body is thrusting.
            )mydelimiter")
        .def_readwrite("covariance", &IntegBody::covariance, R"mydelimiter(
            Covariance of the body's initial state.
            )mydelimiter")
        .def_readwrite("ngParams", &IntegBody::ngParams, R"mydelimiter(
            Non-gravitational parameters of the body.
            )mydelimiter");

    py::class_<Event>(m, "Event", R"mydelimiter(
        The Event class contains the properties of an integration event.
        )mydelimiter")
        .def(py::init<>())
        .def_readwrite("t", &ImpulseEvent::t, R"mydelimiter(
            MJD TDB Time of the event.
            )mydelimiter")
        .def_readwrite("bodyName", &ImpulseEvent::bodyName, R"mydelimiter(
            Name of the integration body to apply the event to.
            )mydelimiter")
        .def_readwrite("bodyIndex", &ImpulseEvent::bodyIndex, R"mydelimiter(
            Index of the integration body to apply the event to.
            )mydelimiter");

    py::class_<ImpulseEvent, Event>(m, "ImpulseEvent", R"mydelimiter(
        The ImpulseEvent class contains the properties of an impulsive delta-V event.
        )mydelimiter")
        .def(py::init<>())
        .def_readwrite("deltaV", &ImpulseEvent::deltaV, R"mydelimiter(
            Delta-V of the event.
            )mydelimiter")
        .def_readwrite("multiplier", &ImpulseEvent::multiplier, R"mydelimiter(
            Multiplier on the delta-V the event.
            )mydelimiter");

    py::class_<propSimulation>(m, "propSimulation", R"mydelimiter(
        Class to perform an orbit propagation simulation.
        )mydelimiter")
        .def(py::init<std::string, real, const int, std::string>(),
             py::arg("name"), py::arg("t0"), py::arg("defaultSpiceBodies"),
             py::arg("DEkernelPath"), R"mydelimiter(
            Constructor for the propSimulation class.

            name : str
                Name of the simulation.
            t0 : real
                Initial MJD TDB time of the simulation.
            defaultSpiceBodies : int
                Version of the DE kernel to get the default SPICE bodies from.
            DEkernelPath : str
                Path to the SPICE DE kernel.
            )mydelimiter")
        .def(py::init<std::string, const propSimulation &>(), py::arg("name"),
             py::arg("simRef"), R"mydelimiter(
            Constructor for the propSimulation class.

            name : str
                Name of the simulation.
            simRef : propSimulation
                Simulation to copy.
            )mydelimiter")
        .def_readwrite("name", &propSimulation::name, R"mydelimiter(
            Name of the simulation.
            )mydelimiter")
        .def_readwrite("DEkernelPath", &propSimulation::DEkernelPath,
                       R"mydelimiter(
            Path to the SPICE DE kernel.
            )mydelimiter")
        .def_readwrite("consts", &propSimulation::consts, R"mydelimiter(
            Constants of the simulation. propSimulation.Constants object.
            )mydelimiter")
        .def_readwrite("integParams", &propSimulation::integParams,
                       R"mydelimiter(
            Integration parameters of the simulation. propSimulation.IntegParams object.
            )mydelimiter")
        .def_readwrite("tStep", &propSimulation::tStep, R"mydelimiter(
            Epochs of each step taken by the propagator.
            )mydelimiter")
        .def_readwrite("xIntegStep", &propSimulation::xIntegStep, R"mydelimiter(
            States of each step taken by the propagator.
            )mydelimiter")
        .def_readwrite("spiceBodies", &propSimulation::spiceBodies,
                       R"mydelimiter(
            SPICE bodies of the simulation. List of propSimulation.SpiceBodies objects.
            )mydelimiter")
        .def_readwrite("integBodies", &propSimulation::integBodies,
                       R"mydelimiter(
            Integration bodies of the simulation. List of propSimulation.IntegBody objects.
            )mydelimiter")
        .def_readwrite("events", &propSimulation::events, R"mydelimiter(
            Events of the simulation. List of propSimulation.Event objects.
            )mydelimiter")
        .def_readwrite("t", &propSimulation::t, R"mydelimiter(
            Current time of the simulation.
            )mydelimiter")
        .def_readwrite("xInteg", &propSimulation::xInteg, R"mydelimiter(
            Current states of each integration body in the simulation.
            )mydelimiter")
        .def_readwrite("forceParams", &propSimulation::forceParams,
                       R"mydelimiter(
            Force parameters of the simulation. propSimulation.ForceParams object.
            )mydelimiter")
        .def_readwrite("evalApparentState", &propSimulation::evalApparentState,
                       R"mydelimiter(
            Whether to evaluate the apparent state of the integration bodies.
            )mydelimiter")
        .def_readwrite("convergedLightTime",
                       &propSimulation::convergedLightTime, R"mydelimiter(
            Whether to use converged Newtonian light time correction.
            )mydelimiter")
        .def_readwrite("tEvalUTC", &propSimulation::tEvalUTC, R"mydelimiter(
            Whether the MJD evaluation time is in UTC for each value in propSimulation.tEval,
            as opposed to TDB.
            )mydelimiter")
        .def_readwrite("xObserver", &propSimulation::xObserver, R"mydelimiter(
            State of the observer for each value in propSimulation.tEval.
            )mydelimiter")
        .def_readwrite("observerInfo", &propSimulation::observerInfo,
                       R"mydelimiter(
            Observer information for each value in propSimulation.tEval.
            )mydelimiter")
        .def_readwrite("tEvalMargin", &propSimulation::tEvalMargin,
                       R"mydelimiter(
            Margin for allowing evaluation past the propagation start and end times.
            )mydelimiter")
        .def_readwrite("tEval", &propSimulation::tEval, R"mydelimiter(
            MJD Times to evaluate the states of the integrated bodies at.
            Can be TDB or UTC based on propSimulation.tEvalUTC.
            )mydelimiter")
        .def_readwrite("radarObserver", &propSimulation::radarObserver,
                       R"mydelimiter(
            Whether the observer for each value in propSimulation.tEval is for radar.
            )mydelimiter")
        .def_readwrite("lightTimeEval", &propSimulation::lightTimeEval,
                       R"mydelimiter(
            Light time from the observer to each integration body for each value in propSimulation.tEval.
            )mydelimiter")
        .def_readwrite("xIntegEval", &propSimulation::xIntegEval, R"mydelimiter(
            States of each integration body in the simulation for each value in propSimulation.tEval.
            )mydelimiter")
        .def_readwrite("radarObsEval", &propSimulation::radarObsEval,
                       R"mydelimiter(
            Radar observation of each integration body in the simulation for each value in propSimulation.tEval.
            )mydelimiter")
        .def("add_spice_body",
             static_cast<void (propSimulation::*)(std::string, std::string, int,
                                                  real, real, real, Constants)>(
                 &propSimulation::add_spice_body),
             py::arg("DEkernelPath"), py::arg("name"), py::arg("spiceId"),
             py::arg("t0"), py::arg("mass"), py::arg("radius"),
             py::arg("constants"), R"mydelimiter(
            Adds a SPICE body to the simulation.

            DEkernelPath : str
                Path to the SPICE DE kernel.
            name : str
                Name of the body.
            spiceId : int
                SPICE ID of the body.
            t0 : real
                Initial MJD epoch of the body. Must be in TDB. Same as the initial epoch of the simulation.
            mass : real
                Mass of the body.
            radius : real
                Radius of the body.
            constants : propSimulation.Constants
                Constants of the simulation.
            )mydelimiter")
        .def("add_spice_body",
             static_cast<void (propSimulation::*)(SpiceBody)>(
                 &propSimulation::add_spice_body),
             py::arg("body"), R"mydelimiter(
            Adds a SPICE body to the simulation.

            body : propSimulation.SpiceBody
                SPICE body to add to the simulation.
            )mydelimiter")
        .def("add_integ_body",
             static_cast<void (propSimulation::*)(
                 std::string, std::string, real, real, real, std::vector<real>,
                 std::vector<std::vector<real>>, NongravParamaters, Constants)>(
                 &propSimulation::add_integ_body),
             py::arg("DEkernelPath"), py::arg("name"), py::arg("t0"),
             py::arg("mass"), py::arg("radius"), py::arg("cometaryState"),
             py::arg("covariance"), py::arg("ngParams"), py::arg("constants"),
                R"mydelimiter(
            Adds an integration body to the simulation.

            DEkernelPath : str
                Path to the SPICE DE kernel.
            name : str
                Name of the body.
            t0 : real
                Initial MJD epoch of the body. Must be in TDB. Same as the initial epoch of the simulation.
            mass : real
                Mass of the body.
            radius : real
                Radius of the body.
            cometaryState : list of real
                Initial Heliocentric Ecliptic Cometary state of the body.
            covariance : list of list of real
                Covariance of the initial state.
            ngParams : propSimulation.NongravParamaters
                Nongravitational parameters of the body.
            constants : propSimulation.Constants
                Constants of the simulation.
            )mydelimiter")
        .def(
            "add_integ_body",
            static_cast<void (propSimulation::*)(
                std::string, real, real, real, std::vector<real>,
                std::vector<real>, std::vector<std::vector<real>>,
                NongravParamaters, Constants)>(&propSimulation::add_integ_body),
            py::arg("name"), py::arg("t0"), py::arg("mass"), py::arg("radius"),
            py::arg("pos"), py::arg("vel"), py::arg("covariance"),
            py::arg("ngParams"), py::arg("constants"), R"mydelimiter(
            Adds an integration body to the simulation.

            name : str
                Name of the body.
            t0 : real
                Initial MJD epoch of the body. Must be in TDB. Same as the initial epoch of the simulation.
            mass : real
                Mass of the body.
            radius : real
                Radius of the body.
            pos : list of real
                Initial barycentric Cartesian position of the body.
            vel : list of real
                Initial barycentric Cartesian velocity of the body.
            covariance : list of list of real
                Covariance of the initial state.
            ngParams : propSimulation.NongravParamaters
                Nongravitational parameters of the body.
            constants : propSimulation.Constants
                Constants of the simulation.
            )mydelimiter")
        .def("add_integ_body",
             static_cast<void (propSimulation::*)(IntegBody)>(
                 &propSimulation::add_integ_body),
             py::arg("body"), R"mydelimiter(
            Adds an integration body to the simulation.

            body : propSimulation.IntegBody
                Integration body to add to the simulation.
            )mydelimiter")
        .def("remove_body", &propSimulation::remove_body, py::arg("name"),
             R"mydelimiter(
            Removes a body from the simulation.

            Parameters
            ----------
            name : str
                Name of the body to remove.
            )mydelimiter")
        .def("add_event", &propSimulation::add_event, py::arg("body"),
             py::arg("tEvent"), py::arg("deltaV"), py::arg("multiplier") = 1.0L,
             R"mydelimiter(
            Adds an impulsive delta-V event to the simulation.

            Parameters
            ----------
            body : str
                Name of the body to apply the delta-V to.
            tEvent : real
                MJD Epoch of the event. Must be in TDB.
            deltaV : list of real
                Delta-V to apply to the body.
            multiplier : real
                Multiplier to apply to the delta-V.
            )mydelimiter")
        .def("set_sim_constants", &propSimulation::set_sim_constants,
             py::arg("du2m") = 149597870700.0L, py::arg("tu2sec") = 86400.0L,
             py::arg("G") = 6.6743e-11L /
                 (149597870700.0L * 149597870700.0L * 149597870700.0L) *
                 86400.0L * 86400.0L,
             py::arg("clight") = 299792458.0L / 149597870700.0L * 86400.0L,
             R"mydelimiter(
            Sets the constants of the simulation.

            Parameters
            ----------
            du2m : real
                Conversion factor from distance units to meters.
            tu2sec : real
                Conversion factor from time units to seconds.
            G : real
                Gravitational constant.
            clight : real
                Speed of light in a vacuum.
            )mydelimiter")
        .def("set_integration_parameters",
             &propSimulation::set_integration_parameters, py::arg("tf"),
             py::arg("tEval") = std::vector<real>(),
             py::arg("tEvalUTC") = false, py::arg("evalApparentState") = false,
             py::arg("convergedLightTims") = false,
             py::arg("observerInfo") = std::vector<std::vector<real>>(),
             py::arg("adaptiveTimestep") = true, py::arg("dt0") = 0.0L,
             py::arg("dtMax") = 6.0L, py::arg("dtMin") = 5.0e-3L,
             py::arg("dtChangeFactor") = 0.25L, py::arg("tolInteg") = 1.0e-6L,
             py::arg("tolPC") = 1.0e-16L, R"mydelimiter(
            Sets the integration parameters.

            Parameters
            ----------
            tf : real
                Final time of integration (MJD TDB).
            tEval : list of real
                MJD Times to evaluate the states of the integrated bodies at.
                Can be TDB or UTC based on tEvalUTC.
            tEvalUTC : bool
                Whether the MJD evaluation time is in UTC for each value in
                propSimulation.tEval, as opposed to TDB.
            evalApparentState : bool
                Whether to evaluate the apparent state of the integration bodies.
            convergedLightTimes : bool
                Whether to use converged Newtonian light time correction.
            observerInfo : list of list of real
                Observer information. Each list at least contains the central body SPICE ID
                (e.g., 399 for Earth) and the body-fixed longitude, latitude, and distance.
                This information might be repeated for bistatic radar observations.
            adaptiveTimestep : bool
                Flag to use adaptive time step for the propagation.
            dt0 : real
                Initial time step.
            dtMax : real
                Maximum time step.
            dtMin : real
                Minimum time step.
            dtChangeFactor : real
                Factor by which to limit the change in time step.
            tolInteg : real
                Tolerance for integration.
            tolPC : real
                Tolerance for predictor-corrector within IAS15.
            )mydelimiter")
        .def("get_sim_constants", &propSimulation::get_sim_constants,
                R"mydelimiter(
                Gets the constants of the simulation.

                Returns
                -------
                du2m : real
                    Conversion factor from distance units to meters.
                tu2sec : real
                    Conversion factor from time units to seconds.
                G : real
                    Gravitational constant.
                clight : real
                    Speed of light in a vacuum.
                j2000Jd : real
                    Julian date of J2000 epoch.
                JdMinusMjd : real
                    Difference between Julian date and modified Julian date.
                )mydelimiter")
        .def("get_integration_parameters",
             &propSimulation::get_integration_parameters, R"mydelimiter(
            Gets the integration parameters.

            Returns
            -------
            nInteg : int
                Number of integrated bodies.
            nSpice : int
                Number of bodies with SPICE ephemerides.
            nTotal : int
                Total number of bodies. nTotal = nInteg + nSpice.
            t0 : real
                Initial time of integration (MJD TDB).
            tf : real
                Final time of integration (MJD TDB).
            adaptiveTimestep : bool
                Flag to use adaptive time step for the propagation.
            dt0 : real
                Initial time step.
            dtMax : real
                Maximum time step.
            dtMin : real
                Minimum time step.
            dtChangeFactor : real
                Factor by which to limit the change in time step.
            tolInteg : real
                Tolerance for integration.
            tolPC : real
                Tolerance for predictor-corrector within IAS15.
            )mydelimiter")
        .def("preprocess", &propSimulation::preprocess, R"mydelimiter(
            Preprocesses the simulation and assembles the forceParams attribute.
            )mydelimiter")
        .def("integrate", &propSimulation::integrate, R"mydelimiter(
            Propagates the simulation using the Gauss-Radau integrator.
            )mydelimiter")
        .def("extend", &propSimulation::extend, py::arg("tf"),
             py::arg("tEvalNew") = std::vector<real>(),
             py::arg("xObserverNew") = std::vector<std::vector<real>>(), R"mydelimiter(
            Extends the simulation to a new final time.

            Parameters
            ----------
            tf : real
                New final time of integration (MJD TDB).
            tEvalNew : list of real
                Extra MJD Times to evaluate the states of the integrated bodies at.
                Can be TDB or UTC based on tEvalUTC.
            xObserverNew : list of list of real
                Extra observer information. Each list at least contains the central body SPICE ID
                (e.g., 399 for Earth) and the body-fixed longitude, latitude, and distance.
                This information might be repeated for bistatic radar observations.
            )mydelimiter");
}
