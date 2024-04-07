#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "grss.h"

namespace py = pybind11;

PYBIND11_MODULE(libgrss, m) {
    m.doc() = "GRSS Python bindings for C++ library";

    py::class_<Constants>(m, "Constants", R"mydelimiter(
        The Constants class contains physical constants and conversion factors
        used in the GRSS orbit propagation code.
        )mydelimiter")
        .def(py::init<>())
        .def_readwrite("du2m", &Constants::du2m, R"mydelimiter(
            Conversion factor from distance units to meters.
            )mydelimiter")
        .def_readwrite("tu2s", &Constants::tu2s, R"mydelimiter(
            Conversion factor from time units to seconds.
            )mydelimiter")
        .def_readwrite("duptu2mps", &Constants::duptu2mps, R"mydelimiter(
            Conversion factor from distance units per time units to meters per second.
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
        .def_readwrite("n2Derivs", &IntegrationParameters::n2Derivs,
                       R"mydelimiter(
            Number of second derivatives.
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

    py::class_<NongravParameters>(m, "NongravParameters", R"mydelimiter(
        The NongravParameters class contains constants used for calculating the
        non-gravitational accelerations on integrated bodies.
        )mydelimiter")
        .def(py::init<>())
        .def_readwrite("a1", &NongravParameters::a1, R"mydelimiter(
            Radial non-gravitational parameter a1.
            )mydelimiter")
        .def_readwrite("a2", &NongravParameters::a2, R"mydelimiter(
            Transverse non-gravitational parameter a2.
            )mydelimiter")
        .def_readwrite("a3", &NongravParameters::a3, R"mydelimiter(
            Normal non-gravitational parameter a3.
            )mydelimiter")
        .def_readwrite("a1Est", &NongravParameters::a1Est, R"mydelimiter(
            Flag for whether to estimate a1 (used when propagating STM).
            )mydelimiter")
        .def_readwrite("a2Est", &NongravParameters::a2Est, R"mydelimiter(
            Flag for whether to estimate a2 (used when propagating STM).
            )mydelimiter")
        .def_readwrite("a3Est", &NongravParameters::a3Est, R"mydelimiter(
            Flag for whether to estimate a3 (used when propagating STM).
            )mydelimiter")
        .def_readwrite("alpha", &NongravParameters::alpha, R"mydelimiter(
            Non-gravitational parameter alpha from Marsden et al. (1973).
            )mydelimiter")
        .def_readwrite("k", &NongravParameters::k, R"mydelimiter(
            Non-gravitational parameter k from Marsden et al. (1973).
            )mydelimiter")
        .def_readwrite("m", &NongravParameters::m, R"mydelimiter(
            Non-gravitational parameter m from Marsden et al. (1973).
            )mydelimiter")
        .def_readwrite("n", &NongravParameters::n, R"mydelimiter(
            Non-gravitational parameter n from Marsden et al. (1973).
            )mydelimiter")
        .def_readwrite("r0_au", &NongravParameters::r0_au, R"mydelimiter(
            Non-gravitational parameter r0 in AU from Marsden et al. (1973).
            )mydelimiter");

    py::class_<InterpolationParameters>(m, "InterpolationParameters",
                                        R"mydelimiter(
        The InterpolationParameters class contains parameters used for
        interpolation of the states of the integrated bodies.
        )mydelimiter")
        .def(py::init<>())
        .def_readwrite("tStack", &InterpolationParameters::tStack,
                       R"mydelimiter(
            Stack of times used for interpolation at steps taken by the integrator.
            )mydelimiter")
        .def_readwrite("xIntegStack", &InterpolationParameters::xIntegStack,
                       R"mydelimiter(
            Stack of states of the integrated bodies used for interpolation at steps taken by the integrator.
            )mydelimiter")
        .def_readwrite("bStack", &InterpolationParameters::bStack,
                       R"mydelimiter(
            Stack of b matrices used for the interpolating coefficients at steps taken by the integrator.
            )mydelimiter")
        .def_readwrite("accIntegStack", &InterpolationParameters::accIntegStack,
                       R"mydelimiter(
            Stack of accelerations of the integrated bodies at steps taken by the integrator.
            )mydelimiter");

    py::class_<BPlaneParameters>(m, "BPlaneParameters", R"mydelimiter(
        The BPlaneParameters class contains parameters used for calculating the
        B-plane of a close approach.
        )mydelimiter")
        .def(py::init<>())
        .def_readwrite("x", &BPlaneParameters::x, R"mydelimiter(
            X coordinate of the B-plane.
            )mydelimiter")
        .def_readwrite("y", &BPlaneParameters::y, R"mydelimiter(
            Y coordinate of the B-plane.
            )mydelimiter")
        .def_readwrite("z", &BPlaneParameters::z, R"mydelimiter(
            Z coordinate of the B-plane.
            )mydelimiter")
        .def_readwrite("dx", &BPlaneParameters::dx, R"mydelimiter(
            Partial of X coordinate of the B-plane with respect to state.
            )mydelimiter")
        .def_readwrite("dy", &BPlaneParameters::dy, R"mydelimiter(
            Partial of Y coordinate of the B-plane with respect to state.
            )mydelimiter");

    py::class_<CloseApproachParameters>(m, "CloseApproachParameters",
                                        R"mydelimiter(
        The CloseApproachParameters class contains parameters used for calculating
        the close approach between two integrated or SPICE bodies.
        )mydelimiter")
        .def(py::init<>())
        .def_readwrite("t", &CloseApproachParameters::t, R"mydelimiter(
            Time of the close approach.
            )mydelimiter")
        .def_readwrite("xRel", &CloseApproachParameters::xRel,
                       R"mydelimiter(
            Relative state of the close approach.
            )mydelimiter")
        .def_readwrite("tMap", &CloseApproachParameters::tMap, R"mydelimiter(
            Time of the mapping point.
            )mydelimiter")
        .def_readwrite("xRelMap", &CloseApproachParameters::xRelMap,
                       R"mydelimiter(
            Relative state of the mapping point.
            )mydelimiter")
        .def_readwrite("dist", &CloseApproachParameters::dist, R"mydelimiter(
            Distance of the close approach.
            )mydelimiter")
        .def_readwrite("vel", &CloseApproachParameters::vel, R"mydelimiter(
            Velocity of the close approach.
            )mydelimiter")
        .def_readwrite("vInf", &CloseApproachParameters::vInf,
                       R"mydelimiter(
            Hyperbolic excess velocity of the close approach.
            )mydelimiter")
        .def_readwrite("flybyBody", &CloseApproachParameters::flybyBody,
                       R"mydelimiter(
            Name of the flyby body.
            )mydelimiter")
        .def_readwrite("flybyBodyIdx", &CloseApproachParameters::flybyBodyIdx,
                       R"mydelimiter(
            Index of the flyby body.
            )mydelimiter")
        .def_readwrite("centralBody", &CloseApproachParameters::centralBody,
                       R"mydelimiter(
            Name of the central body.
            )mydelimiter")
        .def_readwrite("centralBodyIdx", &CloseApproachParameters::centralBodyIdx,
                       R"mydelimiter(
            Index of the central body.
            )mydelimiter")
        .def_readwrite("centralBodySpiceId",
                       &CloseApproachParameters::centralBodySpiceId,
                       R"mydelimiter(
            SPICE ID of the central body.
            )mydelimiter")
        .def_readwrite("impact", &CloseApproachParameters::impact,
                       R"mydelimiter(
            Whether the close approach is an impact when accounting for gravitational focusing.
            )mydelimiter")
        .def_readwrite("tPeri", &CloseApproachParameters::tPeri, R"mydelimiter(
            Time of periapsis (according to Keplerian hyperbolic motion).
            )mydelimiter")
        .def_readwrite("tLin", &CloseApproachParameters::tLin, R"mydelimiter(
            Linearized time of periapsis
            )mydelimiter")
        .def_readwrite("bVec", &CloseApproachParameters::bVec, R"mydelimiter(
            B-vector of the close approach.
            )mydelimiter")
        .def_readwrite("bMag", &CloseApproachParameters::bMag, R"mydelimiter(
            Magnitude of the B-vector of the close approach (Impact parameter).
            )mydelimiter")
        .def_readwrite("gravFocusFactor",
                       &CloseApproachParameters::gravFocusFactor, R"mydelimiter(
            Lambda parameter of the close approach (gravitational focusing).
            )mydelimiter")
        .def_readwrite("kizner", &CloseApproachParameters::kizner,
                       R"mydelimiter(
            Kizner B-plane parameters of the close approach.
            )mydelimiter")
        .def_readwrite("opik", &CloseApproachParameters::opik, R"mydelimiter(
            Öpik B-plane parameters of the close approach.
            )mydelimiter")
        .def_readwrite("scaled", &CloseApproachParameters::scaled,
                       R"mydelimiter(
            Scaled B-plane parameters of the close approach.
            )mydelimiter")
        .def_readwrite("mtp", &CloseApproachParameters::mtp, R"mydelimiter(
            Modified Target Plane (MTP) B-plane parameters of the close approach.
            )mydelimiter")
        .def_readwrite("dTLinMinusT", &CloseApproachParameters::dTLinMinusT,
                       R"mydelimiter(
            Partials of difference between linearized time of periapsis and time of periapsis with respect to CA state.
            )mydelimiter")
        .def_readwrite("dt", &CloseApproachParameters::dt, R"mydelimiter(
            Partials of time of periapsis with respect to CA state.
            )mydelimiter")
        .def("get_ca_parameters", &CloseApproachParameters::get_ca_parameters,
             py::arg("propSim"), py::arg("tMap"), R"mydelimiter(
            Calculate the close approach parameters.

            Parameters
            ----------
            propSim : PropSimulation
                Simulation containing the close approach.
            tMap : real
                Time of the mapping point.

            Returns
            -------
            None : NoneType
                None.
            )mydelimiter")
        .def("print_summary", &CloseApproachParameters::print_summary,
             py::arg("prec") = 8, R"mydelimiter(
            Print a summary of the close approach parameters.

            Parameters
            ----------
            prec : int, optional
                Precision of the printed values, by default 8.

            Returns
            -------
            None : NoneType
                None.
            )mydelimiter");

    py::class_<ImpactParameters, CloseApproachParameters>(m, "ImpactParameters",
                                                          R"mydelimiter(
        The ImpactParameters class contains parameters used for calculating
        the impact between two integrated or SPICE bodies.
        )mydelimiter")
        .def(py::init<>())
        .def_readwrite("xRelBodyFixed", &ImpactParameters::xRelBodyFixed,
                       R"mydelimiter(
            Relative state of the impact in the body-fixed frame of the central body.
            )mydelimiter")
        .def_readwrite("lon", &ImpactParameters::lon, R"mydelimiter(
            Longitude of the impact.
            )mydelimiter")
        .def_readwrite("lat", &ImpactParameters::lat, R"mydelimiter(
            Latitude of the impact.
            )mydelimiter")
        .def_readwrite("alt", &ImpactParameters::alt, R"mydelimiter(
            Altitude of the impact.
            )mydelimiter")
        .def("print_summary", &ImpactParameters::print_summary,
             py::arg("prec") = 8, R"mydelimiter(
            Print a summary of the impact parameters.

            Parameters
            ----------
            prec : int, optional
                Precision of the printed values, by default 8.

            Returns
            -------
            None : NoneType
                None.
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

    m.def(
        "matrix_inverse",
        [](std::vector<std::vector<real>> mat, const real &tol) {
            std::vector<std::vector<real>> invMat(
                mat.size(), std::vector<real>(mat[0].size()));
            mat_inv(mat, invMat, tol);
            return invMat;
        },
        py::arg("mat"), py::arg("tol") = 1.0e-16L, R"mydelimiter(
        Calculate the inverse of a matrix using LU decomposition.

        Parameters
        ----------
        mat : list of list of real
            Matrix to invert.
        tol : real, optional
            Tolerance for the matrix inversion, by default 1.0e-16L.

        Returns
        -------
        invMat : list of list of real
            Inverse of the matrix.
        )mydelimiter");

    m.def("delta_at_utc", &delta_at_utc, py::arg("mjdUtc"), R"mydelimiter(
        Calculate the difference between TAI and UTC time with UTC as the input.

        Parameters
        ----------
        mjdUtc : real
            Modified Julian date in UTC.

        Returns
        -------
        deltaAt : real
            Difference between TAI and UTC time.
        )mydelimiter");

    m.def("delta_at_tai", &delta_at_tai, py::arg("mjdTai"), R"mydelimiter(
        Calculate the difference between TAI and UTC time with TAI as the input.

        Parameters
        ----------
        mjdTai : real
            Modified Julian date in TAI.

        Returns
        -------
        deltaAt : real
            Difference between TAI and UTC time.
        )mydelimiter");

    m.def("delta_et_utc", &delta_et_utc, py::arg("mjdUtc"), R"mydelimiter(
        Calculate the difference between TDB and UTC time with UTC as the input.

        Parameters
        ----------
        mjdUtc : real
            Modified Julian date in UTC.

        Returns
        -------
        delta_et : real
            Difference between TDB and UTC time.
        )mydelimiter");

    m.def("delta_et_tdb", &delta_et_tdb, py::arg("mjdTdb"), R"mydelimiter(
        Calculate the difference between TDB and UTC time with TDB as the input.

        Parameters
        ----------
        mjdTdb : real
            Modified Julian date in TDB.

        Returns
        -------
        delta_et : real
            Difference between TDB and UTC time.
        )mydelimiter");

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
        .def_readwrite("caTol", &Body::caTol, R"mydelimiter(
            Distance tolerance for close approaches.
            )mydelimiter")
        .def("set_J2", &Body::set_J2, py::arg("J2"), py::arg("poleRA"),
             py::arg("poleDec"), R"mydelimiter(
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
        .def(py::init<std::string, int, real, real, real>(), py::arg("name"),
             py::arg("spiceId"), py::arg("t0"), py::arg("mass"),
             py::arg("radius"),
             R"mydelimiter(
            Constructor for the SpiceBody class.

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
        .def(py::init<std::string, real, real, real, std::vector<real>,
                      NongravParameters>(),
             py::arg("name"), py::arg("t0"), py::arg("mass"), py::arg("radius"),
             py::arg("cometaryState"), py::arg("ngParams"),
             R"mydelimiter(
            Constructor for the IntegBody class.

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
            ngParams : PropSimulation.NongravParameters
                Non-gravitational parameters of the body.
            )mydelimiter")

        .def(py::init<std::string, real, real, real, std::vector<real>,
                      std::vector<real>, NongravParameters>(),
             py::arg("name"), py::arg("t0"), py::arg("mass"), py::arg("radius"),
             py::arg("pos"), py::arg("vel"), py::arg("ngParams"), R"mydelimiter(
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
            ngParams : PropSimulation.NongravParameters
                Non-gravitational parameters of the body.
            )mydelimiter")
        .def_readwrite("spiceId", &IntegBody::spiceId, R"mydelimiter(
            SPICE ID of the body.
            )mydelimiter")
        .def_readwrite("isCometary", &IntegBody::isCometary, R"mydelimiter(
            Whether the body is a cometary body.
            )mydelimiter")
        .def_readwrite("initState", &IntegBody::initState, R"mydelimiter(
            Initial input state of the body (Cometary heliocentric/Cartesian barycentric).
            )mydelimiter")
        .def_readwrite("initCart", &IntegBody::initCart, R"mydelimiter(
            Initial barycentric Cartesian state of the body.
            )mydelimiter")
        .def_readwrite("isInteg", &IntegBody::isInteg, R"mydelimiter(
            Whether the body is an integrated body. Always True.
            )mydelimiter")
        .def_readwrite("isThrusting", &IntegBody::isThrusting, R"mydelimiter(
            Whether the body is thrusting.
            )mydelimiter")
        .def_readwrite("ngParams", &IntegBody::ngParams, R"mydelimiter(
            Non-gravitational parameters of the body.
            )mydelimiter")
        .def_readwrite("n2Derivs", &IntegBody::n2Derivs, R"mydelimiter(
            Number of second derivatives of the body.
            )mydelimiter")
        .def_readwrite("propStm", &IntegBody::propStm, R"mydelimiter(
            Boolean for whether to propagate the state transition matrix of the body.
            )mydelimiter")
        .def_readwrite("stm", &IntegBody::stm, R"mydelimiter(
            State transition matrix of the body.
            )mydelimiter")
        .def_readwrite("dCartdState", &IntegBody::dCartdState, R"mydelimiter(
            Partials of initial cartesian state with repect to initial input state of the body.
            )mydelimiter")
        .def("prepare_stm", &IntegBody::prepare_stm, R"mydelimiter(
            Prepare the state transition matrix of the body for propagation.

            Returns
            -------
            None : NoneType
                None.
            )mydelimiter");

    m.def("reconstruct_stm", &reconstruct_stm, py::arg("stm"), R"mydelimiter(
        Reconstruct the state transition matrix from the flattened vector.

        Parameters
        ----------
        stm : list of real
            Flattened state transition matrix.

        Returns
        -------
        stmMat : list of list of real
            Reconstructed state transition matrix.
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

    m.def("propSim_parallel_omp", &propSim_parallel_omp, py::arg("refSim"),
          py::arg("allBodies"), py::arg("isCometary"), R"mydelimiter(
        Propagate a simulation in parallel using OpenMP.

        Parameters
        ----------
        refSim : PropSimulation
            Reference simulation to copy.
        allBodies : list of list of real
            List of all bodies to propagate. Each list contains the initial MJD TDB time,
            mass, radius, initial state, and list of non-gravitational parameters of the body.
            The initial state is either the initial Heliocentric Ecliptic Cometary state
            or the initial barycentric Cartesian state (position and velocity separated).
        isCometary : bool
            Whether the bodies are cometary bodies.

        Returns
        -------
        allSims : list of PropSimulation
            List of all simulations propagated in parallel.
        )mydelimiter");

    py::class_<PropSimulation>(m, "PropSimulation", R"mydelimiter(
        The PropSimulation class contains the orbit propagation simulation for intgrating solar system small bodies.
        )mydelimiter")
        .def(py::init<std::string, real, const int, std::string>(),
             py::arg("name"), py::arg("t0"), py::arg("defaultSpiceBodies"),
             py::arg("DEkernelPath"), R"mydelimiter(
            Constructor for the PropSimulation class.

            name : str
                Name of the simulation.
            t0 : real
                Initial MJD TDB time of the simulation.
            defaultSpiceBodies : int
                Version of the DE kernel to get the default SPICE bodies from.
            DEkernelPath : str
                Path to the SPICE DE kernel.
            )mydelimiter")
        .def(py::init<std::string, const PropSimulation &>(), py::arg("name"),
             py::arg("simRef"), R"mydelimiter(
            Constructor for the PropSimulation class.

            name : str
                Name of the simulation.
            simRef : PropSimulation
                Simulation to copy.
            )mydelimiter")
        .def_readwrite("name", &PropSimulation::name, R"mydelimiter(
            Name of the simulation.
            )mydelimiter")
        .def_readwrite("DEkernelPath", &PropSimulation::DEkernelPath,
                       R"mydelimiter(
            Path to the SPICE DE kernel.
            )mydelimiter")
        .def_readwrite("consts", &PropSimulation::consts, R"mydelimiter(
            Constants of the simulation. PropSimulation.Constants object.
            )mydelimiter")
        .def_readwrite("integParams", &PropSimulation::integParams,
                       R"mydelimiter(
            Integration parameters of the simulation. PropSimulation.IntegParams object.
            )mydelimiter")
        .def_readwrite("spiceBodies", &PropSimulation::spiceBodies,
                       R"mydelimiter(
            SPICE bodies of the simulation. List of PropSimulation.SpiceBodies objects.
            )mydelimiter")
        .def_readwrite("integBodies", &PropSimulation::integBodies,
                       R"mydelimiter(
            Integration bodies of the simulation. List of PropSimulation.IntegBody objects.
            )mydelimiter")
        .def_readwrite("events", &PropSimulation::events, R"mydelimiter(
            Events of the simulation. List of PropSimulation.Event objects.
            )mydelimiter")
        .def_readwrite("caParams", &PropSimulation::caParams, R"mydelimiter(
            Close approach parameters of the simulation. List of PropSimulation.CloseApproachParameters objects.
            )mydelimiter")
        .def_readwrite("impactParams", &PropSimulation::impactParams,
                       R"mydelimiter(
            Impact parameters of the simulation. List of PropSimulation.ImpactParameters objects.
            )mydelimiter")
        .def_readwrite("t", &PropSimulation::t, R"mydelimiter(
            Current time of the simulation.
            )mydelimiter")
        .def_readwrite("xInteg", &PropSimulation::xInteg, R"mydelimiter(
            Current states of each integration body in the simulation.
            )mydelimiter")
        .def_readwrite("interpParams", &PropSimulation::interpParams,
                       R"mydelimiter(
            Interpolation parameters of the simulation. PropSimulation.InterpolationParameters object.
            )mydelimiter")
        .def_readwrite("tEvalUTC", &PropSimulation::tEvalUTC, R"mydelimiter(
            Whether the MJD evaluation time is in UTC for each value in PropSimulation.tEval,
            as opposed to TDB.
            )mydelimiter")
        .def_readwrite("evalApparentState", &PropSimulation::evalApparentState,
                       R"mydelimiter(
            Whether to evaluate the apparent state of the integration bodies.
            )mydelimiter")
        .def_readwrite("evalMeasurements", &PropSimulation::evalMeasurements,
                       R"mydelimiter(
            Whether to evaluate the measurements of the integration bodies.
            )mydelimiter")
        .def_readwrite("convergedLightTime",
                       &PropSimulation::convergedLightTime, R"mydelimiter(
            Whether to use converged Newtonian light time correction.
            )mydelimiter")
        .def_readwrite("xObserver", &PropSimulation::xObserver, R"mydelimiter(
            State of the observer for each value in PropSimulation.tEval.
            )mydelimiter")
        .def_readwrite("observerInfo", &PropSimulation::observerInfo,
                       R"mydelimiter(
            Observer information for each value in PropSimulation.tEval.
            )mydelimiter")
        .def_readwrite("tEvalMargin", &PropSimulation::tEvalMargin,
                       R"mydelimiter(
            Margin for allowing evaluation past the propagation start and end times.
            )mydelimiter")
        .def_readwrite("tEval", &PropSimulation::tEval, R"mydelimiter(
            MJD Times to evaluate the states of the integrated bodies at.
            Can be TDB or UTC based on PropSimulation.tEvalUTC.
            )mydelimiter")
        .def_readwrite("radarObserver", &PropSimulation::radarObserver,
                       R"mydelimiter(
            Whether the observer for each value in PropSimulation.tEval is for radar.
            )mydelimiter")
        .def_readwrite("lightTimeEval", &PropSimulation::lightTimeEval,
                       R"mydelimiter(
            Light time from the observer to each integration body for each value in PropSimulation.tEval.
            )mydelimiter")
        .def_readwrite("xIntegEval", &PropSimulation::xIntegEval, R"mydelimiter(
            States of each integration body in the simulation for each value in PropSimulation.tEval.
            )mydelimiter")
        .def_readwrite("opticalObs", &PropSimulation::opticalObs,
                       R"mydelimiter(
            Optical observation of each integration body in the simulation for each value in PropSimulation.tEval.
            )mydelimiter")
        .def_readwrite("opticalPartials", &PropSimulation::opticalPartials,
                       R"mydelimiter(
            Optical observation partials of each integration body in the simulation for each value in PropSimulation.tEval.
            )mydelimiter")
        .def_readwrite("radarObs", &PropSimulation::radarObs,
                       R"mydelimiter(
            Radar observation of each integration body in the simulation for each value in PropSimulation.tEval.
            )mydelimiter")
        .def_readwrite("radarPartials", &PropSimulation::radarPartials,
                       R"mydelimiter(
            Radar observation partials of each integration body in the simulation for each value in PropSimulation.tEval.
            )mydelimiter")
        .def("interpolate", &PropSimulation::interpolate, py::arg("t"),
             R"mydelimiter(
            Interpolates the states of the integrated bodies to a given time.

            Parameters
            ----------
            t : real
                Time to interpolate to.

            Returns
            -------
            xIntegInterp : list of real
                Interpolated states of the integration bodies.
            )mydelimiter")
        .def("add_spice_body", &PropSimulation::add_spice_body, py::arg("body"),
             R"mydelimiter(
            Adds a SPICE body to the simulation.

            Parameters
            ----------
            body : PropSimulation.SpiceBody
                SPICE body to add to the simulation.
            )mydelimiter")
        .def("map_ephemeris", &PropSimulation::map_ephemeris, R"mydelimiter(
            Memory maps the ephemeris of the simulation.

            Returns
            -------
            None : NoneType
                None.
            )mydelimiter")
        .def("unmap_ephemeris", &PropSimulation::unmap_ephemeris, R"mydelimiter(
            Unmaps the ephemeris of the simulation.

            Returns
            -------
            None : NoneType
                None.
            )mydelimiter")
        .def("get_spiceBody_state", &PropSimulation::get_spiceBody_state,
             py::arg("t"), py::arg("bodyName"), R"mydelimiter(
            Gets the state of a SPICE body at a given time.

            Parameters
            ----------
            t : real
                Time to get the state at.
            bodyName : str
                Name of the SPICE body in the simulation.
            )mydelimiter")
        .def("add_integ_body", &PropSimulation::add_integ_body, py::arg("body"),
             R"mydelimiter(
            Adds an integration body to the simulation.

            Parameters
            ----------
            body : PropSimulation.IntegBody
                Integration body to add to the simulation.
            )mydelimiter")
        .def("remove_body", &PropSimulation::remove_body, py::arg("name"),
             R"mydelimiter(
            Removes a body from the simulation.

            Parameters
            ----------
            name : str
                Name of the body to remove.
            )mydelimiter")
        .def("add_event", &PropSimulation::add_event, py::arg("body"),
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
        .def("set_sim_constants", &PropSimulation::set_sim_constants,
             py::arg("du2m") = 149597870700.0L, py::arg("tu2s") = 86400.0L,
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
            tu2s : real
                Conversion factor from time units to seconds.
            G : real
                Gravitational constant.
            clight : real
                Speed of light in a vacuum.
            )mydelimiter")
        .def("set_integration_parameters",
             &PropSimulation::set_integration_parameters, py::arg("tf"),
             py::arg("tEval") = std::vector<real>(),
             py::arg("tEvalUTC") = false, py::arg("evalApparentState") = false,
             py::arg("convergedLightTims") = false,
             py::arg("observerInfo") = std::vector<std::vector<real>>(),
             py::arg("adaptiveTimestep") = true, py::arg("dt0") = 0.0L,
             py::arg("dtMax") = 21.0L, py::arg("dtMin") = 5.0e-3L,
             py::arg("dtChangeFactor") = 0.25L, py::arg("tolInteg") = 1.0e-11L,
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
                PropSimulation.tEval, as opposed to TDB.
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
        .def("get_sim_constants", &PropSimulation::get_sim_constants,
             R"mydelimiter(
                Gets the constants of the simulation.

                Returns
                -------
                du2m : real
                    Conversion factor from distance units to meters.
                tu2s : real
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
             &PropSimulation::get_integration_parameters, R"mydelimiter(
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
        .def("integrate", &PropSimulation::integrate, R"mydelimiter(
            Propagates the simulation using the Gauss-Radau integrator.
            )mydelimiter")
        .def("extend", &PropSimulation::extend, py::arg("tf"),
             py::arg("tEvalNew") = std::vector<real>(),
             py::arg("xObserverNew") = std::vector<std::vector<real>>(),
             R"mydelimiter(
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
            )mydelimiter")
        .def("save", &PropSimulation::save, py::arg("filename"),
             R"mydelimiter(
            Saves the simulation to a file.

            Parameters
            ----------
            filename : str
                Name of the file to save the simulation to.
            )mydelimiter");
}
