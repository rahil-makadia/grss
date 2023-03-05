#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#include "utilities.h"
#include "force.h"
#include "simulation.h"
#include "interpolate.h"
#include "gr15.h"
#include "grss.h"

namespace py = pybind11;

PYBIND11_MODULE(cppgrss, m) {
    m.doc() = "pybind11 binding for grss library";

    // from utilities.h
    py::class_<Constants>(m, "Constants")
        .def(py::init<>())
        .def_readwrite("du2m", &Constants::du2m)
        .def_readwrite("tu2sec", &Constants::tu2sec)
        .def_readwrite("G", &Constants::G)
        .def_readwrite("clight", &Constants::clight)
        .def_readwrite("j2000Jd", &Constants::j2000Jd)
        .def_readwrite("JdMinusMjd", &Constants::JdMinusMjd);
    
    py::class_<IntegrationParameters>(m, "IntegrationParameters")
        .def(py::init<>())
        .def_readwrite("nInteg", &IntegrationParameters::nInteg)
        .def_readwrite("nSpice", &IntegrationParameters::nSpice)
        .def_readwrite("nTotal", &IntegrationParameters::nTotal)
        .def_readwrite("t0", &IntegrationParameters::t0)
        .def_readwrite("tf", &IntegrationParameters::tf)
        .def_readwrite("dt0", &IntegrationParameters::dt0)
        .def_readwrite("dtMax", &IntegrationParameters::dtMax)
        .def_readwrite("dtMin", &IntegrationParameters::dtMin)
        .def_readwrite("dtChangeFactor", &IntegrationParameters::dtChangeFactor)
        .def_readwrite("adaptiveTimestep", &IntegrationParameters::adaptiveTimestep)
        .def_readonly("timestepCounter", &IntegrationParameters::timestepCounter)
        .def_readwrite("tolPC", &IntegrationParameters::tolPC)
        .def_readwrite("tolInteg", &IntegrationParameters::tolInteg);

    py::class_<NongravParamaters>(m, "NongravParamaters")
        .def(py::init<>())
        .def_readwrite("a1", &NongravParamaters::a1)
        .def_readwrite("a2", &NongravParamaters::a2)
        .def_readwrite("a3", &NongravParamaters::a3)
        .def_readwrite("alpha", &NongravParamaters::alpha)
        .def_readwrite("k", &NongravParamaters::k)
        .def_readwrite("m", &NongravParamaters::m)
        .def_readwrite("n", &NongravParamaters::n)
        .def_readwrite("r0_au", &NongravParamaters::r0_au);

    py::class_<Event>(m, "Event")
        .def(py::init<>())
        .def_readwrite("t", &ImpulseEvent::t)
        .def_readwrite("bodyName", &ImpulseEvent::bodyName)
        .def_readwrite("bodyIndex", &ImpulseEvent::bodyIndex);
    
    py::class_<ImpulseEvent, Event>(m, "ImpulseEvent")
        .def(py::init<>())
        .def_readwrite("deltaV", &ImpulseEvent::deltaV)
        .def("apply", &ImpulseEvent::apply, py::arg("t"), py::arg("xInteg"), py::arg("propDir"));

    // // from force.h
    py::class_<ForceParameters>(m, "ForceParameters")
        .def(py::init<>())
        .def_readwrite("masses", &ForceParameters::masses)
        .def_readwrite("radii", &ForceParameters::radii)
        .def_readwrite("spiceIdList", &ForceParameters::spiceIdList)
        .def_readwrite("ngParamsList", &ForceParameters::ngParamsList)
        .def_readwrite("isPPNList", &ForceParameters::isPPNList)
        .def_readwrite("isJ2List", &ForceParameters::isJ2List)
        .def_readwrite("J2List", &ForceParameters::J2List)
        .def_readwrite("obliquityList", &ForceParameters::obliquityList)
        .def_readwrite("isNongravList", &ForceParameters::isNongravList);

    // from simulation.h
    py::class_<Body>(m, "Body")
        .def(py::init<>())
        .def_readwrite("t0", &Body::t0)
        .def_readwrite("mass", &Body::mass)
        .def_readwrite("radius", &Body::radius)
        .def_readwrite("J2", &Body::J2)
        .def_readwrite("obliquityToEcliptic", &Body::obliquityToEcliptic)
        .def_readwrite("name", &Body::name)
        .def_readwrite("pos", &Body::pos)
        .def_readwrite("vel", &Body::vel)
        .def_readwrite("isPPN", &Body::isPPN)
        .def_readwrite("isJ2", &Body::isJ2)
        .def_readwrite("isNongrav", &Body::isNongrav)
        .def("set_J2", &Body::set_J2, py::arg("J2"), py::arg("obliquityToEcliptic"));

    py::class_<SpiceBody, Body>(m, "SpiceBody")
        .def(py::init<std::string, std::string, int, real, real, real, Constants>(), py::arg("DEkernelPath"), py::arg("name"), py::arg("spiceId"), py::arg("t0"), py::arg("mass"), py::arg("radius"), py::arg("constants"))
        .def_readwrite("spiceId", &SpiceBody::spiceId)
        .def_readwrite("isSpice", &SpiceBody::isSpice);

    py::class_<IntegBody, Body>(m, "IntegBody")
        .def(py::init<std::string, std::string, real, real, real, std::vector<real>, std::vector< std::vector<real> >, NongravParamaters, Constants>(), py::arg("DEkernelPath"), py::arg("name"), py::arg("t0"), py::arg("mass"), py::arg("radius"), py::arg("cometaryState"), py::arg("covariance"), py::arg("ngParams"), py::arg("constants"))
        .def(py::init<std::string, real, real, real, std::vector<real>, std::vector<real>, std::vector< std::vector<real> >, NongravParamaters, Constants>(), py::arg("name"), py::arg("t0"), py::arg("mass"), py::arg("radius"), py::arg("pos"), py::arg("vel"), py::arg("covariance"), py::arg("ngParams"), py::arg("constants"))
        .def_readwrite("isInteg", &IntegBody::isInteg)
        .def_readwrite("covariance", &IntegBody::covariance)
        .def_readwrite("ngParams", &IntegBody::ngParams);

    py::class_<propSimulation>(m, "propSimulation")
        .def(py::init<std::string, real, const int, std::string>(), py::arg("name"), py::arg("t0"), py::arg("defaultSpiceBodies"), py::arg("DEkernelPath"))
        .def(py::init<std::string, const propSimulation&>(), py::arg("name"), py::arg("simRef"))
        .def_readwrite("name", &propSimulation::name)
        .def_readwrite("DEkernelPath", &propSimulation::DEkernelPath)
        .def_readwrite("consts", &propSimulation::consts)
        .def_readwrite("integParams", &propSimulation::integParams)
        .def_readwrite("spiceBodies", &propSimulation::spiceBodies)
        .def_readwrite("integBodies", &propSimulation::integBodies)
        .def_readwrite("events", &propSimulation::events)
        .def_readwrite("t", &propSimulation::t)
        .def_readwrite("xInteg", &propSimulation::xInteg)
        .def_readwrite("forceParams", &propSimulation::forceParams)
        .def_readwrite("evalApparentState", &propSimulation::evalApparentState)
        .def_readwrite("convergedLightTime", &propSimulation::convergedLightTime)
        .def_readwrite("tEvalUTC", &propSimulation::tEvalUTC)
        .def_readwrite("xObserver", &propSimulation::xObserver)
        .def_readwrite("observerInfo", &propSimulation::observerInfo)
        .def_readwrite("tEvalMargin", &propSimulation::tEvalMargin)
        .def_readwrite("tEval", &propSimulation::tEval)
        .def_readwrite("radarObserver", &propSimulation::radarObserver)
        .def_readwrite("lightTimeEval", &propSimulation::lightTimeEval)
        .def_readwrite("xIntegEval", &propSimulation::xIntegEval)
        .def_readwrite("radarObsEval", &propSimulation::radarObsEval)
        .def("add_spice_body", static_cast<void (propSimulation::*)(std::string, std::string, int, real, real, real, Constants)>(&propSimulation::add_spice_body), py::arg("DEkernelPath"), py::arg("name"), py::arg("spiceId"), py::arg("t0"), py::arg("mass"), py::arg("radius"), py::arg("constants"))
        .def("add_spice_body", static_cast<void (propSimulation::*)(SpiceBody)>(&propSimulation::add_spice_body), py::arg("body"))
        .def("add_integ_body", static_cast<void (propSimulation::*)(std::string, std::string, real, real, real, std::vector<real>, std::vector< std::vector<real> >, NongravParamaters, Constants)>(&propSimulation::add_integ_body), py::arg("DEkernelPath"), py::arg("name"), py::arg("t0"), py::arg("mass"), py::arg("radius"), py::arg("cometaryState"), py::arg("covariance"), py::arg("ngParams"), py::arg("constants"))
        .def("add_integ_body", static_cast<void (propSimulation::*)(std::string, real, real, real, std::vector<real>, std::vector<real>, std::vector< std::vector<real> >, NongravParamaters, Constants)>(&propSimulation::add_integ_body), py::arg("name"), py::arg("t0"), py::arg("mass"), py::arg("radius"), py::arg("pos"), py::arg("vel"), py::arg("covariance"), py::arg("ngParams"), py::arg("constants"))
        .def("add_integ_body", static_cast<void (propSimulation::*)(IntegBody)>(&propSimulation::add_integ_body), py::arg("body"))
        .def("remove_body", &propSimulation::remove_body, py::arg("name"))
        .def("add_event", &propSimulation::add_event, py::arg("body"), py::arg("tEvent"), py::arg("deltaV"))
        .def("set_sim_constants", &propSimulation::set_sim_constants, py::arg("du2m")=149597870700.0L, py::arg("tu2sec")=86400.0L, py::arg("G")=6.6743e-11L/(149597870700.0L*149597870700.0L*149597870700.0L)*86400.0L*86400.0L, py::arg("clight")=299792458.0L/149597870700.0L*86400.0L)
        .def("set_integration_parameters", &propSimulation::set_integration_parameters, py::arg("tf"), py::arg("tEval")=std::vector<real>(), py::arg("tEvalUTC")=false, py::arg("evalApparentState")=false, py::arg("convergedLightTims")=false, py::arg("observerInfo")=std::vector< std::vector<real> >(), py::arg("adaptiveTimestep")=true, py::arg("dt0")=0.0L, py::arg("dtMax")=6.0L, py::arg("dtMin")=7.0e-3L, py::arg("dtChangeFactor")=0.25L, py::arg("tolInteg")=1.0e-6L, py::arg("tolPC")=1.0e-16L)
        .def("get_sim_constants", &propSimulation::get_sim_constants)
        .def("get_integration_parameters", &propSimulation::get_integration_parameters)
        .def("preprocess", &propSimulation::preprocess)
        .def("integrate", &propSimulation::integrate)
        .def("extend", &propSimulation::extend, py::arg("tf"), py::arg("tEvalNew")=std::vector<real>(), py::arg("xObserverNew")=std::vector< std::vector<real> >());

    #ifdef VERSION_INFO
        m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
    #else
        m.attr("__version__") = "dev";
    #endif
}
