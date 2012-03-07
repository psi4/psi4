#include <boost/python.hpp>
#include <libfunctional/superfunctional.h>

using namespace boost;
using namespace boost::python;
using namespace psi::functional;

void export_functional()
{
    class_<SuperFunctional, boost::shared_ptr<SuperFunctional> >("SuperFunctional", "docstring").
        def("create_superfunctional", &SuperFunctional::createSuperFunctional, "docstring").
        staticmethod("create_superfunctional").
        def("build_superfunctional", &SuperFunctional::buildSuperFunctional, "docstring").
        staticmethod("build_superfunctional").
        def("available_superfunctionals", &SuperFunctional::availableSuperFunctionals, "docstring").
        staticmethod("available_superfunctionals").
        def("available_names", &SuperFunctional::availableNames, "docstring").
        staticmethod("available_names").
        //def("test_superfunctionals", &SuperFunctional::testSuperFunctionals).
        //staticmethod("test_superfunctionals").
        //def("test_superfunctional", &SuperFunctional::testSuperFunctional).

        def("get_name", &SuperFunctional::getName, "docstring").
        def("get_description", &SuperFunctional::getDescription, "docstring").
        def("get_citation", &SuperFunctional::getCitation, "docstring").
        def("get_composition", &SuperFunctional::getComposition, "docstring").
        def("get_exact_exchange", &SuperFunctional::getExactExchange, "docstring").
        def("get_pt2", &SuperFunctional::getPT2, "docstring").
        def("get_omega", &SuperFunctional::getOmega, "docstring").
        def("get_dash_d", &SuperFunctional::getDashD, "docstring").
        def("get_npoints", &SuperFunctional::getNPoints, "docstring").
        def("get_deriv", &SuperFunctional::getDeriv, "docstring").
        def("get_functional", &SuperFunctional::getFunctional, "docstring").
        def("get_weight", &SuperFunctional::getWeight, "docstring").
        def("get_size", &SuperFunctional::size, "docstring").

        def("set_name", &SuperFunctional::setName, "docstring").
        def("set_description", &SuperFunctional::setDescription, "docstring").
        def("set_citation", &SuperFunctional::setCitation, "docstring").
        def("set_parameter", &SuperFunctional::setParameter, "docstring").
        def("set_exact_exchange", &SuperFunctional::setExactExchange, "docstring").
        def("set_pt2", &SuperFunctional::setPT2, "docstring").
        def("set_omega", &SuperFunctional::setOmega, "docstring").
        def("set_dash_d", &SuperFunctional::setDashD, "docstring").
        def("set_npoints", &SuperFunctional::setNPoints, "docstring").
        def("set_deriv", &SuperFunctional::setDeriv, "docstring").
        def("set_size", &SuperFunctional::size, "docstring").

        def("is_gga", &SuperFunctional::isGGA, "docstring").
        def("is_meta", &SuperFunctional::isMeta, "docstring").
        def("is_hybrid", &SuperFunctional::isHybrid, "docstring").
        def("is_double_hybrid", &SuperFunctional::isDoubleHybrid, "docstring").
        def("is_range_corrected", &SuperFunctional::isRangeCorrected, "docstring").
        def("is_dash_d", &SuperFunctional::isDashD, "docstring").

        //def("add_functional", &SuperFunctional::addFunctional).

        def("computeRKSFunctional", &SuperFunctional::computeRKSFunctional, "docstring").
        def("computeUKSFunctional", &SuperFunctional::computeUKSFunctional, "docstring");

    class_<Functional, boost::shared_ptr<Functional>, boost::noncopyable >("Functional", "docstring", no_init).
        def("create_functional", &Functional::createFunctional, "docstring").
        staticmethod("create_functional").
        def("available_functionals", &Functional::availableFunctionals, "docstring").
        staticmethod("available_functionals").
        def("available_names", &Functional::availableNames, "docstring").
        staticmethod("available_names").
        //def("test_functionals", &Functional::testFunctionals).
        //staticmethod("test_functionals").
        //def("test_functional", &Functional::testFunctional).

        def("get_name", &Functional::getName, "docstring").
        def("get_description", &Functional::getDescription, "docstring").
        def("get_citation", &Functional::getCitation, "docstring").
        def("get_parameters_string", &Functional::getParametersString, "docstring").
        def("get_parameters", &Functional::getParameters, "docstring").
        def("get_npoints", &Functional::getNPoints, "docstring").
        def("get_deriv", &Functional::getDeriv, "docstring").
        def("get_density_cutoff", &Functional::getDensityCutoff, "docstring").

        def("set_name", &Functional::setName, "docstring").
        def("set_description", &Functional::setDescription, "docstring").
        def("set_citation", &Functional::setCitation, "docstring").
        def("set_parameter", &Functional::setParameter, "docstring").
        def("set_parameters", &Functional::setParameters, "docstring").
        def("set_npoints", &Functional::setNPoints, "docstring").
        def("set_deriv", &Functional::setDeriv, "docstring").
        def("set_density_cutoff", &Functional::setDensityCutoff, "docstring").

        def("is_gga", &Functional::isGGA, "docstring").
        def("is_meta", &Functional::isMeta, "docstring").

        def("computeRKSFunctional", &Functional::computeRKSFunctional, "docstring").
        def("computeUKSFunctional", &Functional::computeUKSFunctional, "docstring");
}
