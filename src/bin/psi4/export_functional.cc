#include <boost/python.hpp>
#include <libfunctional/superfunctional.h>

using namespace boost;
using namespace boost::python;
using namespace psi::functional;

void export_functional()
{
    class_<SuperFunctional, boost::shared_ptr<SuperFunctional> >("SuperFunctional").
        def("create_superfunctional", &SuperFunctional::createSuperFunctional).
        staticmethod("create_superfunctional").
        def("build_superfunctional", &SuperFunctional::buildSuperFunctional).
        staticmethod("build_superfunctional").
        def("available_superfunctionals", &SuperFunctional::availableSuperFunctionals).
        staticmethod("available_superfunctionals").
        def("available_names", &SuperFunctional::availableNames).
        staticmethod("available_names").
        //def("test_superfunctionals", &SuperFunctional::testSuperFunctionals).
        //staticmethod("test_superfunctionals").
        //def("test_superfunctional", &SuperFunctional::testSuperFunctional).

        def("get_name", &SuperFunctional::getName).
        def("get_description", &SuperFunctional::getDescription).
        def("get_citation", &SuperFunctional::getCitation).
        def("get_composition", &SuperFunctional::getComposition).
        def("get_exact_exchange", &SuperFunctional::getExactExchange).
        def("get_pt2", &SuperFunctional::getPT2).
        def("get_omega", &SuperFunctional::getOmega).
        def("get_dash_d", &SuperFunctional::getDashD).
        def("get_npoints", &SuperFunctional::getNPoints).
        def("get_deriv", &SuperFunctional::getDeriv).
        def("get_functional", &SuperFunctional::getFunctional).
        def("get_weight", &SuperFunctional::getWeight).
        def("get_size", &SuperFunctional::size).

        def("set_name", &SuperFunctional::setName).
        def("set_description", &SuperFunctional::setDescription).
        def("set_citation", &SuperFunctional::setCitation).
        def("set_parameter", &SuperFunctional::setParameter).
        def("set_exact_exchange", &SuperFunctional::setExactExchange).
        def("set_pt2", &SuperFunctional::setPT2).
        def("set_omega", &SuperFunctional::setOmega).
        def("set_dash_d", &SuperFunctional::setDashD).
        def("set_npoints", &SuperFunctional::setNPoints).
        def("set_deriv", &SuperFunctional::setDeriv).
        def("set_size", &SuperFunctional::size).

        def("is_gga", &SuperFunctional::isGGA).
        def("is_meta", &SuperFunctional::isMeta).
        def("is_hybrid", &SuperFunctional::isHybrid).
        def("is_double_hybrid", &SuperFunctional::isDoubleHybrid).
        def("is_range_corrected", &SuperFunctional::isRangeCorrected).
        def("is_dash_d", &SuperFunctional::isDashD).

        //def("add_functional", &SuperFunctional::addFunctional).

        def("computeRKSFunctional", &SuperFunctional::computeRKSFunctional).
        def("computeUKSFunctional", &SuperFunctional::computeUKSFunctional);

    class_<Functional, boost::shared_ptr<Functional>, boost::noncopyable >("Functional", no_init).
        def("create_functional", &Functional::createFunctional).
        staticmethod("create_functional").
        def("available_functionals", &Functional::availableFunctionals).
        staticmethod("available_functionals").
        def("available_names", &Functional::availableNames).
        staticmethod("available_names").
        //def("test_functionals", &Functional::testFunctionals).
        //staticmethod("test_functionals").
        //def("test_functional", &Functional::testFunctional).

        def("get_name", &Functional::getName).
        def("get_description", &Functional::getDescription).
        def("get_citation", &Functional::getCitation).
        def("get_parameters_string", &Functional::getParametersString).
        def("get_parameters", &Functional::getParameters).
        def("get_npoints", &Functional::getNPoints).
        def("get_deriv", &Functional::getDeriv).
        def("get_density_cutoff", &Functional::getDensityCutoff).

        def("set_name", &Functional::setName).
        def("set_description", &Functional::setDescription).
        def("set_citation", &Functional::setCitation).
        def("set_parameter", &Functional::setParameter).
        def("set_parameters", &Functional::setParameters).
        def("set_npoints", &Functional::setNPoints).
        def("set_deriv", &Functional::setDeriv).
        def("set_density_cutoff", &Functional::setDensityCutoff).

        def("is_gga", &Functional::isGGA).
        def("is_meta", &Functional::isMeta).

        def("computeRKSFunctional", &Functional::computeRKSFunctional).
        def("computeUKSFunctional", &Functional::computeUKSFunctional);
}
