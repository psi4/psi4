#include <boost/python.hpp>
#include <libfunctional/superfunctional.h>
#include <libfunctional/functional.h>
#include <libdisp/dispersion.h>

using namespace boost;
using namespace boost::python;
using namespace psi;

void export_functional()
{
    class_<SuperFunctional, boost::shared_ptr<SuperFunctional>, boost::noncopyable >("SuperFunctional", "docstring", no_init).
        def("build", &SuperFunctional::build, "docstring").
        staticmethod("build").
        def("name", &SuperFunctional::name, "docstring").
        def("description", &SuperFunctional::description, "docstring").
        def("citation", &SuperFunctional::citation, "docstring").
        def("print_out",&SuperFunctional::py_print, "docstring");

    class_<Functional, boost::shared_ptr<Functional>, boost::noncopyable >("Functional", "docstring", no_init).
        def("build", &Functional::build, "docstring").
        staticmethod("build").
        def("build_base", &Functional::build_base, "docstring").
        staticmethod("build_base").
        def("name", &Functional::name, "docstring").
        def("description", &Functional::description, "docstring").
        def("citation", &Functional::citation, "docstring").
        def("print_out", &Functional::py_print, "docstring");

    class_<Dispersion, boost::shared_ptr<Dispersion>, boost::noncopyable >("Dispersion", "docstring", no_init).
        def("build", &Dispersion::build, "docstring").
        staticmethod("build").
        def("name", &Dispersion::name, "docstring").
        def("description", &Dispersion::description, "docstring").
        def("citation", &Dispersion::citation, "docstring").
        def("set_name", &Dispersion::set_name, "docstring").
        def("set_description", &Dispersion::set_description, "docstring").
        def("set_citation", &Dispersion::set_citation, "docstring").
        def("print_energy", &Dispersion::print_energy, "docstring").
        def("print_gradient", &Dispersion::print_gradient, "docstring").
        def("print_hessian", &Dispersion::print_hessian, "docstring").
        //def("compute_energy", &Dispersion::compute_energy, "docstring").
        //def("compute_gradient", &Dispersion::compute_gradient, "docstring").
        //def("compute_hessian", &Dispersion::compute_hessian, "docstring").
        def("print_out",&Dispersion::py_print, "docstring");

}
