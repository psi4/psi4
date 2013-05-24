#include <boost/python.hpp>
#include <libefp_solver/efp_solver.h>
#include <liboptions/liboptions.h>
#include <libmints/vector.h>
#include <libmints/coordentry.h>

using namespace boost::python;
using namespace psi;
using namespace psi::efp;

void export_efp()
{
    // because there is no default constructor for libefp,
    // need the "no_init" flag and the definition of the constructor:
    // def(init<Options&>())
    class_<EFP, boost::shared_ptr<EFP> >("EFP", "docstring", no_init).
        def(init<Options&>()).
        def("Compute", &EFP::Compute, "docstring").
        def("get_electrostatic_gradient", &EFP::get_electrostatic_gradient, "docstring").
        def("set_qm_atoms", &EFP::set_qm_atoms, "docstring").
        def("print_out", &EFP::print_out, "docstring");
        //def("SetGeometry", &EFP::SetGeometry, "docstring");
}

