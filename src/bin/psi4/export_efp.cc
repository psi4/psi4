#include <boost/python.hpp>
#include <libefp_solver/efp_solver.h>
#include <liboptions/liboptions.h>

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
        def("SetGeometry", &EFP::SetGeometry, "docstring");
}

