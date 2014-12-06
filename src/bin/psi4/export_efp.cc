#include <boost/python.hpp>
#include <libefp_solver/efp_solver.h>
#include <liboptions/liboptions.h>

using namespace boost::python;
using namespace psi;
using namespace psi::efp;

void export_efp()
{
    // because there is no default constructor for libefp, need flag
    // "no_init" and the constructor definition, def(init<Options&>())
    class_<EFP, boost::shared_ptr<EFP> >("EFP", "Class interfacing with libefp", no_init).
        def(init<Options&>()).
        def("compute", &EFP::compute, "Computes libefp energies and, if active, torque").
        def("set_qm_atoms", &EFP::set_qm_atoms, "Provides libefp with QM fragment information").
        def("nfragments", &EFP::get_frag_count, "Returns the number of EFP fragments in the molecule").
        def("print_out", &EFP::print_out, "Prints options settings and EFP and QM geometries");
}

